#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>


#include "Elements/TriangleElement.hpp"
#include "Elements/Section.hpp"
#include "Elements/Wall.hpp"
#include "FortunesAlgo/Fortunes.hpp"

#include "Visualization/CreateOpenGlWindow.hpp"


//#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseLU>


void normalizeLogarithmically(std::vector<double>& vec) {
    // Find the maximum element in the vector
    double maxElement = *std::max_element(vec.begin(), vec.end());
    std::cout<<"log range: ["<<0<<", "<<maxElement<<"]\n";
    // Logarithmically normalize each element
    for (double& element : vec) {
        if (maxElement != 0.0) {
            // Avoid log(0) by checking if maxElement is not zero
            element = log(element) / log(maxElement);
        } else {
            // Handle the case when maxElement is zero (avoid division by zero)
            element = 0.0;
        }
    }
}

std::vector<double> normalizeC(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &c){
    std::vector<double> out;

    double maxElement= 0;
    for(auto &v:c){
        if(abs(v) > maxElement){
            maxElement = abs(v);
        }
    }
    for(auto &v:c){
        out.push_back(abs(v)/maxElement);
    }


    std::cout<<"abs result range: ["<<0<<", "<<maxElement<<"]\n";
//    normalizeLogarithmically(out);
    return out;
}


bool isBorder(int k, int w, int h){
    for(int i = 0; i<w; i++){
        if(i==0 || i==w-1){
            for(int j=0;j<h;j++)
                if(k == j*w + i) return true;

        }else
        {
            if(k == i || k == w*(h - 1) + i) return true;
        }
    }
    return false;
}



int main() {


    double liczbafalowa = 51.9988;;
    fem::TriangleElement::k_ = liczbafalowa;

    fem::baseFuncType fType = fem::QUAD;

    const int nVerWidth = 181; //vertices number
    const int nVerHeight = 181; //vertices number


    double width = 5.;
    double height = 5.;
    //vertices number

    double widthEleLen = width / (nVerWidth - 1);
    double heightEleLen = height / (nVerHeight - 1);

    std::cout<<"Vertices CREATION BEGIN"<<std::endl;
    std::vector<Vertex2D> points;


    const int nDampLayers = 4;

    for(int k =0; k<(nVerWidth + 2*nDampLayers)*(nVerHeight + 2*nDampLayers);k++){
        int i = k%(nVerWidth + 2*nDampLayers);
        int j = k/(nVerWidth + 2*nDampLayers);

        points.emplace_back(- width/2 - (nDampLayers*widthEleLen) + (widthEleLen * i), - height/2  - (nDampLayers*heightEleLen)+ (heightEleLen * j ), isBorder(k, nVerWidth + 2*nDampLayers, nVerHeight + 2*nDampLayers));
    }




    //MAKING WALLS AND CONSTRAINTS
    //some tomfoolery
    std::vector<Wall> walls{};
    //
    walls.emplace_back(- width/2 - (nDampLayers*widthEleLen) , - height/2  - (nDampLayers*heightEleLen),
                      nDampLayers*widthEleLen, height + 2.*nDampLayers*heightEleLen,
                       fem::DAMP);

    walls.emplace_back(- width/2 - (nDampLayers*widthEleLen) , - height/2  - (nDampLayers*heightEleLen),
                       width + 2.*nDampLayers*widthEleLen, nDampLayers*heightEleLen,
                       fem::DAMP);

    walls.emplace_back(- width/2 - (nDampLayers*widthEleLen), height/2.,
                       width + 2.*nDampLayers*widthEleLen, nDampLayers*heightEleLen,
                       fem::DAMP);

    walls.emplace_back(width/2.,  -height/2 - nDampLayers*heightEleLen,
                       nDampLayers*widthEleLen, height + 2.*nDampLayers*heightEleLen,
                       fem::DAMP);


//    walls.emplace_back(0,  0, widthEleLen*3, height/2, fem::BRICK);




    std::cout<<"Vertices CREATION END"<<std::endl<<std::endl;


    std::vector<fem::ElementIndices> elementsIdx;

    // with points rotation implemented inside
    build(points, elementsIdx, walls, fType);

    // quad elements intersection
    if(fType == fem::QUAD){
        addQuadVertElements(points, elementsIdx);
    }


    const int N = static_cast<int>(points.size());

    std::vector<fem::TriangleElement> Elements;

    std::cout<<"Elements CREATION BEGIN"<<std::endl;

    for(auto &idx:elementsIdx){
        Elements.emplace_back(0, points, idx);
    }


    std::cout<<"Elements CREATION END"<<std::endl<<std::endl;


    std::cout<<"Tiplets CREATION BEGIN"<<std::endl;

    std::vector<Eigen::Triplet<std::complex<double>>> tripletListS{};

    for(auto &ele:Elements){
        std::vector<std::vector<std::complex<double>>> &E = ele.E_;
        std::vector<int> &glob = ele.globalVectorIdx.indices;
        for (int k=0 ; k < glob.size();k++){
            for(int l = 0; l < glob.size(); l++ ){
                int i = glob[k];
                int j = glob[l];
                if(points[i].isBorder){
//                    std::cout<<i<<" "<<j<<" "<<m<<std::endl;
                }
                else{
                    tripletListS.emplace_back(i, j, E[k][l]);
                }

            }
        }
    }

    for(int i = 0; i < points.size(); i++){
        if(points[i].isBorder){
            tripletListS.emplace_back(i, i, 1.);
        }
    }
    std::cout<<"Triplets CREATION END"<<std::endl<<std::endl;


    std::vector<double> F(N);
    for(auto &ele:Elements){
        std::vector<double> &Fm = ele.F_;
        std::vector<int> &glob = ele.globalVectorIdx.indices;
        for (int k = 0; k < glob.size(); k++){
                int i = glob[k];
                if(points[i].isBorder){
                    F[i] = 0;
                }else{
                    F[i] = Fm[k];
                }
        }
    }

//    F[(nVerWidth*nVerHeight)/2] = 1; //F[(row*row)/2 - row] = 1; F[(row*row)/2 + row] = 1;
//    F[(row*row)/2 + 1] = 1; F[(row*row)/2 - row + 1] = 1; F[(row*row)/2 + row + 1] = 1;
//    F[(row*row)/2 - 1] = 1; F[(row*row)/2 - row- 1] = 1; F[(row*row)/2 + row- 1] = 1;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> EF(N);


    int minIdx = 0;
    double minDist = sqrt(pow(points[minIdx].x - (-1.5),2) +  pow(points[minIdx].y -(-.5), 2));
    for(int i = 1; i < N; i++) {
        double temp = sqrt(pow(points[i].x - (-.5),2) +  pow(points[i].y -(-.5), 2));
        if(temp < minDist){
            minDist = temp;
            minIdx = i;
        }
        EF[i]  = 0;

        //    EF[i] = F[i];
    }
    std::cout<< points[minIdx].x << " " << points[minIdx].x<<std::endl;
    EF[minIdx]  = 1;



    //some debug for saving values to txt files
    bool debug = false;
    if(debug){
        std::ofstream outputFilep("points.txt");
        for(auto & vertex : points){
            outputFilep <<vertex.x << ' ' << vertex.y << ' ' << vertex.isBorder <<std::endl;
        }
        outputFilep.close();

//        std::vector< std::vector<std::complex<double>> > St (N, std::vector<std::complex<double>>(N));

//        for(auto t:tripletListS){
//            St[t.row()][t.col()] += t.value();
//        }

        double maxElement = *std::max_element(F.begin(), F.end());
        std::cout<<maxElement<<std::endl;
//        std::ofstream fileS("resultS.txt");
//        std::ofstream fileF("resultF.txt");
//        for(int i = 0; i< N; i++){
//            if(F[i]> 1e-3)
//            fileF << points[i].x<< " " << points[i].y << " " <<  F[i] <<"\n";
//            for(int j = 0; j<N; j++)
//                fileS << St[i][j] << "\t";
//            fileS << "\n";
//        }
//        fileS.close();
//        fileF.close();

        std::ofstream outputFileEleAir("elementsAirVerticesIdxes.txt");
        std::ofstream outputFileEleWall("elementsWallVerticesIdxes.txt");
        for(auto &e:elementsIdx){
            for(auto &i:e.indices){
                if(e.et == fem::AIR){
                    outputFileEleAir<<i<<" ";
                }else{
                    outputFileEleWall<<i<<" ";
                }
            }
            if(e.et == fem::AIR){
                outputFileEleAir<<"\n";
            }else{
                outputFileEleWall<<"\n";
            }

        }
        outputFileEleAir.close();
        outputFileEleWall.close();

    }



    std::cout<<"SparseMatrix CREATION BEGIN"<<std::endl;
    Eigen::SparseMatrix<std::complex<double>>  S(N, N);
    S.setFromTriplets(tripletListS.begin(), tripletListS.end());
//    tripletListS.clear();
    std::cout<<"SparseMatrix CREATION END"<<std::endl<<std::endl;

    std::cout<<"Solver CREATION BEGIN"<<std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>  solver;

    solver.compute(S);
    ////equivalent to compute above?
    //    solver.analyzePattern(S);
    //    solver.factorize(S);

    std::cout<<"Solver CREATION END"<<std::endl<<std::endl;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> c;


    std::cout<< sizeof(solver)<<std::endl;
    std::cout<< sizeof(S)<<std::endl;


    std::cout<<"Solving BEGIN"<<std::endl;
    c = solver.solve(EF);
    std::cout<<"Solving END"<<std::endl<<std::endl;


if(debug){
    std::ofstream outputFileC("c.txt");
    for (int i = 0; i<N; i ++){
        outputFileC<< abs(c[i]) <<std::endl;
    }
    outputFileC.close();
}




    std::cout<<"normalizing results vip..."<<std::endl;
    //TODO here should be some results handling
    auto normc = normalizeC(c); //placeholder
    std::cout<<"normalizing results END"<<std::endl;

    std::cout<<"OpenGL BEGIN"<<std::endl;
    CreateOpenGlWindow(points, Elements, normc, width, height);
    std::cout<<"OpenGL END"<<std::endl;


    return 0;
}
