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


std::vector<double> normalizeC(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> c){
    std::vector<double> out;

    double max= 0;
    for(auto &v:c){
        if(abs(v.real()) > max){
            max = v.real();
        }
    }
    for(auto &v:c){
        out.push_back((v.real()+max)/(2.*max));
    }
    std::cout<<"result range: ["<<-max<<", "<<max<<"]\n";
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


//    Eigen::nbThreads();
//    Eigen::initParallel();


    double liczbafalowa = 51.9988;;
    fem::TriangleElement::k_ = liczbafalowa;

    const int nVerWidth = 721; //vertices number
    const int nVerHeight = 641; //vertices number


    double width = 4.2;
    double height = 3.2;
    //vertices number

    double widthEleLen = width / (nVerWidth - 1);
    double heightEleLen = height / (nVerHeight - 1);

    std::cout<<"Vertices CREATION BEGIN"<<std::endl;
    std::vector<Vertex2D> points;


    for(int k =0; k<(nVerWidth)*(nVerHeight);k++){
        int i = k%(nVerWidth);
        int j = k/(nVerWidth);

        points.emplace_back(- width/2 + (widthEleLen * i), - height/2 + (heightEleLen * j ), isBorder(k, nVerWidth, nVerHeight));
    }

    std::vector<Wall> walls{};


    //MAKING WALLS AND CONSTRAINTS
    //some tomfoolery
    double thicness = .1;

    walls.emplace_back(-width/2., -height/2., thicness, height, fem::CONCRETE);
    walls.emplace_back(-width/2., -height/2., width, thicness, fem::CONCRETE);
    walls.emplace_back(-width/2., height/2. - thicness, width,  thicness, fem::CONCRETE);
    walls.emplace_back(width/2. - thicness,  - height/2, thicness, height, fem::CONCRETE);
    walls.emplace_back(width/4.,  0, thicness, height, fem::CONCRETE);




    std::cout<<"Vertices CREATION END"<<std::endl<<std::endl;


    std::vector<fem::ElementIndices> elementsIdx;

    // with points rotation implemented inside
    build(points, elementsIdx, walls);


//    // quad elements intersection
//    addQuadVertElements(points, elementsIdx);



    const int N = static_cast<int>(points.size());

    std::vector<fem::TriangleElement> Elements;

    std::cout<<"Elements CREATION BEGIN"<<std::endl;

    for(auto &idx:elementsIdx){
        Elements.emplace_back(0, points, idx.indices, idx.n);
    }


    std::cout<<"Elements CREATION END"<<std::endl<<std::endl;


    std::cout<<"Tiplets CREATION BEGIN"<<std::endl;

    std::vector<Eigen::Triplet<std::complex<double>>> tripletListS{};

    for(auto &ele:Elements){
        std::vector<std::vector<std::complex<double>>> &E = ele.E_;
        std::vector<int> &glob = ele.globalVectorIdx;
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
        std::vector<int> &glob = ele.globalVectorIdx;
        for (int k = 0; k < glob.size(); k++){
                int i = glob[k];
                if(points[i].isBorder){
                    F[i] = 0;
                }else{
                    F[i] += Fm[k];
                }
        }
    }


//    F[(nVerWidth*nVerHeight)/2] = 1; //F[(row*row)/2 - row] = 1; F[(row*row)/2 + row] = 1;
//    F[(row*row)/2 + 1] = 1; F[(row*row)/2 - row + 1] = 1; F[(row*row)/2 + row + 1] = 1;
//    F[(row*row)/2 - 1] = 1; F[(row*row)/2 - row- 1] = 1; F[(row*row)/2 + row- 1] = 1;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> EF(N);
    for(int i = 0; i < N; i++) EF[i] = F[i];



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

//        std::ofstream fileS("resultS.txt");
//        std::ofstream fileF("resultF.txt");
//
//        for(int i = 0; i< N; i++){
//            fileF << F[i] << "\n";
////            for(int j = 0; j<N; j++)
////                fileS << St[i][j] << "\t";
////            fileS << "\n";
//        }
//        fileS.close();
//        fileF.close();

//        std::ofstream outputFileEleMain("elementsVerticesIdxes.txt");
        std::ofstream outputFileEleAir("elementsAirVerticesIdxes.txt");
        std::ofstream outputFileEleWall("elementsWallVerticesIdxes.txt");
        for(auto &e:elementsIdx){
            for(auto &i:e.indices){
                if(e.n == fem::AIR){
                    outputFileEleAir<<i<<" ";
                }else{
                    outputFileEleWall<<i<<" ";
                }

//                outputFileEleMain<<i<<" ";
            }
            if(e.n == fem::AIR){
                outputFileEleAir<<"\n";
            }else{
                outputFileEleWall<<"\n";
            }

//            outputFileEleMain<<"\n";
        }
//        outputFileEleMain.close();
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

    std::cout<<"Solving BEGIN"<<std::endl;
    c = solver.solve(EF);
    std::cout<<"Solving END"<<std::endl<<std::endl;


//if(debug){
//    std::ofstream outputFileC("c.txt");
//    for (int i = 0; i<N; i ++){
//        outputFileC<< c[i].real() <<std::endl;
//    }
//    outputFileC.close();
//}




    std::cout<<"normalizing results vip..."<<std::endl;
    //TODO here should be some results handling
    auto normc = normalizeC(c); //placeholder
    std::cout<<"normalizing results END"<<std::endl;

    std::cout<<"OpenGL BEGIN"<<std::endl;
    CreateOpenGlWindow(points, Elements, normc, width, height);
    std::cout<<"OpenGL END"<<std::endl;


    return 0;
}
