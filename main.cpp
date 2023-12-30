#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>


#include "Elements/TriangleElement.hpp"
#include "Elements/Section.hpp"
#include "Elements/Wall.hpp"
#include "FortunesAlgo/Fortunes.hpp"


#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/SparseLU>



bool isBorder(int k, int row){
    for(int i = 0; i<row; i++){
        if(i==0 || i==row-1){
            for(int j=0;j<row;j++)
                if(k == j*row + i) return true;

        }else
        {
            if(k == i || k == row*(row - 1) + i) return true;
        }
    }

    return false;
}



int main() {

    double liczbafalowa = 51.9988;
    const int row = 201; //vertices number
    double L = 3.2;
    //vertices number


    double horLen = L / (row - 1);

    std::cout<<"Vertices CREATION BEGIN"<<std::endl;
    std::vector<Vertex2D> points;


    for(int k =0; k<(row)*(row);k++){
        int i = k/(row);
        int j = k%(row);

        points.emplace_back(- L/2 + (horLen * i), - L/2 + (horLen * j ), isBorder(k,row));
    }

    std::vector<Wall> walls{};


    //some tomfoolery
    double thicness = .1;
    Vertex2D p(-L/2., -L/2.);
    walls.emplace_back(p, thicness, L);
    walls.emplace_back(p, L, thicness);
    p.y = L/2. - thicness;
    walls.emplace_back(p, L,  thicness);
    p.x = L/2. - thicness; p.y =  - L/2;
    walls.emplace_back(p, thicness, L);

    for(auto &w:walls) w.rot(1);


    std::ofstream outputFilep("points.txt");
    for(auto & vertex : points){
        outputFilep <<vertex.x << ' ' << vertex.y << std::endl;
    }
    outputFilep.close();
    std::cout<<"Vertices CREATION END"<<std::endl<<std::endl;


    std::vector<std::vector<int>> elementsIdxWall;
    std::vector<std::vector<int>> elementsIdxAir;

    // with points rotation implemented inside
    build(points, elementsIdxWall, elementsIdxAir, walls);


//    // quad elements intersection
//    addQuadVertElements(points, elementsIdx);


    std::ofstream outputFileEleMain("elementsVerticesIdxes.txt");
    std::ofstream outputFileEle("elementsWallVerticesIdxes.txt");
    for(auto &e:elementsIdxWall){
        for(auto &i:e){
            outputFileEle<<i<<" ";
            outputFileEleMain<<i<<" ";
        }
        outputFileEle<<"\n";
        outputFileEleMain<<"\n";
    }
    outputFileEle.close();

    outputFileEle.open("elementsAirVerticesIdxes.txt");
    for(auto &e:elementsIdxAir){
        for(auto &i:e){
            outputFileEle<<i<<" ";
            outputFileEleMain<<i<<" ";
        }
        outputFileEle<<"\n";
        outputFileEleMain<<"\n";
    }
    outputFileEle.close();

    const int N = static_cast<int>(points.size());

    std::vector<fem::TriangleElement> Elements;

    std::cout<<"Elements CREATION BEGIN"<<std::endl;

    for(auto &idxs:elementsIdxAir){
        Elements.emplace_back(0, points, idxs, liczbafalowa, fem::AIR);
    }

    for(auto &idxs:elementsIdxWall){
        Elements.emplace_back(0, points, idxs, liczbafalowa, fem::CONCRETE);
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

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> EF(N);
    for(int i = 0; i < N; i++) EF[i] = F[i];

    F[(row*row)/2] = 1;// F[(row*row)/2 - row] = 1; F[(row*row)/2 + row] = 1;
//    F[(row*row)/2 + 1] = 1; F[(row*row)/2 - row + 1] = 1; F[(row*row)/2 + row + 1] = 1;
//    F[(row*row)/2 - 1] = 1; F[(row*row)/2 - row- 1] = 1; F[(row*row)/2 + row- 1] = 1;

    bool debug = false;
    if(debug){
//        std::vector< std::vector<std::complex<double>> > St (N, std::vector<std::complex<double>>(N));

//        for(auto t:tripletListS){
//            St[t.row()][t.col()] += t.value();
//        }

        std::ofstream fileS("resultS.txt");
        std::ofstream fileF("resultF.txt");

        for(int i = 0; i< N; i++){
            fileF << F[i] << "\n";
//            for(int j = 0; j<N; j++)
//                fileS << St[i][j] << "\t";
//            fileS << "\n";
        }
        fileS.close();
        fileF.close();
    }



    std::cout<<"SparseMatrix CREATION BEGIN"<<std::endl;
    Eigen::SparseMatrix<std::complex<double>>  S(N, N);
    S.setFromTriplets(tripletListS.begin(), tripletListS.end());
//    tripletListS.clear();
    std::cout<<"SparseMatrix CREATION END"<<std::endl<<std::endl;

    std::cout<<"Solver CREATION BEGIN"<<std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>  solver;

    solver.analyzePattern(S);

    solver.factorize(S);

    std::cout<<"Solver CREATION END"<<std::endl<<std::endl;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> c;

    std::cout<<"Solving BEGIN"<<std::endl;
    c = solver.solve(EF);
    std::cout<<"Solving END"<<std::endl<<std::endl;

    std::ofstream outputFileC("c.txt");
    for (int i = 0; i<N; i ++){
        outputFileC<< c[i].real() <<std::endl;
    }
    outputFileC.close();

    return 0;
}
