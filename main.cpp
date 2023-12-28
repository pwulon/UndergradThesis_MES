#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdio>
#include <set>
#include <unordered_set>
#include <functional>
#include <initializer_list>

#include "Elements/TriangleElement.hpp"
#include "Elements/Vertex2D.hpp"
#include "Elements/Wall.hpp"
#include "FortunesAlgo/Fortunes.hpp"


#include <eigen3/Eigen/SparseCholesky>
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

void pointsRot(std::vector<Point2D> &points, double ang){
    double alfa = ang*M_PI/180.;
    double sin1 = sin(alfa);
    double cos1 = cos(alfa);

    for (int i = 0; i < points.size(); i++) {
        double x = points[i].x*cos1 - points[i].y*sin1;
        double y = points[i].x*sin1 + points[i].y*cos1;
        points[i].x = x;
        points[i].y = y;
    }

}

struct segment{
    int i1_, i2_, i3_;

    segment(int i1, int i2, int i3 = 0): i3_{i3}{
        if(i1 < i2 ){
            i1_ = i1;
            i2_ = i2;
        }else{
            i1_ = i2;
            i2_ = i1;
        }
    }
};

struct comppair{
    constexpr bool operator()(const segment &p1, const segment &p2) const{
        return (p1.i1_ == p2.i1_  && p1.i2_  == p2.i2_);
    }
};

struct SegmentHash {
    std::size_t operator () (const segment& p) const {
        auto h1 = std::hash<int>{}(p.i1_);
        auto h2 = std::hash<int>{}(p.i2_);
//        auto h3 = std::hash<T1>{}(p[2]);

        // A simple combination hash function
        return h1 ^ h2;
    }
};

int main() {

    double liczbafalowa = 9;
    const int row = 5; //vertices number
    double L = 2.4;
    //vertices number


    double horLen = L / (row - 1);

    std::cout<<"Vertices CREATION BEGIN"<<std::endl;
    std::vector<Vertex2D> vertices;
    std::vector<Point2D> points;


    for(int k =0; k<(row)*(row);k++){
        int i = k/(row);
        int j = k%(row);
//        Point2D p(- L/2 + (horLen * i), - L/2 + (horLen * j ));
        points.emplace_back(- L/2 + (horLen * i), - L/2 + (horLen * j ));
        vertices.emplace_back(points[k], isBorder(k,row));
    }

    Point2D a(- L/2,- L/2), b(L/2,L/2), aa(.001,.001);
    int count=0;
    Wall w(a,b);
    for(auto &p: points )
    if(w.isInsideWall(p)){
        count++;
    }
    std::cout<<count<<std::endl;


    std::ofstream outputFilep("points.txt");
    for(auto & vertex : points){
        outputFilep <<vertex.x << ' ' << vertex.y << std::endl;
    }
    outputFilep.close();
    std::cout<<"Vertices CREATION END"<<std::endl<<std::endl;


    std::vector<std::vector<int>> elementsIdx;

//    pointsRot(points, 1);
//    points[0].y -= 0.1;
    build(points, elementsIdx);


//    std::cout<<elementsIdx.size()<<std::endl;

    if(0){

//TODO przechodzac po vectoerze indeksow elementu jako boki sprawdzamy czy element jest w secie jesli nie liczymy srodek tego boku dodajemy wierzcholek i indeksy tego boku do setu
// set powinien miec indeks wierzcholka pierwszego, drugiego i utworzonego;

        std::unordered_set<segment, SegmentHash, comppair> test_set;
        for(auto &idxs:elementsIdx){
            for(int i=0;i<3;i++){
                if (auto iter = test_set.find({idxs[i],idxs[(i+1)%3]}); iter != test_set.end()){
                    idxs.push_back(iter->i3_);
                }else{

                    Point2D v((vertices[idxs[i]].x() + vertices[idxs[(i+1)%3]].x())/2, (vertices[idxs[i]].y() + vertices[idxs[(i+1)%3]].y())/2);
                    vertices.push_back({v, vertices[idxs[i]].isBorder && vertices[idxs[(i+1)%3]].isBorder});
                    int i3 = static_cast<int>(vertices.size() - 1);
                    test_set.insert({idxs[i],idxs[(i+1)%3],i3});
                    idxs.push_back(i3);
                }
            }
        }
//        for (const auto& element : test_set) {
//            std::cout << element.i1_ << " " << element.i2_ << " " << element.i3_ << std::endl;
//        }

    }

    std::ofstream outputFile("vertices.txt");
    for(auto & vertex : vertices){
        outputFile <<vertex.x() << ' ' << vertex.y() <<  ' ' << vertex.isBorder << std::endl;
    }
    outputFile.close();

    std::ofstream outputFileEle("elementsVerticesIdxes.txt");
    for(auto &e:elementsIdx){
        for(auto &i:e){
            outputFileEle<<i<<" ";
        }
        outputFileEle<<"\n";
    }
    outputFileEle.close();

    const int N = vertices.size();

    std::vector<fem::TriangleElement> Elements;

    std::cout<<"Elements CREATION BEGIN"<<std::endl;

    for(auto &idxs:elementsIdx){
        Elements.emplace_back(0, vertices, idxs, liczbafalowa, fem::LIN);
    }
    std::cout<<"Elements CREATION END"<<std::endl<<std::endl;

//    for(auto &ele:Elements){
//        double x = -1, y =-1;
//        std::cout<<ele.linearBaseFunc[0](x,y)<<" "<<ele.linearBaseFunc[1](x,y)<<" "<<ele.linearBaseFunc[2](x,y)<<" "<<
//        ele.linearBaseFunc[3](x,y)<<" "<<ele.linearBaseFunc[4](x,y)<<" "<<ele.linearBaseFunc[5](x,y)<<std::endl;
//        x = 1; y =-1;
//        std::cout<<ele.linearBaseFunc[0](x,y)<<" "<<ele.linearBaseFunc[1](x,y)<<" "<<ele.linearBaseFunc[2](x,y)<<" "<<
//                 ele.linearBaseFunc[3](x,y)<<" "<<ele.linearBaseFunc[4](x,y)<<" "<<ele.linearBaseFunc[5](x,y)<<std::endl;
//        x = -1, y =1;
//        std::cout<<ele.linearBaseFunc[0](x,y)<<" "<<ele.linearBaseFunc[1](x,y)<<" "<<ele.linearBaseFunc[2](x,y)<<" "<<
//                 ele.linearBaseFunc[3](x,y)<<" "<<ele.linearBaseFunc[4](x,y)<<" "<<ele.linearBaseFunc[5](x,y)<<std::endl;
//        x = 0, y =-1;
//        std::cout<<ele.linearBaseFunc[0](x,y)<<" "<<ele.linearBaseFunc[1](x,y)<<" "<<ele.linearBaseFunc[2](x,y)<<" "<<
//                 ele.linearBaseFunc[3](x,y)<<" "<<ele.linearBaseFunc[4](x,y)<<" "<<ele.linearBaseFunc[5](x,y)<<std::endl;
//        x = 0, y =0;
//        std::cout<<ele.linearBaseFunc[0](x,y)<<" "<<ele.linearBaseFunc[1](x,y)<<" "<<ele.linearBaseFunc[2](x,y)<<" "<<
//                 ele.linearBaseFunc[3](x,y)<<" "<<ele.linearBaseFunc[4](x,y)<<" "<<ele.linearBaseFunc[5](x,y)<<std::endl;
//        x = -1, y =0;
//        std::cout<<ele.linearBaseFunc[0](x,y)<<" "<<ele.linearBaseFunc[1](x,y)<<" "<<ele.linearBaseFunc[2](x,y)<<" "<<
//                 ele.linearBaseFunc[3](x,y)<<" "<<ele.linearBaseFunc[4](x,y)<<" "<<ele.linearBaseFunc[5](x,y)<<std::endl;
//        std::cout<<std::endl;
//    }

    std::cout<<"Tiplets CREATION BEGIN"<<std::endl;

    std::vector<Eigen::Triplet<double>> tripletListS{};

    for(auto &ele:Elements){
        std::vector<std::vector<double>> &E = ele.E_;
        std::vector<int> &glob = ele.globalVectorIdx;
        for (int k=0 ; k < glob.size();k++){
            for(int l = 0; l < glob.size(); l++ ){
                int i = glob[k];
                int j = glob[l];
                if(vertices[i].isBorder){
//                    std::cout<<i<<" "<<j<<" "<<m<<std::endl;
                }
                else{
                    tripletListS.emplace_back(i, j, E[k][l]);
                }

            }
        }
    }

    for(int i = 0; i < vertices.size(); i++){
        if(vertices[i].isBorder){
            tripletListS.emplace_back(i,i,1.);
        }
    }
    std::cout<<"Triplets CREATION END"<<std::endl<<std::endl;



//    Eigen::VectorXd F(N);
    std::vector<double> F(N);
    for(auto &ele:Elements){
        std::vector<double> &Fm = ele.F_;
        std::vector<int> &glob = ele.globalVectorIdx;
        for (int k = 0; k < glob.size(); k++){
                int i = glob[k];
                if(vertices[i].isBorder){
                    F[i] = 0;
                }else{
                    F[i] += Fm[k];
                }
        }
    }

    Eigen::VectorXd EF(N);
    for(int i = 0; i<row*row; i++) EF[i] = F[i];



//    for(int i = 0; i<row; i++){
//        if(i==0 || i==row-1){
//            for(int j=0;j<row;j++){
//                tripletListS.emplace_back(j*row  + i, j*row + i, 1.);
//            }
//        }else
//        {
//            tripletListS.emplace_back(i, i, 1.);
//            tripletListS.emplace_back(row*(row - 1) + i, row*(row - 1) + i, 1.);
//        }
//    }

//    F[(row*row)/2] = 1;// F[(row*row)/2 - row] = 1; F[(row*row)/2 + row] = 1;
//    F[(row*row)/2 + 1] = 1; F[(row*row)/2 - row + 1] = 1; F[(row*row)/2 + row + 1] = 1;
//    F[(row*row)/2 - 1] = 1; F[(row*row)/2 - row- 1] = 1; F[(row*row)/2 + row- 1] = 1;

    bool debug = true;
    if(debug){
        std::vector< std::vector<double> > St (N, std::vector<double>(N));

        for(auto t:tripletListS){
            St[t.row()][t.col()] += t.value();
        }

        FILE* fileS;
        FILE* fileF;
        fileS = fopen("resultS.txt","w");
        fileF = fopen("resultF.txt","w");

        for(int i = 0; i< N; i++){
            fprintf(fileF,"%g\n", F[i]);
            for(int j = 0; j<N; j++)
                fprintf(fileS,"%.2f\t", St[i][j]);
            fprintf(fileS,"\n");
        }
        fclose(fileS);
        fclose(fileF);
    }



    std::cout<<"SparseMatrix CREATION BEGIN"<<std::endl;
    Eigen::SparseMatrix<double> S(N, N);
    S.setFromTriplets(tripletListS.begin(), tripletListS.end());
//    tripletListS.clear();
    std::cout<<"SparseMatrix CREATION END"<<std::endl<<std::endl;

    std::cout<<"Solver CREATION BEGIN"<<std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  solver;

    solver.analyzePattern(S);

    solver.factorize(S);

    std::cout<<"Solver CREATION END"<<std::endl<<std::endl;

    Eigen::VectorXd c(N);

    std::cout<<"Solving BEGIN"<<std::endl;
    c = solver.solve(EF);
    std::cout<<"Solving END"<<std::endl<<std::endl;

    std::ofstream outputFileC("c.txt");
    for (int i = 0; i<N; i ++){
        outputFileC<<c[i]<<std::endl;
    }
    outputFileC.close();



//    double dx = 0.01;
//    int ni = int(L/dx + 1);
//    for(int i = 0; i <ni; i++ ){
//        for(int j = 0; j <ni; j++ ){
//            double x = i*dx - L/2.;
//            double y = j*dx - L/2.;
//            double u = 0;
//            for(auto e:Elements){
//                if(e.ismember(x,y)){
//                    u = e.u_m(x,y, c);
//                    continue;
//                }
//
//            }
//
//            fprintf(file,"%g\t%g\t%g\n",x , y, fabs(u));
//        }
//        fprintf(file,"\n");
//    }
//    fclose(file);
    return 0;
}
