#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include <thread>
#include <chrono>

#include "Solver/Solver.hpp"


int main() {

    fem::baseFuncType fType = fem::LIN;

    const int nVerWidth = 521; //vertices number
    const int nVerHeight = 361; //vertices number


    double width = 5.2;
    double height = 3.6;
    //vertices number

    Vertex2D sPoint(0, 0);
    std::cout<<"Constructor"<<std::endl;
    fem::solve::Solver s1(width, height, nVerWidth, nVerHeight, sPoint, fType);

    std::cout<<"setNumberOfDampLayers"<<std::endl;
    s1.setNumberOfDampLayers(4);

    std::cout<<"generateSimpleMesh"<<std::endl;
    s1.generateSimpleMesh();

    std::cout<<"initDampWalls"<<std::endl;
    s1.initDampWalls();
    Wall w1(1., .33, .2, height, fem::DAMP);
    Wall w3(1.75, .33, width, .2, fem::DAMP);
    Wall w4(.5, -1., width, .2, fem::DAMP);
    Wall w2(-width*.5, -height*.5, .2, height, fem::DAMP);
    w2.rot(45);
    s1.addWall(w1).addWall(w2).addWall(w3).addWall(w4);


    std::cout<<"divideIntoElements"<<std::endl;
    s1.divideIntoElements();

    std::cout<<"buildElements"<<std::endl;
    s1.buildElements();

    std::cout<<"buildStiffnessMatrix"<<std::endl;
    s1.buildStiffnessMatrix();

    std::cout<<"buildSolver"<<std::endl;
    s1.buildSolver();

    std::cout<<"buildLoadVector"<<std::endl;
    s1.buildLoadVector(1., 0.);

    std::cout<<"solve"<<std::endl;
    s1.solve();

    std::cout<<"draw"<<std::endl;
    s1.draw();

//    std::cout<<"buildLoadVector"<<std::endl;

    for(int t = 1; t<25;t++){
        double x = 1.*cos(t*2*M_PI/25.);
        double y = 1.*sin(t*2*M_PI/25.);
        s1.buildLoadVector(x, y).solve().draw();
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }


    return 0;
}
