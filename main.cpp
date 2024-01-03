#include <iostream>

#include <thread>

#include "MES.h"

int main() {
    const int nVerWidth = 271; //vertices number
    const int nVerHeight = 271; //vertices number

    const double width = 5.2;
    const double height = 5.2;

    std::cout<<"Constructor"<<std::endl;
    mes::solver::Solver s1(width, height, nVerWidth, nVerHeight);
    s1.setImageSize(1600,1600).setBaseFunctionType(mes::LIN);


    std::cout<<"initDampWalls"<<std::endl;
    Wall w1(1., .33, .2, height, mes::DAMP);
    Wall w3(1.75, .33, width, .2, mes::DAMP);
    Wall w4(.5, -1., width, .2, mes::DAMP);
    Wall w5(-.5, -.5, .3, 1., mes::DAMP);
    Wall w2(-width*.5, -height*.5, .2, height, mes::DAMP);
    w2.rot(45);
    s1.addWall(w1).addWall(w2).addWall(w3).addWall(w4).addWall(w5);

    std::cout<<"generateSimpleMesh"<<std::endl;
    s1.generateSimpleMesh();

    std::cout<<"divideIntoElements"<<std::endl;
    s1.divideIntoElements();

    std::cout<<"buildElements"<<std::endl;
    s1.buildElements();

    std::cout<<"buildStiffnessMatrix"<<std::endl;
    s1.buildStiffnessMatrix();

    std::cout<<"buildSolver"<<std::endl;
    s1.buildSolver();

    std::cout<<"buildLoadVector"<<std::endl;
    s1.buildLoadVector(1.5, 0.);

    std::cout<<"solve"<<std::endl;
    s1.solve();

    std::cout<<"draw"<<std::endl;
    s1.draw();


    return 0;
}
