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
    s1.setImageSize(1600,1600).setBaseFunctionType(mes::LIN).setFrequency(2.4);


    mes::Wall w1 = mes::Wall::createFromDimensions(1., .33, .2, height, mes::DAMP);
    mes::Wall w2 = mes::Wall::createFromDimensions(-width*.5, -height*.5, .2, height, mes::DAMP).rotSelf(45);
    mes::Wall w3 = mes::Wall::createFromDimensions(1.75, .33, width, .2, mes::DAMP);

    s1.addWall(w1)
        .addWall(w2)
        .addWall(w3)
        .addWall(mes::Wall::createFromDimensions(.5, -1., width, .2))
        .addWall(mes::Wall::createFromCorners(-.5, -.5, .5, .5));

    s1.buildStructure(true).computeSolver(true).solve(1.5, 0);
    std::cout << "abs result range: [" << 0 << ", " << s1.getMaxValue() << "]\n";
    std::cout << "draw" <<std::endl;
    s1.draw();


    return 0;
}
