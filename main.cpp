#include <iostream>

#include <thread>
#include <chrono>

#include "MES.h"


int main() {



    const double width = 3; // user input
    const double height = 3;

    mes::solver::Solver s1(width, height, 0.02);
    s1.setImageSize(8020,8020).setBaseFunctionType(mes::QUAD);

    s1
//    .addWall(mes::Wall::createFromDimensions(-1.75, -1.75, .25, 3.5, mes::CONCRETE))
//    .addWall(mes::Wall::createFromDimensions(-1.75, -1.75, 3.5, .25, mes::CONCRETE))
//    .addWall(mes::Wall::createFromDimensions(-1.75, 1.5, 3.5, .25, mes::CONCRETE))
//    .addWall(mes::Wall::createFromDimensions(1.5, -1.75, .25, 3.5, mes::CONCRETE))
    .addWall(mes::Wall::createFromDimensions(-1.5,  0, 1.5, 1.5, mes::CONCRETE));
//            .addWall(mes::Wall::createFromCorners(-1., -.35, -.9, .35,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(-1., -3., -.9, -1.05,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(1.5, -3., 1.6, -1.05,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(-1.2, -1.05, .3, -.95,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(-.4, -1.05, 1., -.95,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(-1., 1.05, -.9, 3.,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(-3.2, -3, -1.5, -2,mes::BRICK))
//            .addWall(mes::Wall::createFromCorners(-3.2, -3, -1.5, -2,mes::BRICK))

//            .addWall(mes::Wall::createFromCorners(1., -3., 1.2, 3.,mes::DAMP))
//            .addWall(mes::Wall::createFromCorners(-2.2, -3., -2., 3.,mes::DAMP))

//            .addWall(mes::Wall::createFromCorners(2.7, 1.7, 3.7, 2.,mes::DAMP))
//            .addWall(mes::Wall::createFromCorners(3.05, 1.35, 3.35, 2.35,mes::DAMP))
            ;

    s1.buildStructure(true);
    s1.computeSolver(true).solve(0.5,-0.5).draw();
//    std::this_thread::sleep_for(std::chrono::seconds(1));
//    s1.drawReal();
//    std::this_thread::sleep_for(std::chrono::seconds(1));
//    s1.drawImag();

    return 0;
}
