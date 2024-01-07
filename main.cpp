#include <iostream>

#include <thread>
#include <chrono>

#include "MES.h"

//preview 3
int main() {
    const double width = 5;
    const double height = 5;
    double frequency = 5;
    mes::baseFuncType b = mes::baseFuncType::LIN; // LIN, QUAD

    mes::solver::Solver s1(b, frequency);
    s1.setMeshSize(width, height, .01).setDampLayersDepth(0.25).setImageSize(1200,1200);

    mes::elementType walltype = mes::elementType::CONCRETE;

    for(int i = 0; i<180; i++){
        s1.addWall(mes::Wall::createFromCorners(-1, -0.05, 1, .05, walltype).rotSelf(i));
    }

    s1.addWall(mes::Wall::createFromCorners(-3, 2, -2, 3, walltype).rotSelf(45));
    s1.addWall(mes::Wall::createFromCorners(2, -3, 3, -2, walltype).rotSelf(45));

    s1.buildStructure(true);
    s1.computeSolver(true);

    auto filename = s1.solve(-2.,-2., false).draw();

    return 0;
}

// preview 1
//{
//s1.addWall(mes::Wall::createFromCorners(-1, -3.1, 0., -2.9, walltype));
//s1.addWall(mes::Wall::createFromCorners(-4, -.1, -.1, .1, walltype))
//        .addWall(mes::Wall::createFromCorners(-.3, -.4, -.1, .4, walltype))
//        .addWall(mes::Wall::createFromCorners(-.3, -3., -.1, - 1.4, walltype))
//        .addWall(mes::Wall::createFromCorners(-.3, 1.4, -.1, 4., walltype))
//
//        .addWall(mes::Wall::createFromCorners(2., 1., 3., 2., walltype))
//
//        .addWall(mes::Wall::createFromCorners(1., -3., 1.2, 0., walltype))
//        .addWall(mes::Wall::createFromCorners(1, -.2, 1.5, 0., walltype))
//        .addWall(mes::Wall::createFromCorners(2.5, -.2, 4, 0., walltype))
//
//        .addWall(mes::Wall::createFromCorners(2., -3, 3.5, -2.25, walltype))
//
//        .addWall(mes::Wall::createFromCorners(-5., 2, -3., 4, walltype).rotSelf(45))
//
//        .addWall(mes::Wall::createFromCorners(-2, -3, -1.8, -2.55, walltype))
//        .addWall(mes::Wall::createFromCorners(-2, -1.55, -1.8, -.1, walltype))
//
//s1.addWall(mes::Wall::  (-1, -3.1, 0., -2.9, walltype));
//s1.addWall(mes::Wall::createFromCorners(0., -2.1, 1., -1.9, walltype));
//
//s1.addWall(mes::Wall::createFromCorners(-1, -1.1, 0., -0.9, walltype));
//s1.addWall(mes::Wall::createFromCorners(0., -0.1, 1., .1, walltype));
//
//s1.addWall(mes::Wall::createFromCorners(-1, .9, 0., 1.1, walltype));
//s1.addWall(mes::Wall::createFromCorners(0., 1.9, 1., 2.1, walltype));
//
//s1.addWall(mes::Wall::createFromCorners(-1, 2.9, 0., 3.1, walltype));
//s1.addWall(mes::Wall::createFromCorners(-1, -0.1, 0., .1, walltype));
//}