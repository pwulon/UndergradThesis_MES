//
// Created by Pawulon on 31/12/2023.
//

#ifndef MES_SOLVER_HPP
#define MES_SOLVER_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "../Elements/TriangleElement.hpp"
#include "../Elements/Section.hpp"
#include "../Elements/Wall.hpp"

#include "../FortunesAlgo/Fortunes.hpp"

#include "../Visualization/CreateOpenGlWindow.hpp"

#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseLU>


namespace fem {
    namespace solve{
class Solver {
private:


public:
    std::vector<Vertex2D> points;

    unsigned int canvasWidth, canvasHeight;
    double width = 3.2;
    double height = 3.2;
    const int nVerWidth = 641; //vertices number
    const int nVerHeight = 641; //vertices number

    bool isOnEdge(int k);
    void generateSimpleMesh();


};


#endif //MES_SOLVER_HPP
} }
