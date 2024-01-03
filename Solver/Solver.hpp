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
        int nDampLayers;
        std::vector<Vertex2D> points;
        std::vector<Wall> walls;
        std::vector<fem::ElementIndices> elementsIdx;
        std::vector<fem::TriangleElement> Elements;
        fem::baseFuncType fType;
        int nVertices;
        Eigen::SparseMatrix<std::complex<double>> stiffnessMatrix;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> loadVector;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> solutions;
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>  solverLU;

        double frequency;
        std::vector<double> normalizeSolution();
        bool isOnEdge(int k);

        unsigned int canvasWidth, canvasHeight;
        double width ;
        double height ;
        int nVerWidth; //vertices number
        int nVerHeight; //vertices number

        double widthEleLen ;
        double heightEleLen ;

    public:
        Vertex2D sourcePoint;
//        Solver();

        Solver(double width, double height, int nVerWidth, int nVerHeight, const Vertex2D &sourcePoint,
               fem::baseFuncType fType = fem::LIN,  double frequency = 2.4 * pow(10,9));

        Solver& setNumberOfDampLayers(int i);

        Solver& generateSimpleMesh();

        Solver& initDampWalls();
        Solver& addWall(double leftDownX, double leftDownY, double w, double h, fem::elementType elt);
        Solver& addWall(Wall &w);

        Solver& divideIntoElements();
        Solver& buildElements();
        Solver& buildStiffnessMatrix();
        Solver& buildLoadVector();
        Solver& buildLoadVector(double sx, double sy);
        Solver& buildSolver();
        Solver& solve();
        Solver& draw();
        Solver& doAll();
    };


#endif //MES_SOLVER_HPP
} }
