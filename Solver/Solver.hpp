//
// Created by Pawulon on 31/12/2023.
//

#ifndef MES_SOLVER_HPP
#define MES_SOLVER_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>

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
        int nDampLayers = 4;
        double frequency = 2.4 ;
        fem::baseFuncType fType = fem::LIN;
        Vertex2D sourcePoint = Vertex2D(0., 0.);

        double widthEleLen ;
        double heightEleLen ;

        unsigned int canvasWidth, canvasHeight;
        double width ;
        double height ;
        int nVerWidth; //vertices number
        int nVerHeight; //vertices number


        std::vector<Vertex2D> points;
        std::vector<Wall> walls;
        std::vector<fem::ElementIndices> elementsIdx;
        std::vector<fem::TriangleElement> Elements;
        int nVertices;
        Eigen::SparseMatrix<std::complex<double>> stiffnessMatrix;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> loadVector;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> solutions;
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>  solverLU;


        std::vector<double> normalizeSolution();
        bool isOnEdge(int k);
        Solver& initDampWalls();



    public:


        Solver(double width, double height, int nVerWidth, int nVerHeight);

        Solver& setNumberOfDampLayers(int i);
        Solver& setFrequency(double f);
        Solver& setBaseFunctionType( fem::baseFuncType fType);
        Solver& setSourcePoint(const Vertex2D &p);
        Solver& setSourcePoint(double x, double y);
        Solver& setImageSize(unsigned int x, unsigned int y);



        Solver& generateSimpleMesh();


        Solver& addWall(double leftDownX, double leftDownY, double w, double h, fem::elementType elt);
        Solver& addWall(Wall &w);

        Solver& divideIntoElements();
        Solver& buildElements();
        Solver& buildStiffnessMatrix();
        Solver& buildLoadVector();
        Solver& buildLoadVector(double sx, double sy);
        Solver& buildLoadVector(const Vertex2D &p);
        Solver& buildSolver();
        Solver& solve();
        std::string draw();
        Solver& doAll();
    };


#endif //MES_SOLVER_HPP
} }
