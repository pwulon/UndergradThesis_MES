//
// Created by Pawulon on 31/12/2023.
//

#ifndef MES_SOLVER_HPP
#define MES_SOLVER_HPP
#define MEASURE_TIME(func) \
        {auto start_time = std::chrono::high_resolution_clock::now(); \
        func; \
        auto end_time = std::chrono::high_resolution_clock::now(); \
        auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count(); \
        std::cout << "Time taken by " #func ": " << elapsed_time << " microseconds." << std::endl;}


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


namespace mes {
    namespace solver{

    class Solver {
    private:
        int nDampLayers = 4;
        double frequency = 2.4 ; // in GHz
        mes::baseFuncType fType = mes::LIN;
        Vertex2D sourcePoint = Vertex2D(0., 0.);
        double maxSolutionValue = 0;

        double widthEleLen ;
        double heightEleLen ;

        unsigned int canvasWidth, canvasHeight;
        double width ;
        double height ;
        int nVerWidth; //vertices number
        int nVerHeight; //vertices number



        std::vector<Vertex2D> points;
        std::vector<Wall> walls;
        std::vector<mes::ElementIndices> elementsIdx;
        std::vector<mes::TriangleElement> Elements;
        int nVertices, nVerticesQuad = 0;
        Eigen::SparseMatrix<std::complex<double>> stiffnessMatrix;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> loadVector;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> solutions;
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>  solverLU;


        std::vector<double> normalizeSolution();
        bool isOnEdge(int &i, int &j);
        Solver& generateSimpleMesh();
        Solver& initDampWalls();
        Solver& divideIntoElements();
        Solver& buildElements();
        Solver& buildStiffnessMatrix();
        Solver& buildLoadVector();
        Solver& buildSolver();

    public:
        Solver(double width, double height, int nVerWidth, int nVerHeight);

        Solver& setNumberOfDampLayers(int i);
        Solver& setFrequency(double f);
        Solver& setBaseFunctionType(mes::baseFuncType fType);
        Solver& setSourcePoint(double x, double y);
        Solver& setImageSize(unsigned int x, unsigned int y);

        double getMaxValue() const;

        Solver& addWall(Wall w);

        Solver& buildStructure(bool compTime = false);
        Solver& computeSolver(bool compTime = false);

        Solver& solve();
        Solver& solve(double sx, double sy);
        Solver& draw();
    };


#endif //MES_SOLVER_HPP
} }
