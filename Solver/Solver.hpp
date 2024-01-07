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
#include <eigen3/Eigen/SparseCholesky>


namespace mes {
    namespace solver{

    class Solver {
    private:
        double dx;
        int nDampLayers = 4;
        int minIdx = 0;
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
        int nVertices, nVerticesLin = 0;
        Eigen::SparseMatrix<std::complex<double>> stiffnessMatrix;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> loadVector;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> solutions;
        std::vector<double> sum_solutions;
        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>  solverLU;


        std::vector<double> normalizeSolutionAbs();
        std::vector<double> normalizeSolutionReal();
        std::vector<double> normalizeSolutionImg();
        bool isOnEdge(int &i, int &j);
        Solver& generateSimpleMesh();
        Solver& initDampWalls();
        Solver& divideIntoElements();
        Solver& buildElements();
        Solver& buildStiffnessMatrix();
        Solver& buildLoadVector();
        Solver& buildSolver();

    public:
        Solver(mes::baseFuncType _bft, double _f);
        Solver& setMeshSize(double _width, double _height, int _nVerWidth, int _nVerHeight);
        Solver& setMeshSize(double _width, double _height, double _dx);

        Solver& setDampLayersDepth(double d);

        Solver& setSourcePoint(double x, double y);

        Solver& setImageSize(unsigned int x, unsigned int y);

        double getMaxValue() const;

        Solver& addWall(Wall w);

        Solver& buildStructure(bool compTime = false);
        Solver& computeSolver(bool compTime = false);

        Solver& solve();
        Solver& solve(double sx, double sy, bool add = false);
        std::string draw();
        std::string drawImag();
        std::string drawReal();
    };


#endif //MES_SOLVER_HPP
} }
