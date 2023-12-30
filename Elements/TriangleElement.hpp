//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_TRIANGLEELEMENT_HPP
#define MES_TRIANGLEELEMENT_HPP

#include <iostream>
#include <vector>
#include <functional>
#include <initializer_list>
#include <cmath>
#include <complex>

#include "Vertex2D.hpp"

namespace fem{
    enum baseFuncType{
        LIN = 1, QUAD = 2
    };

    enum elementType{
        AIR, BRICK, CONCRETE
    };

    class ElementIndices{
    public:
        elementType n;
        std::vector<int> indices;
        ElementIndices(std::vector<int> &_indices, elementType _n);
    };


    double lin_phi0(double &zeta, double &eta);
    double lin_phi1(double &zeta, double &eta);
    double lin_phi2(double &zeta, double &eta);

    double quad_phi0(double &zeta, double &eta);
    double quad_phi1(double &zeta, double &eta);
    double quad_phi2(double &zeta, double &eta);
    double quad_phi3(double &zeta, double &eta);
    double quad_phi4(double &zeta, double &eta);
    double quad_phi5(double &zeta, double &eta);

    double rho(double x, double y);

    double diffQuotient_x(const std::function<double(double&, double&)> &phi, double &x, double &y);
    double diffQuotient_y(const std::function<double(double&, double&)> &phi, double &x, double &y);

    class TriangleElement {
        private:
            static constexpr std::complex<double> air {1,0};
            static constexpr std::complex<double> brick {3.31,5.64923e-09};
            static constexpr std::complex<double> concrete {5.24002,-0.00752353};

            void initE();
            void initF();
            void initJacob();
            std::complex<double> getRefIdx();
            double Jacobian(double &zeta, double &eta);
            double map_x(double &zeta, double &eta, int diffFlag = 0);
            double map_y(double &zeta, double &eta, int diffFlag = 0);
            std::pair<double, double> nablaPhik(double &zeta, double &eta, double &jacob, int &k);

    public:
            int m_;
            Vertex2D globalVector(int &i);
            elementType n_;
            static double k_;
            std::vector<int> &globalVectorIdx;
            std::vector<Vertex2D> *vertices_;
            std::vector<std::vector<std::complex<double>>> E_;
            std::vector<double> F_;
            std::vector<double> jacob_;
            std::vector<std::function<double(double&, double&)>> baseFunc;

            TriangleElement(int m, std::vector<Vertex2D> &ver, std::vector<int> &globIndx, elementType et = AIR, baseFuncType bft = LIN);
    };
}


#endif //MES_TRIANGLEELEMENT_HPP
