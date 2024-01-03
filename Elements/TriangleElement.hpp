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

namespace mes{
    enum baseFuncType{
        LIN = 1, QUAD = 2
    };

    enum elementType{
        AIR, BRICK, CONCRETE, DAMP
    };

    class ElementIndices{
    public:
        elementType et;
        baseFuncType bft;
        std::vector<int> indices;
        ElementIndices(std::vector<int> &_indices, elementType _et, baseFuncType _ft = mes::LIN);
    };


    inline double lin_phi0(double &zeta, double &eta);
    inline double lin_phi1(double &zeta, double &eta);
    inline double lin_phi2(double &zeta, double &eta);

    inline double quad_phi0(double &zeta, double &eta);
    inline double quad_phi1(double &zeta, double &eta);
    inline double quad_phi2(double &zeta, double &eta);
    inline double quad_phi3(double &zeta, double &eta);
    inline double quad_phi4(double &zeta, double &eta);
    inline double quad_phi5(double &zeta, double &eta);

    inline double rho(double x, double y);

    double diffQuotient_x(const std::function<double(double&, double&)> &phi, double &x, double &y);
    double diffQuotient_y(const std::function<double(double&, double&)> &phi, double &x, double &y);

    class TriangleElement {
        private:
            static constexpr std::complex<double> air {1,0};
            static constexpr std::complex<double> brick {3.31,5.64923e-09};
            static constexpr std::complex<double> concrete {1, 1};
            static constexpr std::complex<double> damp {1, 1};

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
            static Vertex2D source_;
            mes::ElementIndices &globalVectorIdx;
            std::vector<Vertex2D> *vertices_;
            std::vector<std::vector<std::complex<double>>> E_;
            std::vector<double> F_;
            std::vector<double> jacob_;
            std::vector<std::function<double(double&, double&)>> baseFunc;

            TriangleElement(int m, std::vector<Vertex2D> &ver, ElementIndices &globIndx);
    };
}


#endif //MES_TRIANGLEELEMENT_HPP
