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


    double diffQuotient_x(const std::function<double(double&, double&)> &phi, double &x, double &y);
    double diffQuotient_y(const std::function<double(double&, double&)> &phi, double &x, double &y);

    class TriangleElement {
        private:
            static constexpr std::complex<double> air {1,0};
            static constexpr std::complex<double> brick {1.97805,-0.0518468};
            static constexpr std::complex<double> concrete {2.29399,-0.149624};
            static constexpr std::complex<double> damp {2.29399,-0.149624};
            static constexpr std::complex<double> metal {6120.32,-6120.32};

            void initE();
            void initJacob();
            std::complex<double> getRefIdx(elementType &_et);
            double Jacobian(double &zeta, double &eta);
            double map_x(double &zeta, double &eta, int diffFlag = 0);
            double map_y(double &zeta, double &eta, int diffFlag = 0);
            std::pair<double, double> nablaPhik(double &zeta, double &eta, double &jacob, int &k);

    public:
            Vertex2D& globalVector(int &i);
            static double k_;

            mes::ElementIndices &globalVectorIdx;
            std::vector<Vertex2D> &vertices_;
            std::vector<std::vector<std::complex<double>>> E_;
            std::vector<double> jacob_;
            std::vector<std::function<double(double&, double&)>> baseFunc;

            TriangleElement( std::vector<Vertex2D> &ver, ElementIndices &globIndx);
    };
}


#endif //MES_TRIANGLEELEMENT_HPP
