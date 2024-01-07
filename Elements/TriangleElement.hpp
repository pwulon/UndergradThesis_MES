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
        LIN = 0, QUAD = 1
    };

    enum elementType{
        DAMP = -1,
        AIR = 0,
        BRICK = 1,
        CONCRETE = 2,
        STEEL = 3,
        GLASS = 4
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

    static std::complex<double> air;
    static std::complex<double> brick;
    static std::complex<double> concrete;
    static std::complex<double> damp;
    static std::complex<double> steel;
    static std::complex<double> glass;

    void setRefIdx(double f);


    class TriangleElement {
        private:
            static int ngaus;

            static std::vector<std::function<double(double&, double&)>> baseFunc;
            static std::vector<std::vector<double>>  baseFunc_Values;
            static std::vector<std::vector<double>>  baseFuncDiffQuotient_xValues;
            static std::vector<std::vector<double>>  baseFuncDiffQuotient_yValues;


            void initE();
            void initJacob();
            std::complex<double> getRefIdx(elementType &_et);
            double Jacobian(double &zeta, double &eta);
            double map_x(double &zeta, double &eta, int diffFlag = 0);
            double map_y(double &zeta, double &eta, int diffFlag = 0);
            double map_x(int k, int diffFlag = 0);
            double map_y(int k, int diffFlag = 0);
            std::pair<double, double> nablaPhik(double &zeta, double &eta, double &jacob, int &k);
            std::pair<double, double> nablaPhik(int &i, int &k);

    public:
        static std::vector<double> xgaus, ygaus, wgaus;
            Vertex2D& globalVector(int &i);
            static double k_;

            mes::ElementIndices &globalVectorIdx;
            std::vector<Vertex2D> &vertices_;
            std::vector<std::vector<std::complex<double>>> E_;
            std::vector<double> jacob_;


            TriangleElement( std::vector<Vertex2D> &ver, ElementIndices &globIndx);
            static void staticMembersInit(baseFuncType T);
    };
}


#endif //MES_TRIANGLEELEMENT_HPP
