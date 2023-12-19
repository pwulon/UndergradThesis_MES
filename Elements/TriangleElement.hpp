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

#include "Vertex2D.hpp"

namespace fem{
    enum type{
        LIN = 1, QUAD = 2
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
            void initE();
            void initF();
            void initJacob();
            double Jacobian(double &zeta, double &eta);
            double map_x(double &zeta, double &eta, int diffFlag = 0);
            double map_y(double &zeta, double &eta, int diffFlag = 0);
            std::pair<double, double> nablaPhik(double &zeta, double &eta, double &jacob, int &k);
            Vertex2D globalVector(int &i);
    public:
            int m_;
            double k_;
            std::vector<int> &globalVectorIdx;
            std::vector<Vertex2D> *vertices_;
            std::vector<std::vector<double>> E_;
            std::vector<double> F_;
            std::vector<double> jacob_;
            std::vector<std::function<double(double&, double&)>> linearBaseFunc;

            TriangleElement(int m, std::vector<Vertex2D> &ver, std::vector<int> &globIndx, double k = 1, type t = LIN);
    };
}


#endif //MES_TRIANGLEELEMENT_HPP
