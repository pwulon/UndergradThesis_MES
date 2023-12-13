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
    double phi0(double &zeta, double &eta);
    double phi1(double &zeta, double &eta);
    double phi2(double &zeta, double &eta);
    double rho(double x, double y);

    double diffQuotient_x(const std::function<double(double&, double&)> &phi, double &x, double &y);
    double diffQuotient_y(const std::function<double(double&, double&)> &phi, double &x, double &y);

    class TriangleElement {
        private:
            void initE();
            void initF();
            double Jacobian(double &zeta, double &eta);
            double map_x(double &zeta, double &eta, int diffFlag = 0);
            double map_y(double &zeta, double &eta, int diffFlag = 0);
            std::pair<double, double> nablaPhik(double &zeta, double &eta, int &k);
            Vertex2D globalVector(int &i);
    public:
            int m_;
            double k_;
            std::vector<int> globalVectorIdx{};
            const std::vector<Vertex2D> *vertices_;
            std::vector<std::vector<double>> E_;
            std::vector<double> F_;
            std::vector<std::function<double(double&, double&)>> linearBaseFunc {&phi0, &phi1, &phi2};

            TriangleElement(int m, std::vector<Vertex2D> &ver, std::initializer_list<int> globIndx, double k = 1);
    };
}


#endif //MES_TRIANGLEELEMENT_HPP
