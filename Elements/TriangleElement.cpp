//
// Created by Pawulon on 13/12/2023.
//

#include "TriangleElement.hpp"


namespace fem {
    std::vector<double> xgaus = {-0.3333333, -0.0597158717, -0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853,
                                 0.5948539707};
    std::vector<double> ygaus  = {-0.3333333, -0.0597158717, -0.8805682564, -0.0597158717, -0.7974269853, 0.5948539707,
                                  -0.7974269853};
    std::vector<double> wgaus  = {0.45, 0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};


    TriangleElement::TriangleElement(int m, std::vector<Vertex2D> &ver,
                                     std::initializer_list<int> globIndx,
                                     double k) : m_{m}, vertices_{&ver}, globalVectorIdx{globIndx}, k_{k} {
        F_ = std::vector<double>(3);
        E_ = std::vector<std::vector<double>>(3, std::vector<double>(3));
        initE();
        initF();
    }

    void TriangleElement::initE() {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 7; k++) {
                    std::pair<double, double> nablaphi_i = nablaPhik(xgaus[k], ygaus[k], i);
                    std::pair<double, double> nablaphi_j = nablaPhik(xgaus[k], ygaus[k], j);
                    E_[i][j] += wgaus[k] * Jacobian(fem::xgaus[k], fem::ygaus[k])
                                * (-(nablaphi_i.first * nablaphi_j.first
                                     + nablaphi_i.second * nablaphi_j.second) +
                                   pow(k_, 2) *
                                   linearBaseFunc[i](fem::xgaus[k], fem::ygaus[k]) * linearBaseFunc[j](fem::xgaus[k], fem::ygaus[k]));
                }
//                if (m_ < 2)std::cout << E_[i][j] << "\t"; //debug
            }
//            if (m_ < 2)std::cout << std::endl; //debug
        }
//        if (m_ < 2)std::cout << std::endl; //debug
    }

    void TriangleElement::initF() {
        for (int j = 0; j < 3; j++) {
            if(!globalVector(j).isBorder){
                for (int k = 0; k < 7; k++) {
                    F_[j] += wgaus[k] * Jacobian(xgaus[k], ygaus[k]) *
                             linearBaseFunc[j](xgaus[k], ygaus[k]) *
                             rho(map_x(xgaus[k], ygaus[k]), map_y(xgaus[k], ygaus[k]));
                }
            }else{
                F_[j] = 0;
            }

//            if (m_ < 2)std::cout << F_[j] << "\t"; //debug
        }
//        if (m_ < 2)std::cout << std::endl << std::endl; //debug
    }

    double TriangleElement::Jacobian(double &zeta, double &eta) {
        return map_x(zeta, eta, 1) *
               map_y(zeta, eta,2) -
               map_y(zeta, eta,1) *
               map_x(zeta, eta, 2);
    }

    std::pair<double, double> TriangleElement::nablaPhik(double &zeta, double &eta, int &k) {
        return {diffQuotient_x(linearBaseFunc[k], zeta, eta) *
                map_y(zeta, eta,2) / Jacobian(zeta, eta)
                - diffQuotient_y(linearBaseFunc[k], zeta, eta) * map_y(zeta, eta,1) / Jacobian(zeta, eta),
                -diffQuotient_x(linearBaseFunc[k], zeta, eta) * map_x(zeta, eta,2) / Jacobian(zeta, eta)
                + diffQuotient_y(linearBaseFunc[k], zeta, eta) * map_x(zeta, eta, 1) / Jacobian(zeta, eta)};
    }

    double TriangleElement::map_x(double &zeta, double &eta, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            switch (diffFlag) {
                case 0:
                    sum += globalVector(i).x() * linearBaseFunc[i](zeta, eta);
                    break;
                case 1:
                    sum += globalVector(i).x() * diffQuotient_x(linearBaseFunc[i],zeta,eta);
                    break;
                case 2:
                    sum += globalVector(i).x() * diffQuotient_y(linearBaseFunc[i],zeta,eta);
                    break;
            }
        }
        return sum;
    }

    double TriangleElement::map_y(double &zeta, double &eta, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            switch (diffFlag) {
                case 0:
                    sum += globalVector(i).y() * linearBaseFunc[i](zeta, eta);
                    break;
                case 1:
                    sum += globalVector(i).y() * diffQuotient_x(linearBaseFunc[i], zeta, eta);
                    break;
                case 2:
                    sum += globalVector(i).y() * diffQuotient_y(linearBaseFunc[i], zeta, eta);
                    break;
            }
        }
        return sum;
    }

    Vertex2D TriangleElement::globalVector(int &i){
        return (*vertices_)[globalVectorIdx[i]];
    }

    double phi0(double &zeta, double &eta) {
        return -.5 * (zeta + eta);
    }

    double phi1(double &zeta, double &eta) {
        return .5 * (1 + zeta);
    }

    double phi2(double &zeta, double &eta) {
        return .5 * (1 + eta);
    }

    double rho(double x, double y){
        double a = 1./16.;
        return (1./(fabs(a) * sqrt(M_PI)))*exp(-1.*(pow(x/a,2)+pow(y/a,2)));
    }

    double diffQuotient_x(const std::function<double(double&, double&)> &phi, double &x, double &y) {
        double delta = 1e-4;
        double x1 = x + delta;
        double x2 = x - delta;
        return (phi(x1, y) - phi(x2, y)) / (2 * delta);
    }

    double diffQuotient_y(const std::function<double(double&, double&)> &phi, double &x, double &y) {
        double delta = 1e-4;
        double y1 = y + delta;
        double y2 = y - delta;
        return (phi(x, y1) - phi(x, y2)) / (2 * delta);
    }
}