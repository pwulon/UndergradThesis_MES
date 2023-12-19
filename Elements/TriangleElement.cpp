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
                                     std::vector<int> &globIndx,
                                     double k, type t) : m_{m}, vertices_{&ver}, globalVectorIdx{globIndx}, k_{k} {
        if(type(t) == 1){
            std::cout<<"lin"<<std::endl;
            linearBaseFunc  = {&lin_phi0, &lin_phi1, &lin_phi2};
        } else if(type(t) == 2){
            std::cout<<"quad"<<std::endl;
            linearBaseFunc = {&quad_phi0, &quad_phi1, &quad_phi2, &quad_phi3, &quad_phi4, &quad_phi5};
        }
        F_ = std::vector<double>(linearBaseFunc.size());
        E_ = std::vector<std::vector<double>>(linearBaseFunc.size(), std::vector<double>(linearBaseFunc.size()));
        initJacob();
        initE();
        initF();

//        double x,y, x_ = 0, y_=0;
//        Point2D p(map_x( x_, y_),map_y( x_, y_));
//
//        vertices_->push_back(Vertex2D(p));
//        globalVectorIdx.push_back(vertices_->size() - 1);
    }

    void TriangleElement::initE() {
        for (int i = 0; i < linearBaseFunc.size(); i++) {
            for (int j = 0; j < linearBaseFunc.size(); j++) {
                for (int k = 0; k < 7; k++) {
                    std::pair<double, double> nablaphi_i = nablaPhik(xgaus[k], ygaus[k], jacob_[k], i);
                    std::pair<double, double> nablaphi_j = nablaPhik(xgaus[k], ygaus[k], jacob_[k], j);
                    E_[i][j] += wgaus[k] * jacob_[k] * (
                            -(nablaphi_i.first * nablaphi_j.first +
                            nablaphi_i.second * nablaphi_j.second) +
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
        for (int j = 0; j <  linearBaseFunc.size(); j++) {
            if(!globalVector(j).isBorder){
                for (int k = 0; k < 7; k++) {
                    F_[j] += wgaus[k] * jacob_[k] *
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

    std::pair<double, double> TriangleElement::nablaPhik(double &zeta, double &eta, double &jacob, int &k) {
        return {diffQuotient_x(linearBaseFunc[k], zeta, eta) * map_y(zeta, eta,2) / jacob
                - diffQuotient_y(linearBaseFunc[k], zeta, eta) * map_y(zeta, eta,1) / jacob,
                -diffQuotient_x(linearBaseFunc[k], zeta, eta) * map_x(zeta, eta,2) / jacob
                + diffQuotient_y(linearBaseFunc[k], zeta, eta) * map_x(zeta, eta, 1) / jacob};
    }

    double TriangleElement::map_x(double &zeta, double &eta, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < linearBaseFunc.size(); i++) {
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
        for (int i = 0; i < linearBaseFunc.size(); i++) {
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

    void TriangleElement::initJacob() {
        for(int k = 0; k<7; k++){
            jacob_.push_back(Jacobian(fem::xgaus[k],fem::ygaus[k]));
        }
    }

    double lin_phi0(double &zeta, double &eta) {
        return -.5 * (zeta + eta);
    }

    double lin_phi1(double &zeta, double &eta) {
        return .5 * (1 + zeta);
    }

    double lin_phi2(double &zeta, double &eta) {
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

    double quad_phi0(double &zeta, double &eta) {
        double L = lin_phi0(zeta, eta);
        return  L * (2. * L - 1);
    }

    double quad_phi1(double &zeta, double &eta) {
        double L = lin_phi1(zeta, eta);
        return  L * (2. * L - 1);
    }

    double quad_phi2(double &zeta, double &eta) {
        double L = lin_phi2(zeta, eta);
        return  L * (2. * L - 1);
    }

    double quad_phi3(double &zeta, double &eta) {
        return  4. * lin_phi0(zeta, eta)*lin_phi1(zeta, eta);

    }

    double quad_phi4(double &zeta, double &eta) {
        return  4. * lin_phi1(zeta, eta)*lin_phi2(zeta, eta);

    }

    double quad_phi5(double &zeta, double &eta) {
        return  4. * lin_phi0(zeta, eta)*lin_phi2(zeta, eta);
    }

}