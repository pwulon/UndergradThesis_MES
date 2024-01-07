//
// Created by Pawulon on 13/12/2023.
//

#include "../Elements/TriangleElement.hpp"


namespace mes {

    double lin_phi0(double &zeta, double &eta) {
        return -.5 * (zeta + eta);
    }

    double lin_phi1(double &zeta, double &eta) {
        return .5 * (1 + zeta);
    }

    double lin_phi2(double &zeta, double &eta) {
        return .5 * (1 + eta);
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



    ElementIndices::ElementIndices(std::vector<int> &_indices, elementType _et,  baseFuncType _ft):indices{_indices}, et{_et}, bft{_ft} {};

    double TriangleElement::k_ = 1;
    std::vector<std::function<double(double&, double&)>> TriangleElement::baseFunc;
    std::vector<std::vector<double>> TriangleElement::baseFunc_Values;
    std::vector<std::vector<double>> TriangleElement::baseFuncDiffQuotient_xValues;
    std::vector<std::vector<double>> TriangleElement::baseFuncDiffQuotient_yValues;

    int TriangleElement::ngaus = 7;
    std::vector<double> TriangleElement::xgaus = {-0.3333333, -0.0597158717, -0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853,
                                 0.5948539707};
    std::vector<double> TriangleElement::ygaus  = {-0.3333333, -0.0597158717, -0.8805682564, -0.0597158717, -0.7974269853, 0.5948539707,
                                  -0.7974269853};
    std::vector<double> TriangleElement::wgaus  = {0.45, 0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};



    void TriangleElement::staticMembersInit(baseFuncType T) {
        switch (T) {
            case LIN:
                baseFunc  = {lin_phi0, lin_phi1, lin_phi2};
                break;
            case QUAD:
                baseFunc = {&quad_phi0, &quad_phi1, &quad_phi2, &quad_phi3, &quad_phi4, &quad_phi5};
                //                baseFunc = {&quad_phi5, &quad_phi0, &quad_phi3, &quad_phi1, &quad_phi4, &quad_phi2};
                break;
        }
        baseFunc_Values = std::vector<std::vector<double>>(baseFunc.size(), std::vector<double>(ngaus));
        baseFuncDiffQuotient_xValues = std::vector<std::vector<double>>(baseFunc.size(), std::vector<double>(ngaus));
        baseFuncDiffQuotient_yValues = std::vector<std::vector<double>>(baseFunc.size(), std::vector<double>(ngaus));

        for(int i = 0; i< baseFunc.size(); i++){
            for(int k = 0 ; k < ngaus; k++){
                baseFunc_Values[i][k] = baseFunc[i](xgaus[k], ygaus[k]);
                baseFuncDiffQuotient_xValues[i][k] = diffQuotient_x(baseFunc[i], xgaus[k], ygaus[k]);
                baseFuncDiffQuotient_yValues[i][k] = diffQuotient_y(baseFunc[i], xgaus[k], ygaus[k]);
            }
        }
    }


    TriangleElement::TriangleElement(std::vector<Vertex2D> &ver, ElementIndices &globIndx ) : vertices_{ver}, globalVectorIdx{globIndx} {
        initJacob();
        initE();
    }

    void TriangleElement::initE() {
        E_ = std::vector<std::vector<std::complex<double>>>(baseFunc.size(), std::vector<std::complex<double>>(baseFunc.size()));
        for (int i = 0; i < baseFunc.size(); i++) {
            for (int j = 0; j < baseFunc.size(); j++) {
                for (int k = 0; k < ngaus; k++) {
                    std::pair<double, double> nablaphi_i = nablaPhik(i, k);
                    std::pair<double, double> nablaphi_j = nablaPhik(j, k);
                    E_[i][j] += wgaus[k] * jacob_[k] * (
                                -(nablaphi_i.first * nablaphi_j.first
                                + nablaphi_i.second * nablaphi_j.second)
//                            pow(k_, 2) * pow(getRefIdx(globalVector(i).et), 2) *

                                + pow(k_, 2) * getRefIdx(globalVectorIdx.et) *
                                baseFunc_Values[i][k]  * baseFunc_Values[j][k]
                                );
                }
//                std::cout<<E_[i][j]<<" ";
            }
//            std::cout<<std::endl;
        }
    }

    double TriangleElement::Jacobian(double &zeta, double &eta) {
        return map_x(zeta, eta, 1) *
               map_y(zeta, eta, 2) -
               map_y(zeta, eta,1) *
               map_x(zeta, eta, 2);
    }

    void TriangleElement::initJacob() {
        for(int k = 0; k < ngaus; k++){
            jacob_.push_back(Jacobian(xgaus[k], ygaus[k]));
        }
    }

    std::pair<double, double> TriangleElement::nablaPhik(double &zeta, double &eta, double &jacob, int &k) {
        return {diffQuotient_x(baseFunc[k], zeta, eta) * map_y(zeta, eta, 2) / jacob
                - diffQuotient_y(baseFunc[k], zeta, eta) * map_y(zeta, eta, 1) / jacob,
                -diffQuotient_x(baseFunc[k], zeta, eta) * map_x(zeta, eta, 2) / jacob
                + diffQuotient_y(baseFunc[k], zeta, eta) * map_x(zeta, eta, 1) / jacob};
    }

    std::pair<double, double> TriangleElement::nablaPhik(int &i, int &k) {
        return {baseFuncDiffQuotient_xValues[i][k] * map_y(k, 2) / jacob_[k]
                - baseFuncDiffQuotient_yValues[i][k] * map_y(k, 1) / jacob_[k],
                -baseFuncDiffQuotient_xValues[i][k] * map_x(k, 2) / jacob_[k]
                + baseFuncDiffQuotient_yValues[i][k] * map_x(k, 1) / jacob_[k]};
    }

    double TriangleElement::map_x(double &zeta, double &eta, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < baseFunc.size(); i++) {
            switch (diffFlag) {
                case 0:
                    sum += globalVector(i).x * baseFunc[i](zeta, eta);
                    break;
                case 1:
                    sum += globalVector(i).x * diffQuotient_x(baseFunc[i], zeta, eta);
                    break;
                case 2:
                    sum += globalVector(i).x * diffQuotient_y(baseFunc[i], zeta, eta);
                    break;
            }
        }
        return sum;
    }

    double TriangleElement::map_y(double &zeta, double &eta, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < baseFunc.size(); i++) {
            switch (diffFlag) {
                case 0:
                    sum += globalVector(i).y * baseFunc[i](zeta, eta);
                    break;
                case 1:
                    sum += globalVector(i).y * diffQuotient_x(baseFunc[i], zeta, eta);
                    break;
                case 2:
                    sum += globalVector(i).y * diffQuotient_y(baseFunc[i], zeta, eta);
                    break;
            }
        }
        return sum;
    }

    double TriangleElement::map_x(int k, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < baseFunc.size(); i++) {
            switch (diffFlag) {
                case 0:
                    sum += globalVector(i).x * baseFunc_Values[i][k];
                    break;
                case 1:
                    sum += globalVector(i).x * baseFuncDiffQuotient_xValues[i][k];
                    break;
                case 2:
                    sum += globalVector(i).x * baseFuncDiffQuotient_yValues[i][k];
                    break;
            }
        }
        return sum;
    }

    double TriangleElement::map_y(int k, int diffFlag) {
        double sum = 0;
        for (int i = 0; i < baseFunc.size(); i++) {
            switch (diffFlag) {
                case 0:
                    sum += globalVector(i).y * baseFunc_Values[i][k];
                    break;
                case 1:
                    sum += globalVector(i).y * baseFuncDiffQuotient_xValues[i][k];
                    break;
                case 2:
                    sum += globalVector(i).y *  baseFuncDiffQuotient_yValues[i][k];
                    break;
            }
        }
        return sum;
    }

    Vertex2D& TriangleElement::globalVector(int &i){
        return vertices_[globalVectorIdx.indices[i]];
    }
//     std::complex<double> mes::air;
//     std::complex<double> mes::brick;
//     std::complex<double> mes::concrete;
//     std::complex<double> mes::damp;
//     std::complex<double> mes::steel;
//     std::complex<double> mes::glass;

    void setRefIdx(double f){
        mes::air = std::complex<double>(1.,0.);

        {
            std::vector<double> coeff ={3.91, 0.0238, 0.16};
            mes::brick = std::complex<double>(coeff[0], 17.98 * coeff[1] * pow(f, coeff[2] - 1.));
        }
        {
            std::vector<double> coeff ={5.24, 0.0462, 0.7822};
            mes::concrete = std::complex<double>(coeff[0], 17.98 * coeff[1] * pow(f, coeff[2] - 1.));
//            mes::damp = std::complex<double>(coeff[0], 17.98 * coeff[1] * pow(f, coeff[2] - 1.));
        }
        {
            std::vector<double> coeff ={1., pow(10.,7.), 0.};
            mes::steel = std::complex<double>(coeff[0], 17.98 * coeff[1] * pow(f, coeff[2] - 1.));
        }
        {
            std::vector<double> coeff ={6.31, .0036, 1.3394};
            mes::glass = std::complex<double>(coeff[0], 17.98 * coeff[1] * pow(f, coeff[2] - 1.));
        }
        {
            mes::concrete = mes::concrete;
        }
    }


    std::complex<double> TriangleElement::getRefIdx(elementType &_et) {
        switch (_et) {
            case AIR:
                return air;
            case BRICK:
                return brick;
            case CONCRETE:
                return concrete;
            case DAMP:
                return damp;
            case STEEL:
                return steel;
            case GLASS:
                return glass;
            default:
                return air;
        }
    }


}