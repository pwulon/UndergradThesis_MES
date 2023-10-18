#include <iostream>
#include <math.h>
#include <vector>
#include <array>
#include <functional>
#include <fstream>
#include <cstdio>

typedef std::array<std::pair<double, double>, 3> v3;


// DEBUG HELPER FUCTIONS
void printVerticesOfElement(const v3 &e) {
    std::cout << "0: " << e[0].first << " " << e[0].second <<
              ", 1: " << e[1].first << " " << e[1].second <<
              ", 2: " << e[2].first << " " << e[2].second << std::endl;
}

void printPair(const std::pair<double, double> &v) {
    std::cout << v.first << " " << v.second << std::endl;
}

// FUNKCJE KSZTALTU
double phi0(double zeta, double eta) {
    return -.5 * (zeta + eta);
}

double phi1(double zeta, double eta) {
    return .5 * (1 + zeta);
}

double phi2(double zeta, double eta) {
    return .5 * (1 + eta);
}

double rho(double x, double y){
    return exp(-.5*(pow(x,2)+pow(y,2)));
}

class TriEle {
private:
    // pochodne
    double diffQuotient_x(const std::function<double(double, double)> phi, double x, double y) {
        double delta = 1e-10;
        return (phi(x + delta, y) - phi(x - delta, y)) / (2 * delta);
    }

    double diffQuotient_y(const std::function<double(double, double)> phi, double x, double y) {
        double delta = 1e-10;
        return (phi(x, y + delta) - phi(x, y - delta)) / (2 * delta);
    }

//    //    attempts below dont work
//    double map_x(double zeta, double eta) {
//        double sum = 0;
//        for (int i = 0; i < 3; i++) {
//            sum += origVer[i].first * phi_k[i](zeta, eta);
//        }
//        return sum;
//    }
//    double map_y(double zeta, double eta) {
//        double sum = 0;
//        for (int i = 0; i < 3; i++) {
//            sum += origVer[i].second * phi_k[i](zeta, eta);
//        }
//        return sum;
//    }
//    std::function<double(double, double)> map_x ;
//    std::function<double(double, double)> map_y ;
public:
    std::pair<double, double> nablaPhik(double zeta, double eta, int k) {
        auto map_x = [this](double zeta, double eta) {
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += origVer[i].first * phi_k[i](zeta, eta);
            }
            return sum;
        };
        auto map_y = [this](double zeta, double eta) {
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += origVer[i].second * phi_k[i](zeta, eta);
            }
            return sum;
        };
        return {diffQuotient_x(phi_k[k], zeta, eta) *
                diffQuotient_y(map_y, zeta, eta) / Jacobian(zeta, eta) -
                diffQuotient_y(phi_k[k], zeta, eta) *
                diffQuotient_x(map_y, zeta, eta) / Jacobian(zeta, eta),
                -diffQuotient_x(phi_k[k], zeta, eta) *
                diffQuotient_y(map_x, zeta, eta) / Jacobian(zeta, eta) +
                diffQuotient_y(phi_k[k], zeta, eta) *
                diffQuotient_x(map_x, zeta, eta) / Jacobian(zeta, eta)};
    }


    double Jacobian(double zeta, double eta) {
        auto map_x = [this](double zeta, double eta) {
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += origVer[i].first * phi_k[i](zeta, eta);
            }
            return sum;
        };
        auto map_y = [this](double zeta, double eta) {
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += origVer[i].second * phi_k[i](zeta, eta);
            }
            return sum;
        };

        return diffQuotient_x(map_x, zeta, eta) *
               diffQuotient_y(map_y, zeta, eta) -
               diffQuotient_x(map_y, zeta, eta) *
               diffQuotient_y(map_x, zeta, eta);
    }

    void initE() {

        double xgaus[] = {-0.3333333, -0.0597158717, -0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853,
                          0.5948539707};
        double ygaus[] = {-0.3333333, -0.0597158717, -0.8805682564, -0.0597158717, -0.7974269853, 0.5948539707,
                          -0.7974269853};
        double wgaus[] = {0.45, 0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double sum = 0;
                for (int k = 0; k < 7; k++) {
                    std::pair<double, double> nablaphi_i = nablaPhik(xgaus[k], ygaus[k], i);
                    std::pair<double, double> nablaphi_j = nablaPhik(xgaus[k], ygaus[k], j);
                    sum += wgaus[k] * Jacobian(xgaus[k], ygaus[k])
                           * (nablaphi_i.first * nablaphi_j.first
                              + nablaphi_i.second * nablaphi_j.second);
                }
                E[i][j] = sum;
                if(m_<2)std::cout<<sum<<"\t"; //debug
            }
            if(m_<2)std::cout<<std::endl; //debug
        }
        if(m_<2)std::cout<<std::endl; //debug
    }
    void initF(){
        auto map_x = [this](double zeta, double eta) {
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += origVer[i].first * phi_k[i](zeta, eta);
            }
            return sum;
        };
        auto map_y = [this](double zeta, double eta) {
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                sum += origVer[i].second * phi_k[i](zeta, eta);
            }
            return sum;
        };
        double xgaus[] = {-0.3333333, -0.0597158717, -0.0597158717, -0.8805682564, -0.7974269853, -0.7974269853,
                          0.5948539707};
        double ygaus[] = {-0.3333333, -0.0597158717, -0.8805682564, -0.0597158717, -0.7974269853, 0.5948539707,
                          -0.7974269853};
        double wgaus[] = {0.45, 0.2647883055, 0.2647883055, 0.2647883055, 0.251878361, 0.251878361, 0.251878361};

        for(int j =0 ;j<3;j++){
            double sum = 0;
            for (int k = 0; k < 7; k++){
                sum += wgaus[k] * Jacobian(xgaus[k], ygaus[k]) *
                        phi_k[j](xgaus[k], ygaus[k]) *
                        rho(map_x(xgaus[k], ygaus[k]),map_y(xgaus[k], ygaus[k]));
            }
            F[j] = sum;
            if(m_<2)std::cout<<sum<<"\t"; //debug
        }
        if(m_<2)std::cout<<std::endl<<std::endl; //debug
    }

    int m_;
    v3 origVer;
    v3 locVer;
    std::array<int, 3> localLabel;
    std::array<std::array<double, 3>, 3> E;
    std::array<double, 3> F;
    std::array<std::function<double(double, double)>, 3> phi_k;

    TriEle(int m, v3 orig, std::array<int, 3> loc) : m_{m}, localLabel{loc} {
        phi_k[0] = &phi0;
        phi_k[1] = &phi1;
        phi_k[2] = &phi2;
        locVer[0] = {-1., -1.};
        locVer[1] = {1., -1.};
        locVer[2] = {-1., 1.};
        origVer[0] = orig[0];
        origVer[1] = orig[1];
        origVer[2] = orig[2];

//        map_x = [this](double zeta, double eta) {
//            double sum = 0;
//            for (int i = 0; i < 3; i++) {
//                sum += origVer[i].first * phi_k[i](zeta, eta);
//            }
//            return sum;
//        };
//        map_y = [this](double zeta, double eta) {
//            double sum = 0;
//            for (int i = 0; i < 3; i++) {
//                sum += origVer[i].second * phi_k[i](zeta, eta);
//            }
//            return sum;
//        };

        initE();
        initF();
    }


};

int main() {
    const double L = 10;
    const size_t N = 100; //vertices number
    const size_t row = 10; //vertices number

    double horLen = L / (L - 1);

    std::array<std::pair<double, double>, N> verticies{};
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++)
            verticies[i * 10 + j] = {-5.0 + horLen * i, -5.0 + horLen * j};
    }


    const size_t nEle = (row - 1) * (row - 1) * 2;
    std::vector<TriEle> Elements;

    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            int e = 2 * (i * 9 + j);
            std::array<int, 3> iVer = {i * 10 + j, (i + 1) * 10 + j, i * 10 + j + 1};
            v3 ver = {verticies[iVer[0]],
                       verticies[iVer[1]],
                       verticies[iVer[2]]};

            Elements.emplace_back(e, ver, iVer);
            iVer = {(i + 1) * 10 + j, (i + 1) * 10 + j + 1, i * 10 + j + 1};
            ver = {verticies[iVer[0]],
                    verticies[iVer[1]],
                    verticies[iVer[2]]};

            Elements.emplace_back(e + 1, ver, iVer);
        }
    }
    Elements[1].initE();


    std::array<std::array<double, 100>,100> S{};

    for(int m = 0; m<nEle; m++){
        std::array<std::array<double,3>,3>&E = Elements[m].E;
        std::array<int,3>&glob = Elements[m].localLabel;
        for (int k=0;k<3;k++){
            for(int l = 0; l<3; l++ ){
                int i = glob[k];
                int j = glob[l];
                S[i][j] = S[i][j] + E[k][l];
            }
        }
    }
    FILE* file;
    file = fopen("S.txt","w");
    for(int i = 0; i<100; i++){
        for(int j = 0; j<100; j++) fprintf(file,"%.2f\t",S[i][j]);
        fprintf(file,"\n");
    }
    fclose(file);

    std::array<double,100> F{};
    for(int m = 0; m<nEle; m++){
        std::array<double, 3>&Fm = Elements[m].F;
        std::array<int,3>&glob = Elements[m].localLabel;
        for (int k=0;k<3;k++){
                int i = glob[k];
                F[i] = F[i] + Fm[k];
        }
    }
    file = fopen("F.txt","w");
    for(int j = 0; j<100; j++) fprintf(file,"%g\n",F[j]);
    fclose(file);

    return 0;
}
