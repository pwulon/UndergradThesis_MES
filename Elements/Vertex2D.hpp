//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_VERTEX2D_HPP
#define MES_VERTEX2D_HPP

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>


#define POINT_EPSILON 1.0e-6

namespace mes {




    class Vertex2D {

    public:
        double x, y;
        bool isBorder;


        Vertex2D(double x = 0.0, double y = 0.0, bool _b = false);

        Vertex2D(const Vertex2D &point);

//        void setElementType(elementType &_et);

        friend double dotProduct(const Vertex2D &p1, const Vertex2D &p2);

        friend double crossProduct(const Vertex2D &p1, const Vertex2D &p2);

        friend Vertex2D operator+(const Vertex2D &p1, const Vertex2D &p2);

        friend Vertex2D operator-(const Vertex2D &p1, const Vertex2D &p2);

        friend Vertex2D operator/(const Vertex2D &p1, const Vertex2D &p2);

        friend Vertex2D operator*(const Vertex2D &p, double value);

        friend Vertex2D operator*(double value, const Vertex2D &p);

        friend Vertex2D operator/(const Vertex2D &p, double value);

        friend Vertex2D operator-(const Vertex2D &p);

        friend std::ostream &operator<<(std::ostream &stream, const Vertex2D &p);

        friend std::vector<Vertex2D> &operator<<(std::vector<Vertex2D> &v, const Vertex2D &p);

        Vertex2D &operator-=(const Vertex2D &p);

        Vertex2D &operator+=(const Vertex2D &p);

        Vertex2D &operator*=(double value);

        Vertex2D &operator/=(double value);

        Vertex2D normalized();

        double norm();

        double operator[](int i);

    };

    double dotProduct(const Vertex2D &p1, const Vertex2D &p2);

    double crossProduct(const Vertex2D &p1, const Vertex2D &p2);
}

#endif //MES_VERTEX2D_HPP
