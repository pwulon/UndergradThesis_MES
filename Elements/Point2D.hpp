//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_POINT2D_HPP
#define MES_POINT2D_HPP

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

#define POINT_EPSILON 1.0e-6

class Point2D {

public:
    double x, y;


    Point2D(double x = 0.0, double y = 0.0);
    Point2D(const Point2D &point);

    friend double dotProduct(const Point2D &p1, const Point2D &p2);
    friend double crossProduct(const Point2D &p1, const Point2D &p2);

    friend Point2D operator+(const Point2D &p1, const Point2D &p2);
    friend Point2D operator-(const Point2D &p1, const Point2D &p2);
    friend Point2D operator/(const Point2D &p1, const Point2D &p2);
    friend Point2D operator*(const Point2D &p, double value);
    friend Point2D operator*(double value, const Point2D &p);
    friend Point2D operator/(const Point2D &p, double value);
    friend Point2D operator-(const Point2D &p);

    friend std::ostream &operator<<(std::ostream &stream, const Point2D &p);
    friend std::vector<Point2D> &operator<<(std::vector<Point2D> &v, const Point2D &p);

    Point2D &operator-=(const Point2D &p);
    Point2D &operator+=(const Point2D &p);
    Point2D &operator*=(double value);
    Point2D &operator/=(double value);

    Point2D normalized();
    double norm();

    double operator[](int i);



};

double dotProduct(const Point2D &p1, const Point2D &p2);
double crossProduct(const Point2D &p1, const Point2D &p2);



#endif //MES_POINT2D_HPP
