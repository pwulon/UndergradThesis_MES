//
// Created by Pawulon on 03/12/2023.
//

#include "Point2D.hpp"
#include <cmath>

#if defined(_WIN64) || defined(_WIN32)
#define isnan(x) _isnan(x)
#endif

using namespace std;


Point2D::Point2D(double _x, double _y) : x(_x), y(_y) {
}

Point2D::Point2D(const Point2D &point) : x(point.x), y(point.y) {
}

double dotProduct(const Point2D &p1, const Point2D &p2) {
    return p1.x * p2.x + p1.y * p2.y;
}

double crossProduct(const Point2D &p1, const Point2D &p2) {
    return p1.x * p2.y - p1.y * p2.x;
}

Point2D operator+(const Point2D &p1, const Point2D &p2) {
    return Point2D(p1.x + p2.x, p1.y + p2.y);
}

Point2D operator-(const Point2D &p1, const Point2D &p2) {
    return Point2D(p1.x - p2.x, p1.y - p2.y);
}

Point2D operator/(const Point2D &p1, const Point2D &p2) {
    return Point2D(p1.x / p2.x, p1.y / p2.y);
}

Point2D operator*(const Point2D &p, double value) {
    return Point2D(p.x * value, p.y * value);
}

Point2D operator*(double value, const Point2D &p) {
    return Point2D(p.x * value, p.y * value);
}

Point2D operator/(const Point2D &p, double value) {
    return Point2D(p.x / value, p.y / value);
}

Point2D operator-(const Point2D &p) {
    return Point2D(-p.x, -p.y);
}

std::ostream &operator<<(std::ostream &stream, const Point2D &p) {
    stream << "(" << p.x << "," << p.y << ")";
    return stream;
}

std::vector<Point2D> &operator<<(std::vector<Point2D> &v, const Point2D &p) {
    v.push_back(p);
    return v;
}

Point2D &Point2D::operator-=(const Point2D &p) {
    x -= p.x;
    y -= p.y;
    return *this;
}

Point2D &Point2D::operator+=(const Point2D &p) {
    x += p.x;
    y += p.y;
    return *this;
}

Point2D &Point2D::operator*=(double value) {
    x *= value;
    y *= value;
    return *this;
}

Point2D &Point2D::operator/=(double value) {
    x /= value;
    y /= value;
    return *this;
}

double Point2D::operator[](int i) {
    if (i==0) return x;
    else return y;
}

Point2D Point2D::normalized() {
    return (*this) / this->norm();
}


double Point2D::norm() {
    return sqrt(x * x + y * y);
}
