//
// Created by Pawulon on 03/12/2023.
//

#include "../Elements/Vertex2D.hpp"


#if defined(_WIN64) || defined(_WIN32)
#define isnan(x) _isnan(x)
#endif

using namespace std;

namespace mes{
    Vertex2D::Vertex2D(double _x, double _y, bool _b) : x(_x), y(_y), isBorder{_b} {
    }

    Vertex2D::Vertex2D(const Vertex2D &point) : x(point.x), y(point.y), isBorder{point.isBorder} {
    }


    double dotProduct(const Vertex2D &p1, const Vertex2D &p2) {
        return p1.x * p2.x + p1.y * p2.y;
    }

    double crossProduct(const Vertex2D &p1, const Vertex2D &p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }

    Vertex2D operator+(const Vertex2D &p1, const Vertex2D &p2) {
        return Vertex2D(p1.x + p2.x, p1.y + p2.y);
    }

    Vertex2D operator-(const Vertex2D &p1, const Vertex2D &p2) {
        return Vertex2D(p1.x - p2.x, p1.y - p2.y);
    }

    Vertex2D operator/(const Vertex2D &p1, const Vertex2D &p2) {
        return Vertex2D(p1.x / p2.x, p1.y / p2.y);
    }

    Vertex2D operator*(const Vertex2D &p, double value) {
        return Vertex2D(p.x * value, p.y * value);
    }

    Vertex2D operator*(double value, const Vertex2D &p) {
        return Vertex2D(p.x * value, p.y * value);
    }

    Vertex2D operator/(const Vertex2D &p, double value) {
        return Vertex2D(p.x / value, p.y / value);
    }

    Vertex2D operator-(const Vertex2D &p) {
        return Vertex2D(-p.x, -p.y);
    }

    std::ostream &operator<<(std::ostream &stream, const Vertex2D &p) {
        stream << "(" << p.x << "," << p.y << ")";
        return stream;
    }

    std::vector<Vertex2D> &operator<<(std::vector<Vertex2D> &v, const Vertex2D &p) {
        v.push_back(p);
        return v;
    }

    Vertex2D &Vertex2D::operator-=(const Vertex2D &p) {
        x -= p.x;
        y -= p.y;
        return *this;
    }

    Vertex2D &Vertex2D::operator+=(const Vertex2D &p) {
        x += p.x;
        y += p.y;
        return *this;
    }

    Vertex2D &Vertex2D::operator*=(double value) {
        x *= value;
        y *= value;
        return *this;
    }

    Vertex2D &Vertex2D::operator/=(double value) {
        x /= value;
        y /= value;
        return *this;
    }

    double Vertex2D::operator[](int i) {
        if (i == 0) return x;
        else return y;
    }

    Vertex2D Vertex2D::normalized() {
        return (*this) / this->norm();
    }


    double Vertex2D::norm() {
        return sqrt(x * x + y * y);
    }

    void Vertex2D::setElementType(elementType &_et) {
        et = _et;
    }

}