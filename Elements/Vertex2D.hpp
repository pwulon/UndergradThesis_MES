//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_VERTEX2D_HPP
#define MES_VERTEX2D_HPP

#include "Point2D.hpp"
#include <memory>

class Vertex2D{
public:
    int id_;

    double x_, y_;
    bool isBorder;

    Vertex2D(Point2D &pos, bool b = false): x_{pos.x}, y_{pos.y}, isBorder{b}{};

    double x() { return this->x_; }
    double y() { return this->y_; }

};



#endif //MES_VERTEX2D_HPP
