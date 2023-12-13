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

    std::shared_ptr<Point2D> point;
    bool isBorder;

    Vertex2D(Point2D &pos, bool b = false): point{std::make_shared<Point2D>(pos)}, isBorder{b}{};

    double x() { return point->x; }
    double y() { return point->y; }

};



#endif //MES_VERTEX2D_HPP
