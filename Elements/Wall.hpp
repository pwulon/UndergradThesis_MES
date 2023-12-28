//
// Created by Pawulon on 26/12/2023.
//

#ifndef MES_WALL_HPP
#define MES_WALL_HPP

#include "Point2D.hpp"
#include "Vertex2D.hpp"

class Wall {
public:
    Point2D leftDownCorner, rightUpCorner;
    Wall(Point2D &ldc, Point2D &ruc): leftDownCorner{ldc}, rightUpCorner(ruc){}

    bool isInsideWall(double &x, double &y);
    bool isInsideWall(Point2D &p);
    bool isInsideWall(Vertex2D &p);

    bool isRightTo(double &x, double &y,double &x_0, double &y_0, double &x_1, double &y_1);
};


#endif //MES_WALL_HPP
