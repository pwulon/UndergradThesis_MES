//
// Created by Pawulon on 26/12/2023.
//

#ifndef MES_WALL_HPP
#define MES_WALL_HPP

#include "Vertex2D.hpp"
#include "TriangleElement.hpp"


class Wall {
public:
    mes::elementType type;
    Vertex2D leftDownCorner, leftUpCorner, rightDownCorner, rightUpCorner;

    Wall(double leftDownX, double leftDownY, double w, double h, mes::elementType elt);

    void rot(double ang);

    bool isInsideWall(double &x, double &y);

    bool isInsideWall(Vertex2D &p);

    bool isRightTo(double &x, double &y, double &x_0, double &y_0, double &x_1, double &y_1);
};


#endif //MES_WALL_HPP
