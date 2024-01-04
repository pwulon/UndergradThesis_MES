//
// Created by Pawulon on 26/12/2023.
//

#ifndef MES_WALL_HPP
#define MES_WALL_HPP

#include "Vertex2D.hpp"
#include "TriangleElement.hpp"

namespace mes{
class Wall {
private:
    Wall(double leftDownX, double leftDownY, double w, double h, mes::elementType elt);
    Vertex2D leftDownCorner, leftUpCorner, rightDownCorner, rightUpCorner;
    bool isRightTo(double &x, double &y, double &x_0, double &y_0, double &x_1, double &y_1);
public:

    static Wall createFromCorners(double leftDownX, double leftDownY, double RightUpX, double RightUpY, mes::elementType elt = mes::DAMP);
    static Wall createFromDimensions(double leftDownX, double leftDownY, double wi, double he, mes::elementType elt = mes::DAMP);

    mes::elementType type;

    Wall& rot(double ang);
    Wall& rotSelf(double ang);

    bool isInsideWall(double &x, double &y);

    bool isInsideWall(Vertex2D &p);


};

}

#endif //MES_WALL_HPP
