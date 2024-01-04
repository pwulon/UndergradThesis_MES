//
// Created by Pawulon on 26/12/2023.
//

#include "../Elements/Wall.hpp"


namespace mes {
    Wall::Wall(double leftDownX, double leftDownY, double w, double h, mes::elementType elt) :
            leftDownCorner{leftDownX, leftDownY},
            leftUpCorner{leftDownX, leftDownY + h},
            rightUpCorner{leftDownX + w, leftDownY + h},
            rightDownCorner(leftDownX + w, leftDownY),
            type{elt} {}

    bool Wall::isInsideWall(double &x, double &y) {

        return isRightTo(x, y, leftDownCorner.x, leftDownCorner.y, leftUpCorner.x, leftUpCorner.y) &&
               isRightTo(x, y, leftUpCorner.x, leftUpCorner.y, rightUpCorner.x, rightUpCorner.y) &&
               isRightTo(x, y, rightUpCorner.x, rightUpCorner.y, rightDownCorner.x, rightDownCorner.y) &&
               isRightTo(x, y, rightDownCorner.x, rightDownCorner.y, leftDownCorner.x, leftDownCorner.y);
    }

    Wall& Wall::rot(double ang) {
        double alfa = ang * M_PI / 180.;
        double sin1 = sin(alfa);
        double cos1 = cos(alfa);

        double x = leftDownCorner.x * cos1 - leftDownCorner.y * sin1;
        double y = leftDownCorner.x * sin1 + leftDownCorner.y * cos1;
        leftDownCorner.x = x;
        leftDownCorner.y = y;

        x = leftUpCorner.x * cos1 - leftUpCorner.y * sin1;
        y = leftUpCorner.x * sin1 + leftUpCorner.y * cos1;
        leftUpCorner.x = x;
        leftUpCorner.y = y;

        x = rightUpCorner.x * cos1 - rightUpCorner.y * sin1;
        y = rightUpCorner.x * sin1 + rightUpCorner.y * cos1;
        rightUpCorner.x = x;
        rightUpCorner.y = y;

        x = rightDownCorner.x * cos1 - rightDownCorner.y * sin1;
        y = rightDownCorner.x * sin1 + rightDownCorner.y * cos1;
        rightDownCorner.x = x;
        rightDownCorner.y = y;

        return *this;
    }

    Wall& Wall::rotSelf(double ang){
        double cVectorx = (leftDownCorner.x +
                            leftUpCorner.x +
                            rightUpCorner.x +
                            rightDownCorner.x)*.25;

        double cVectory = (leftDownCorner.y +
                           leftUpCorner.y +
                           rightUpCorner.y +
                           rightDownCorner.y)*.25;

        leftDownCorner.x -= cVectorx;
        leftUpCorner.x -= cVectorx;
        rightUpCorner.x -= cVectorx;
        rightDownCorner.x -= cVectorx;

        leftDownCorner.y -= cVectory;
        leftUpCorner.y -= cVectory;
        rightUpCorner.y -= cVectory;
        rightDownCorner.y -= cVectory;

        rot(ang);

        leftDownCorner.x += cVectorx;
        leftUpCorner.x += cVectorx;
        rightUpCorner.x += cVectorx;
        rightDownCorner.x += cVectorx;

        leftDownCorner.y += cVectory;
        leftUpCorner.y += cVectory;
        rightUpCorner.y += cVectory;
        rightDownCorner.y += cVectory;

        return *this;
    }

    bool Wall::isInsideWall(Vertex2D &p) {
        return isInsideWall(p.x, p.y);
    }


    bool Wall::isRightTo(double &x, double &y, double &x_0, double &y_0, double &x_1, double &y_1) {
        return (y - y_0) * (x_1 - x_0) - (x - x_0) * (y_1 - y_0) <= 0;
    }

    Wall Wall::createFromCorners(double leftDownX, double leftDownY, double RightUpX, double RightUpY,
                                 mes::elementType elt) {

        return Wall(leftDownX, leftDownY, RightUpX - leftDownX, RightUpY - leftDownY, elt);
    }

    Wall Wall::createFromDimensions(double leftDownX, double leftDownY, double wi, double he,
                                    mes::elementType elt) {
        return Wall( leftDownX, leftDownY, wi, he, elt);
    }
}