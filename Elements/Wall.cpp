//
// Created by Pawulon on 26/12/2023.
//

#include "Wall.hpp"


bool Wall::isInsideWall(double &x, double &y){

    return isRightTo(x, y, leftDownCorner.x, leftDownCorner.y, leftDownCorner.x, rightUpCorner.y) &&
            isRightTo(x, y, leftDownCorner.x, rightUpCorner.y, rightUpCorner.x, rightUpCorner.y) &&
            isRightTo(x, y, rightUpCorner.x, rightUpCorner.y, rightUpCorner.x, leftDownCorner.y) &&
            isRightTo(x, y, rightUpCorner.x, leftDownCorner.y, leftDownCorner.x, leftDownCorner.y);
}

bool Wall::isInsideWall(Point2D &p){
    return isInsideWall(p.x, p.y);
}

bool Wall::isInsideWall(Vertex2D &p){
    return isInsideWall(p.x_, p.y_);
}

bool Wall::isRightTo(double &x, double &y,double &x_0, double &y_0, double &x_1, double &y_1){
    return (y - y_0)*(x_1 - x_0) - (x - x_0)*(y_1 - y_0) < 0;
};