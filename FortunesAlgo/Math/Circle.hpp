//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_CIRCLE_HPP
#define MES_CIRCLE_HPP

#include "../../Elements/Vertex2D.hpp"

#define CIRCLE_CENTER_EPSILON 1.0e-7

/**

 Find a center of a circle with given three points.
 Returns false if points are collinear.
 Otherwise returns true and updates x- and y-coordinates of the `center` of circle.

 */
namespace mes::fortunes{
    bool findCircleCenter(const Vertex2D &p1, const Vertex2D &p2, const Vertex2D &p3, Vertex2D &center);
}



#endif //MES_CIRCLE_HPP
