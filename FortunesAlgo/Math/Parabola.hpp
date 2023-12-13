//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_PARABOLA_HPP
#define MES_PARABOLA_HPP

#include "../../Elements/Point2D.hpp"
/**

 Calculate number of intersection points between two parabolas with foci `f1` and `f2` and with given `directrix`

 */
int intersectionPointsNum(const Point2D &f1, const Point2D &f2, double directrix);


/**

 Find intersection points of two parabolas with foci `f1` and `f2` and with given `directrix`

 */
std::vector<Point2D> findIntersectionPoints(const Point2D &f1, const Point2D &f2, double directrix);



#endif //MES_PARABOLA_HPP
