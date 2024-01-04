//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_PARABOLA_HPP
#define MES_PARABOLA_HPP

#include "../../Elements/Vertex2D.hpp"

namespace mes::fortunes{/**

 Calculate number of intersection points between two parabolas with foci `f1` and `f2` and with given `directrix`

 */
int intersectionPointsNum(const Vertex2D &f1, const Vertex2D &f2, double directrix);


/**

 Find intersection points of two parabolas with foci `f1` and `f2` and with given `directrix`

 */
std::vector<Vertex2D> findIntersectionPoints(const Vertex2D &f1, const Vertex2D &f2, double directrix);

}

#endif //MES_PARABOLA_HPP
