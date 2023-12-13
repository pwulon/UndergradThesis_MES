//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_FORTUNES_HPP
#define MES_FORTUNES_HPP

#include <queue>
#include "../Elements/Point2D.hpp"
#include "Datastruct/Beachline.hpp"
#include "Datastruct/Event.hpp"
#include "Math/Parabola.hpp"
#include "Math/Circle.hpp"


//namespace bl = beachline;

bool arePointsCollinear(const std::vector<int> &idx,const std::vector<Point2D> &points);

double calculateSlope(const Point2D &p1,const Point2D &p2);

void build(const std::vector<Point2D> &points,
           std::vector<std::vector<int>> &elements);

#endif //MES_FORTUNES_HPP
