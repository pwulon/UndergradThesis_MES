//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_FORTUNES_HPP
#define MES_FORTUNES_HPP

#include <queue>
#include "../Elements/Vertex2D.hpp"
#include "Datastruct/Beachline.hpp"
#include "Datastruct/Event.hpp"
#include "Math/Parabola.hpp"
#include "Math/Circle.hpp"



bool arePointsCollinear(const std::vector<int> &idx,const std::vector<Vertex2D> &points);

double calculateSlope(const Vertex2D &p1, const Vertex2D &p2);

void pointsRot(std::vector<Vertex2D> &points, double ang);

void build(std::vector<Vertex2D> &points,
           std::vector<std::vector<int>> &elements, bool withPointRot = true);

#endif //MES_FORTUNES_HPP
