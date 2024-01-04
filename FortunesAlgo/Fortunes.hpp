//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_FORTUNES_HPP
#define MES_FORTUNES_HPP

#include <queue>
#include "../Elements/Vertex2D.hpp"
#include "../Elements/Wall.hpp"
#include "Datastruct/Beachline.hpp"
#include "Datastruct/Event.hpp"
#include "Math/Parabola.hpp"
#include "Math/Circle.hpp"

namespace mes{
    namespace fortunes{
        void pointsRot(std::vector<Vertex2D> &points, double ang);

        void build(std::vector<Vertex2D> &points,
                   std::vector<mes::ElementIndices> &elements,
                   std::vector<Wall> walls, mes::baseFuncType _bft = mes::LIN, bool withPointRot = true);

    }
}

#endif //MES_FORTUNES_HPP
