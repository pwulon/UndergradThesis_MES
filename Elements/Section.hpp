//
// Created by Pawulon on 28/12/2023.
//

#ifndef MES_SECTION_HPP
#define MES_SECTION_HPP

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "Vertex2D.hpp"
#include "TriangleElement.hpp"

namespace mes{
    struct Section;

    struct SecComp;

    struct SecHash;

    void addQuadVertElements(std::vector<Vertex2D> &ver, std::vector<mes::ElementIndices> &elementsIdx);

}
#endif //MES_SECTION_HPP
