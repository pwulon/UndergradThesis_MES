//
// Created by Pawulon on 28/12/2023.
//

#ifndef MES_SECTION_HPP
#define MES_SECTION_HPP

#include <iostream>
#include <vector>
#include <unordered_set>
#include "Vertex2D.hpp"

struct Section;

struct SecComp;

struct SecHash ;
void addQuadVertElements(std::vector<Vertex2D> &ver, std::vector<std::vector<int>> elementsIdx);

#endif //MES_SECTION_HPP
