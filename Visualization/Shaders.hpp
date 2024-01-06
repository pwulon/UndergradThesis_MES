//
// Created by Pawulon on 30/12/2023.
//

#ifndef MES_SHADERS_HPP
#define MES_SHADERS_HPP

#include <iostream>
#include <vector>
#include <memory>

#include "glad/glad.h"
#include "../Elements/TriangleElement.hpp"
#include <GLFW/glfw3.h>


namespace mes::plot{
    unsigned int createShaderProgram(bool log = true);
    unsigned int createShaderProgramGeometry(bool log = true);

    std::shared_ptr<GLfloat> makeOrto(double l, double r, double b, double t, double f, double n);
}
#endif //MES_SHADERS_HPP
