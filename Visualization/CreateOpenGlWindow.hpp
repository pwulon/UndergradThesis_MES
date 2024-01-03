//
// Created by Pawulon on 30/12/2023.
//

#ifndef MES_CREATEOPENGLWINDOW_HPP
#define MES_CREATEOPENGLWINDOW_HPP


#include "Shaders.hpp"
#include "../Elements/TriangleElement.hpp"
#include "../Elements/Vertex2D.hpp"
#include "base64.h"

#include <iostream>
#include <ctime>
#include <sstream>




void saveScreenshot(const char* filename, GLFWwindow* window);

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

std::string CreateOpenGlWindow(std::vector<Vertex2D> &vertices,
                       std::vector<fem::TriangleElement> &Elements,
                       std::vector<double> &c,
                       double &w, double &h,
                       int resolutionW, int resolutionH);


#endif //MES_CREATEOPENGLWINDOW_HPP
