//
// Created by Pawulon on 30/12/2023.
//

#ifndef MES_CREATEOPENGLWINDOW_HPP
#define MES_CREATEOPENGLWINDOW_HPP


#include "Shaders.hpp"
#include "../Elements/TriangleElement.hpp"
#include "../Elements/Vertex2D.hpp"
//#include <eigen3/Eigen/Eigenvalues>




void saveScreenshot(const char* filename, GLFWwindow* window);

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

int CreateOpenGlWindow(std::vector<Vertex2D> &vertices, std::vector<fem::TriangleElement> &Elements, std::vector<double> &c, double L);


#endif //MES_CREATEOPENGLWINDOW_HPP
