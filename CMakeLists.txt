cmake_minimum_required(VERSION 3.22)
project(MES)

set(CMAKE_CXX_STANDARD 14)

add_compile_options(-O3 -fopenmp)

add_executable(MES main.cpp
        Elements/TriangleElement.cpp
        Elements/TriangleElement.hpp
        Elements/Vertex2D.cpp
        Elements/Vertex2D.hpp
        FortunesAlgo/Datastruct/Beachline.cpp
        FortunesAlgo/Datastruct/Beachline.hpp
        FortunesAlgo/Math/Circle.cpp
        FortunesAlgo/Math/Circle.hpp
        FortunesAlgo/Math/Parabola.cpp
        FortunesAlgo/Math/Parabola.hpp
        FortunesAlgo/Fortunes.cpp
        FortunesAlgo/Fortunes.hpp
        FortunesAlgo/Datastruct/Event.cpp
        FortunesAlgo/Datastruct/Event.hpp
        Elements/Wall.cpp
        Elements/Wall.hpp
        Elements/Section.cpp
        Elements/Section.hpp
        Visualization/glad/glad.h
        Visualization/glad/glad.c
        Visualization/glad/stb_image_write.h
        Visualization/Shaders.hpp
        Visualization/CreateOpenGlWindow.cpp
        Visualization/CreateOpenGlWindow.hpp
        Visualization/Shaders.cpp
        Solver/Solver.cpp
        Solver/Solver.hpp

)

target_link_libraries(MES glfw ${CMAKE_DL_LIBS})

