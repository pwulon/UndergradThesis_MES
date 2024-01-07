cmake_minimum_required(VERSION 3.22)
project(MES)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-O3 -fopenmp)

add_executable(MES main.cpp
        src/TriangleElement.cpp
        Elements/TriangleElement.hpp
        src/Vertex2D.cpp
        Elements/Vertex2D.hpp
        src/Beachline.cpp
        FortunesAlgo/Datastruct/Beachline.hpp
        src/Circle.cpp
        FortunesAlgo/Math/Circle.hpp
        src/Parabola.cpp
        FortunesAlgo/Math/Parabola.hpp
        src/Fortunes.cpp
        FortunesAlgo/Fortunes.hpp
        src/Event.cpp
        FortunesAlgo/Datastruct/Event.hpp
        src/Wall.cpp
        Elements/Wall.hpp
        src/Section.cpp
        Elements/Section.hpp
        Visualization/glad/glad.h
        src/glad.c
        Visualization/glad/stb_image_write.h
        Visualization/Shaders.hpp
        src/CreateOpenGlWindow.cpp
        Visualization/CreateOpenGlWindow.hpp
        src/Shaders.cpp
        src/Solver.cpp
        Solver/Solver.hpp
        Visualization/base64.h
        MES.h
)

target_link_libraries(MES glfw ${CMAKE_DL_LIBS})

