cmake_minimum_required(VERSION 3.22)
project(MES)

set(CMAKE_CXX_STANDARD 14)

add_executable(MES main.cpp
        Elements/TriangleElement.cpp
        Elements/TriangleElement.hpp
        Elements/Point2D.cpp
        Elements/Point2D.hpp
        FortunesAlgo/Datastruct/Beachline.cpp
        FortunesAlgo/Datastruct/Beachline.hpp
        FortunesAlgo/Math/Circle.cpp
        FortunesAlgo/Math/Circle.hpp
        FortunesAlgo/Math/Parabola.cpp
        FortunesAlgo/Math/Parabola.hpp
        FortunesAlgo/Fortunes.cpp
        FortunesAlgo/Fortunes.hpp
        Elements/Vertex2D.hpp
        FortunesAlgo/Datastruct/Event.cpp
        FortunesAlgo/Datastruct/Event.hpp
        Elements/Wall.cpp
        Elements/Wall.hpp)
