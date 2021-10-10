cmake_minimum_required(VERSION 2.8)

project(learnCGAL)

find_package(CGAL COMPONENTS Qt5 Core)
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

include(${CGAL_USE_FILE})


add_executable(demo_10 demo_10.cpp)

target_link_libraries(demo_10 ${CGAL_LIBS})

