cmake_minimum_required(VERSION 2.8)

set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14" )

project(Shortest_path_on_mesh)

find_package(CGAL COMPONENTS Core)
include(${CGAL_USE_FILE})

add_executable(Shortest_path_on_mesh main.cpp)

target_link_libraries(Shortest_path_on_mesh ${CGAL_LIBS})

