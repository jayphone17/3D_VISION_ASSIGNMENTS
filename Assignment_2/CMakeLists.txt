cmake_minimum_required(VERSION 2.8)

set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14" )

project(RegionGrowing)

find_package(CGAL COMPONENTS Core)
include(${CGAL_USE_FILE})

find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(RegionGrowing ${Boost_LIBRARY_DIRS})

add_executable(RegionGrowing RegionGrowing.cpp)
target_link_libraries(RegionGrowing ${CGAL_LIBS} ${Boost_LIBRARIES})