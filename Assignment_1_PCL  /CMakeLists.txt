cmake_minimum_required(VERSION 3.2)
project(pcl_test)

set(Boost_LIBRARY_DIR /usr/lib/x86_64-linux-gnu/)

find_package(PCL 1.7 REQUIRED)

list (REMOVE_ITEM PCL_LIBRARIES "vtkproj4")
#set(PCL_DIR /usr/lib/x86_64-linux-gnu/cmake/pcl)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(pcl_test pcl_test.cpp)

target_link_libraries (pcl_test ${PCL_LIBRARIES})

install(TARGETS pcl_test RUNTIME DESTINATION bin)

