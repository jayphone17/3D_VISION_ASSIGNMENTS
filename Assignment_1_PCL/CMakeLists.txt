cmake_minimum_required(VERSION 3.2)
#cmake 要求的最低版本

set(Qt5_DIR "/Users/jayphone/Qt/5.15.0/clang_64/lib/cmake/Qt5")

project(pcl_test)
#指定项目名称

#set(Boost_LIBRARY_DIR /usr/local/Cellar/boost/1.76.0/lib)

#设置boost包路径。

find_package(PCL 1.7 REQUIRED) 

#找到PCL包

list (REMOVE_ITEM PCL_LIBRARIES "vtkproj4")


include_directories(${PCL_INCLUDE_DIRS})

#宏，让cmake知道我们包含的外部头文件。
#其中PCL_INCLUDE_DIRS包含了我们需要的内容。

link_directories(${PCL_LIBRARY_DIRS})

add_definitions(${PCL_DEFINITIONS})

add_executable(pcl_test pcl_test.cpp)

target_link_libraries (pcl_test ${PCL_LIBRARIES})
#除了包含的头，还要让项目知道包含的库。

install(TARGETS pcl_test RUNTIME DESTINATION bin)
