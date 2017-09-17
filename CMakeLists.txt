cmake_minimum_required(VERSION 3.8)
project(FDM)

include_directories(/home/igor/Eigen/Eigen)

set(CMAKE_CXX_STANDARD 11)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES main.cpp functions.h functions.cpp global.cpp)
add_executable(FDM ${SOURCE_FILES})