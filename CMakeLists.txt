cmake_minimum_required(VERSION 3.8)
project(FDM)

include_directories(/home/igor/Eigen/Eigen)

set(CMAKE_CXX_STANDARD 17)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

set(SOURCE_FILES main.cpp functions_Laplace.h functions_Laplace.cpp global.cpp functions_Poisson.cpp functions_Poisson.h)
add_executable(FDM ${SOURCE_FILES})