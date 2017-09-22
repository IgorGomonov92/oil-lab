//
// Created by igor on 22.09.17.
//

#ifndef FDM_FUNCTIONS_GENERAL_H
#define FDM_FUNCTIONS_GENERAL_H
#include <omp.h>
#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include <chrono>

using namespace Eigen;


void Construct_E(std::vector<double> * E);
void Construct_v(std::vector<double> * v );
void Construct_lamda(std::vector<double> * lamda ,std::vector<double> * E,std::vector<double> * v);
void Construct_G(std::vector<double> * G, std::vector<double> * E,std::vector<double> * v  );


#endif //FDM_FUNCTIONS_GENERAL_H
