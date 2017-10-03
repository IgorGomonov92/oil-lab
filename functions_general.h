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


void Construct_E(VectorXd * E);
void Construct_v(VectorXd * v );
void Construct_lamda(VectorXd * lamda ,VectorXd * E,VectorXd * v);
void Construct_G(VectorXd * G, VectorXd * E,VectorXd * v  );


#endif //FDM_FUNCTIONS_GENERAL_H
