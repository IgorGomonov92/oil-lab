//
// Created by igor on 21.09.17.
//

#ifndef FDM_FUNCTIONS_POISSON_H
#define FDM_FUNCTIONS_POISSON_H
#include <iostream>
#include <omp.h>
#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include <chrono>
#include "functions_general.h"

using namespace Eigen;

void Construct_matrix_Poisson(SparseMatrix<double > * A1);
void Construct_guess_P(VectorXd * initGuess);
void Construct_f( VectorXd * f, VectorXd * uL);
void Construct_BC_Poisson(VectorXd * bc);
void Construct_load_Poisson(VectorXd * b, VectorXd * bc, VectorXd * f);
VectorXd  Solve_Poisson(VectorXd * uL);
void Construct_w_Derivative_z(  VectorXd * w_Derivative_z , VectorXd *  w);
void Construct_w_Derivative_x(  VectorXd * w_Derivative_x , VectorXd *  w);
void Construct_BC_Poisson_XX(VectorXd * bcP);
void Construct_load_Poisson_XX(VectorXd * bP, VectorXd * bcP, VectorXd * f);

#endif //FDM_FUNCTIONS_POISSON_H
