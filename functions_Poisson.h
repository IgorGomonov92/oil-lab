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

using namespace Eigen;

void Construct_matrix_Poisson(SparseMatrix<double > * A1);
void Construct_guess_P(VectorXd * initGuess);
void Construct_f( std::vector<VectorXd> * f, VectorXd * uL);
void Construct_BC_Poisson(VectorXd * bc);
void Construct_load_Poisson(VectorXd * b, VectorXd * bc, VectorXd * f);
std::vector<VectorXd>  Solve_Poissons(VectorXd * uL);

#endif //FDM_FUNCTIONS_POISSON_H
