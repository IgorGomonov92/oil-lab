//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
#include "/home/igor/Eigen/Eigen/SparseCore"
#include "global.cpp"

using namespace Eigen;

void Construct_matrix_Laplace(SparseMatrix<double > * a);
void Construct_matrix_Poisson(SparseMatrix<double > * a1);
VectorXd * Solve_Laplace( );
void Construct_guess(VectorXd * initGuess);
void Construct_BC_Laplace(SparseVector<double> * bc);
void Construct_BC_Poisson(SparseVector<double> * bc);
void Construct_load_Laplace(double h, VectorXd * b, SparseVector<double> * bc);
void Construct_load_Poisson(double h, VectorXd * b, SparseVector<double> * bc, VectorXd * f);


#endif //FDM_FUNCTIONS_H
