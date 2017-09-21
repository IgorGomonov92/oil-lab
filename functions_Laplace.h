//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
#include "/home/igor/Eigen/Eigen/SparseCore"

using namespace Eigen;

void Construct_matrix_Laplace(SparseMatrix<double > * A);
void Construct_guess_L(VectorXd * initGuess);
void Construct_BC_Laplace(SparseVector<double> * bc);
void Construct_load_Laplace(VectorXd * b, SparseVector<double> * bc);
VectorXd  Solve_Laplace();

#endif //FDM_FUNCTIONS_H
