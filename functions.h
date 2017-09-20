//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
#include "/home/igor/Eigen/Eigen/SparseCore"
#include "global.cpp"

using namespace Eigen;

void Construct_matrix_Laplace(SparseMatrix<double > * A);
void Construct_matrix_Poisson(SparseMatrix<double > * A1);
void Construct_f( std::vector<VectorXd> f);
void Construct_guess_L(VectorXd * initGuess);
void Construct_guess_P(VectorXd * initGuess);
void Construct_BC_Laplace(SparseVector<double> * bc);
void Construct_BC_Poisson(VectorXd * bc);
void Construct_load_Laplace(VectorXd * b, SparseVector<double> * bc);
void Construct_load_Poisson(VectorXd * b, VectorXd * bc, std::vector<VectorXd> * f);

VectorXd  Solve_Laplace();
VectorXd  Solve_Poissons();
//VectorXd * Solve_Poissons(SparseMatrix<double > * AL, VectorXd * b, VectorXd * initGuess);



#endif //FDM_FUNCTIONS_H
