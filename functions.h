//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "/home/igor/Eigen/Eigen/SparseCore"


using namespace Eigen;

VectorXd BiCGSTAB(SparseMatrix<double> a, VectorXd b);
void Print_matrix(SparseMatrix<double> a);
SparseMatrix<double> Construct_matrix_Laplace();
SparseMatrix<double> Construct_matrix_Poisson();
VectorXd Construct_BC_Laplace();
VectorXd Construct_BC_Poisson();
VectorXd Construct_load_Laplace(double d, VectorXd vector1);
VectorXd Construct_load_Poisson(double f, VectorXd vector1, VectorXd vector2);
void Print_vectors(VectorXd vector1);

#endif //FDM_FUNCTIONS_H
