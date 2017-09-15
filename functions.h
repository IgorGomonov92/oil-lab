//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "/home/igor/Eigen/Eigen/SparseCore"
#include "global.cpp"

using namespace Eigen;

SpMat BiCGSTAB(SpMat a, SpMat b);
//void Print_matrix(SparseMatrix<double> a);
SparseMatrix<double>   Construct_matrix_Laplace();
SparseMatrix<double>   Construct_matrix_Poisson();

/*SpMat Construct_BC_Laplace();
SpMat Construct_BC_Poisson();
SpMat Construct_load_Laplace(double d, SpMat vector1);
SpMat Construct_load_Poisson(double f, SpMat vector1, SpMat vector2);
void Print_vectors(SpMat vector1);
*/
#endif //FDM_FUNCTIONS_H
