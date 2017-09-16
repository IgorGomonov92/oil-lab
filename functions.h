//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include "/home/igor/Eigen/Eigen/SparseCore"
#include "global.cpp"

using namespace Eigen;

//SpMat BiCGSTAB(SpMat a, SpMat b);
//void Print_matrix(SparseMatrix<double> a);
void Construct_matrix_Laplace(SparseMatrix<double> * a);
void Construct_matrix_Poisson(SparseMatrix<double> * a1);

void Construct_BC_Laplace(SparseVector<double> * b);
void Construct_BC_Poisson(SparseVector<double> * b);
void Construct_load_Laplace(double d, SparseVector<double> * vec, SparseVector<double> * b);
void Construct_load_Poisson(double f, SparseVector<double> * vec, SparseVector<double> b,  SparseVector<double> vector1, SparseVector<double> vector2);
//void Print_vectors(SparseVector<double> vector1);

#endif //FDM_FUNCTIONS_H
