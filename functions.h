//
// Created by gomonov on 11.09.17.
//

#ifndef FDM_FUNCTIONS_H
#define FDM_FUNCTIONS_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


using namespace boost::numeric::ublas;

vector<double> BiCGSTAB(matrix<double> a, vector<double> b);
void Print_matrix(matrix<double> a);
matrix<double> Construct_matrix_Laplace(int q, int i, int i1);
matrix<double> Construct_matrix_Poisson(int q, int i, int i1);
vector<double> Construct_BC_Laplace(int q);
vector<double> Construct_BC_Poisson(int q);
vector<double> Construct_load_Laplace(int q, int h, int bc, double d, vector<double> vector1);
vector<double> Construct_load_Poisson(int q, int h, int bc, double f, vector<double> vector1, vector<double> vector2);
void Print_vectors(int q, int u, int i, vector<double> vector1);

#endif //FDM_FUNCTIONS_H
