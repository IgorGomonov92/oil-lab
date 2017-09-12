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
matrix<double> Construct_matrix_Laplace(int q);
matrix<double> Construct_matrix_Poisson(int q);
vector<double> Construct_BC_Laplace(int q);
vector<double> Construct_BC_Poisson(int q);
vector<double> Construct_load_Laplace(int q, double h, vector<double> bc);
vector<double> Construct_load_Poisson(int q, double h, vector<double> bc, vector<double> f);
void Print_vectors(int q, vector<double> u);

#endif //FDM_FUNCTIONS_H
