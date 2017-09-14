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
matrix<double> Construct_matrix_Laplace();
matrix<double> Construct_matrix_Poisson();
vector<double> Construct_BC_Laplace();
vector<double> Construct_BC_Poisson();
vector<double> Construct_load_Laplace(double d, vector<double> vector1);
vector<double> Construct_load_Poisson(double f, vector<double> vector1, vector<double> vector2);
void Print_vectors(vector<double> vector1);

#endif //FDM_FUNCTIONS_H
