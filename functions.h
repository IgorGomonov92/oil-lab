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
matrix<double> Construct_matrix(int q);
vector<double> Construct_BC(int q);
vector<double> Construct_load(int q, double h, vector<double> bc);
void Print_unknowns( vector<double> u);

#endif //FDM_FUNCTIONS_H
