#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
using namespace boost::numeric::ublas;


int main(int argc, char **argv)
{
    int q=3; // величина дискретизации
    int n=q*q*q;
    vector<double> u(n); //неизв вектор
    double h=1.0; // шаг сетки
    int row = 0;



    matrix<double> a(n,n);
    vector<double> bc(n);
    vector<double> b(n);

    bc = Construct_BC(q);
    b = Construct_load(q, h, bc);
    a = Construct_matrix(q);



    u = BiCGSTAB(a,b);
    Print_matrix(a);
    Print_unknowns(u);



    return 0;
}
