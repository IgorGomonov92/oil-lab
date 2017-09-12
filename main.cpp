#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
using namespace boost::numeric::ublas;


int main(int argc, char **argv)
{
    unsigned long q; // величина дискретизации

    unsigned long n;
    vector<double> u(n); //неизв вектор
    double h=1.0; // шаг сетки

    q = 10;
    n = q * q * q;

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
