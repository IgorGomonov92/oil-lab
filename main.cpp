#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <omp.h>
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

    omp_set_num_threads(8);  // кол во тредов

    matrix<double> a(n,n);
    vector<double> bc(n);
    vector<double> b(n);

    bc = Construct_BC(q);
    b = Construct_load(q, h, bc);
    a = Construct_matrix_Laplace(q);

//реализация с BicGstab c OMP
#pragma omp parallel
    {
        u = BiCGSTAB(a, b); // решаем уравнение Лапласа
    }
    Print_matrix(a);
    Print_unknowns(u);



    return 0;
}
