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

    q = 6;
    n = q * q * q;

    omp_set_num_threads(8);  // кол во тредов

    matrix<double> a(n, n), a1(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв
    vector<double> bc(n);
    vector<double> b(n), b1(n); //вектора нагрузки для ур Лапласа и Пуассона соотв

    bc = Construct_BC(q);
    b = Construct_load(q, h, bc);
    a = Construct_matrix_Laplace(q);

    // настраиваем вывод
    std::cout.precision(3);

    //реализация с BicGstab c OMP
    #pragma omp parallel
    {
        u = BiCGSTAB(a, b); // решаем уравнение Лапласа
    }
    Print_matrix(a);
    Print_vectors(q, u);
    Print_vectors(q, b);

    return 0;
}
