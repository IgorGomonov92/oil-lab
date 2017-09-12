#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <omp.h>
#include "functions.h"
using namespace boost::numeric::ublas;


int main(int argc, char **argv)
{
    int q; // величина дискретизации
    int n;
    q = 3;
    n = q * q * q;

    vector<double> u(n), u1(n); //неизв вектора ур ий Лапласа и Пуассона
    double h=1.0; // шаг сетки



    omp_set_num_threads(8);  // кол во тредов

    matrix<double> a(n, n), a1(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв
    vector<double> bc(n), bc1(n); // граничные условия для задач Лапласа и Дирихле
    vector<double> b(n), b1(n); //вектора нагрузки для ур Лапласа и Пуассона соотв
    vector<double> f(n, .5); // правая часть уравнения пуассона


    bc = Construct_BC_Laplace(q);
    bc1 = Construct_BC_Poisson(q);
    b = Construct_load_Laplace(q, h, bc);
    a = Construct_matrix_Laplace(q);

    a1 = Construct_matrix_Poisson(q);

    // настраиваем вывод
    std::cout.precision(3);

    //реализация с BicGstab c OMP
    #pragma omp parallel
    {
        u = BiCGSTAB(a, b); // решаем уравнение Лапласа
        b1 = Construct_load_Poisson(q, h, bc1, f);
        u1 = BiCGSTAB(a1, b1); // решаем уравнение пуассона
    }
    Print_matrix(a1);
    Print_vectors(q,u1);
    Print_vectors(q,b);
    std::cout << bc1<<std::endl;
    std::cout << b1<<std::endl;
    std::cout << f<<std::endl;

    return 0;
}
