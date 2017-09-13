#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
using namespace boost::numeric::ublas;


int main(int argc, char **argv)
{
    int qx, qy, qz; // кол во узлов по x y z
    int n;
    qx = 3;
    qy = 4;
    qz = 5;

    n = qx * qy * qz;

    vector<double> u(n), u1(n); //неизв вектора ур ий Лапласа и Пуассона
    double h=1.0; // шаг сетки

    matrix<double> a(n, n), a1(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв
    vector<double> bc(n), bc1(n); // граничные условия для задач Лапласа и Дирихле
    vector<double> b(n), b1(n); //вектора нагрузки для ур Лапласа и Пуассона соотв
    vector<double> f(n, .5); // правая часть уравнения пуассона


    bc = Construct_BC_Laplace(n);
    bc1 = Construct_BC_Poisson(n);
    b = Construct_load_Laplace(qx, qy, qz, h, bc);
    a = Construct_matrix_Laplace(qx, qy, qz);

    a1 = Construct_matrix_Poisson(qx, qy, qz);

    // настраиваем вывод
    std::cout.precision(3);

    //реализация с BicGstab c OMP

            u = BiCGSTAB(a, b); // решаем уравнение Лапласа
            b1 = Construct_load_Poisson(qx, qy, qx, h, bc1, f);
            u1 = BiCGSTAB(a1, b1); // решаем уравнение пуассона



    Print_matrix(a1);
    Print_vectors(qx, qy, qz, u1);
    Print_vectors(qx, qy, qz, b);
    std::cout << bc1<<std::endl;
    std::cout << b1<<std::endl;
    std::cout << f<<std::endl;

    return 0;
}
