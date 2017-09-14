#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"

//размеры сетки - глоб константы


using namespace boost::numeric::ublas;

int main(int argc, char **argv)
{
    extern int n;


    vector<double> u(n), u1(n); //неизв векторы ур ий Лапласа и Пуассона
    double h=1.0; // шаг сетки

    matrix<double> a(n, n), a1(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв
    vector<double> bc(n), bc1(n); // граничные условия для задач Лапласа и Дирихле
    vector<double> b(n), b1(n); //векторы нагрузки для ур Лапласа и Пуассона соотв
    vector<double> f(n, .5); // правая часть уравнения пуассона


    bc = Construct_BC_Laplace();
    bc1 = Construct_BC_Poisson();
    b = Construct_load_Laplace(h, bc);
    a = Construct_matrix_Laplace();
    Print_matrix(a);

    a1 = Construct_matrix_Poisson();

    // настраиваем вывод
    std::cout.precision(3);

    //реализация BicGstab

            u = BiCGSTAB(a, b); // решаем уравнение Лапласа
            b1 = Construct_load_Poisson(h, bc1, f);
            u1 = BiCGSTAB(a1, b1); // решаем уравнение пуассона



    Print_matrix(a1);
    Print_vectors(u1);
    Print_vectors(b);

    return 0;
}
