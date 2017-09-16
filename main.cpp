
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include "../../My_linalg/my_linalg.h"

using namespace Eigen;

int main(int argc, char **argv)
{


    RowVectorXd u, u1; //неизв векторы ур ий Лапласа и Пуассона
    double h=1.0; // шаг сетки
    SpMat A1(n, n);
    SpMat  A(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв
    SparseVector<double>  bc(n) , bc1(n); // граничные условия для задач Лапласа и Дирихле
    SparseVector<double> b(n), b1(n); //векторы нагрузки для ур Лапласа и Пуассона соотв
    SparseVector<double> f(n); // правая часть уравнения пуассона



      Construct_matrix_Laplace(&A);
    std::cout<< A;
//    A1 = Construct_matrix_Poisson();

/*
    bc = Construct_BC_Laplace();
    bc1 = Construct_BC_Poisson();
    b = Construct_load_Laplace(h, bc);

    // настраиваем вывод
    std::cout.precision(3);

    //реализация BicGstab

            u = BiCGSTAB(A, b); // решаем уравнение Лапласа
            b1 = Construct_load_Poisson(h, bc1, f);
            u1 = BiCGSTAB(A1, b1); // решаем уравнение пуассона



    Print_matrix(A1);
    Print_vectors(u1);
    Print_vectors(b);
*/
        return 0;
}
