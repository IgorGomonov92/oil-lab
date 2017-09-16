
#include <iostream>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
//#include "../../My_linalg/my_linalg.h"

using namespace Eigen;

int main(int argc, char **argv)
{

    Eigen::setNbThreads(omp_get_max_threads());

    VectorXd u, u1; //неизв векторы ур ий Лапласа и Пуассона
    double h=1.0; // шаг сетки
    SpMat AL(n, n);
    SpMat AP(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв
    Eigen::IncompleteLUT< double > ILUMAT;
    SparseVector<double>  bc(n) , bc1(n); // граничные условия для задач Лапласа и Дирихле
    SparseVector<double> b(n), b1(n); //векторы нагрузки для ур Лапласа и Пуассона соотв
    SparseVector<double> f(n); // правая часть уравнения пуассона


    Construct_matrix_Laplace(&AL);
    Construct_BC_Laplace(&bc);
    Construct_load_Laplace(h, &b, &bc);

    ILUMAT.setFillfactor(7);
    ILUMAT.compute(AL);

    BiCGSTAB< SparseMatrix<double,RowMajor>, Eigen::IncompleteLUT< double > > solver;

    solver.preconditioner();

    solver.setTolerance(2.e-5);
    u = solver.compute( AL ).solve(b);

    std::cout<< AL;
    std::cout<< u;
//    A1 = Construct_matrix_Poisson();

/*

    bc1 = Construct_BC_Poisson();


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
