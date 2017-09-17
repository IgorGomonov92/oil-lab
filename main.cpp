
#include <iostream>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
//#include "../../My_linalg/my_linalg.h"
#include <chrono>

using namespace Eigen;
using namespace std::chrono;

int main(int argc, char **argv)
{


   // omp_set_num_threads(2);
    //Eigen::setNbThreads(2);

    VectorXd u, u1; //неизв векторы ур ий Лапласа и Пуассона
    double h=1.0; // шаг сетки
    SpMat AL(n, n);
    SpMat AP(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв

    VectorXd  bc(n) , bc1(n); // граничные условия для задач Лапласа и Дирихле
    VectorXd b(n), b1(n); //векторы нагрузки для ур Лапласа и Пуассона соотв
    VectorXd f(n); // правая часть уравнения пуассона

    b.fill(1);

    Construct_matrix_Laplace(&AL);
    //Construct_BC_Laplace(&bc);
    //Construct_load_Laplace(h, &b, &bc);


    BiCGSTAB< SparseMatrix<double,RowMajor>, Eigen::IncompleteLUT< double > > solver;

    solver.preconditioner().setFillfactor(7);
    solver.preconditioner().compute(AL);

    solver.setMaxIterations(100000);
    solver.setTolerance(2.e-45);

    solver.compute( AL );

    high_resolution_clock::time_point t1 = high_resolution_clock::now();


   // #pragma omp parallel
    u = solver.solve(b);


    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();

    std::cout << duration;

    std::cout<<' '<< solver.iterations() <<' ' << Eigen::nbThreads();

    //std::cout<< AL;
   // std::cout<< u;
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
