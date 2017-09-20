
#include <iostream>
#include "functions.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include "../../My_linalg/my_linalg.h"
#include <chrono>



using namespace Eigen;
using namespace std::chrono;

int main(int argc, char **argv)
{


    omp_set_num_threads(omp_get_max_threads());
    Eigen::setNbThreads(omp_get_max_threads());

    SpMat AL(n, n);
    SpMat AP(n, n); // матрицы для решения уравнений Лапласа и Пуассона соотв

    SparseVector<double>  bc(n);
    VectorXd bc1(n); // граничные условия для задач Лапласа и Дирихле
    VectorXd b(n), b1(n); //векторы нагрузки для ур Лапласа и Пуассона соотв
    VectorXd f(n); // правая часть уравнения пуассона
    VectorXd initGuess(n); // начальное значение для солвера

//  заполняес параметры системы
    Construct_guess( &initGuess);
    Construct_matrix_Laplace(&AL);
    Construct_matrix_Poisson(&AP);
    Construct_BC_Laplace(&bc);
    Construct_load_Laplace( &b, &bc );
    Construct_load_Poisson( &b1, &bc1, &f );

//--------решаем уравнение Пуассона в каждом слое, где упругие параметры = const

    BiCGSTAB< SparseMatrix<double,RowMajor>, Eigen::IncompleteLUT< double > > solverP;
// устанавливаем треб точность
    solverP.setTolerance(error);
// считаем предобуславливатель
    solverP.preconditioner().setFillfactor(7);
    solverP.preconditioner().compute(AP);
// раскладываем матрицу
    solverP.compute(AP);
// устанавливаем начальное приближение
    solverP.solveWithGuess(b, initGuess);
//запускаем солвер
    high_resolution_clock::time_point tP1 = high_resolution_clock::now();
    u = solverP.solve(b);
    high_resolution_clock::time_point tP2 = high_resolution_clock::now();
//считаем время решения
    auto durationP = duration_cast<seconds>( tP2 - tP1 ).count();
    std::cout << std::endl << durationP << std::endl;

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
