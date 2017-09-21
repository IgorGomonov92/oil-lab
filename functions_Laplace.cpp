//
// Created by gomonov on 11.09.17.
//
#include "functions_Laplace.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include <chrono>
#include "global.cpp"

using namespace std::chrono;
using namespace Eigen;

//--------------------------------------------------------------------------------------


// собираем матрицу СЛАУ для ур ия лапласа
void Construct_matrix_Laplace(SparseMatrix<double> *a)
{

    int row = 0;
    a->reserve(VectorXi::Constant(n, 7));


    for (int k = 1; k <= qz; k++) {
        for (int j = 1; j <= qy; j++) {
            for (int i = 1; i <= qx; i++) {
                if (k > 1) a->insert(row, row - qx * qy) = 1;
                if (j > 1) a->insert(row, row - qx) = 1;
                if (i > 1) a->insert(row, row - 1) = 1;
                a->insert(row, row) = -6;
                if (k < qz && row >= qx * qy) a->insert(row, row + qx * qy) = 1;
                if (j < qy) a->insert(row, row + qx) = 1;
                if (i < qx) a->insert(row, row + 1) = 1;
                //implementing Neumann B.C.
                if (k < qz && row < qx * qy) a->insert(row, row + qx * qy) = 2;
                row++;

            }
        }

    }

std::cout<< *a;
    a->makeCompressed();
}

//--------------------------------------------------------------------------------------


void Construct_guess_L(VectorXd *initGuess)
{
    initGuess->fill(0);
    for (int i = 0; i < n; i++)
        initGuess->coeffRef(i) = -1.0 / 6.0;
}

//--------------------------------------------------------------------------------------


/*
// задаем граничный условия
*/
void Construct_BC_Laplace(SparseVector<double> *bc)
{
    for (int i = 0; i < qx * qy; i++) {
        bc->insert(i) = .2;
    }

}

//--------------------------------------------------------------------------------------


void Construct_load_Laplace(VectorXd *b, SparseVector<double> *bc)
{
    b->fill(0);
    for (int i = 0; i < qx * qy; i++) {
        b->coeffRef(i) = h * h * bc->coeff(i);
    }


}

//--------------------------------------------------------------------------------------


VectorXd Solve_Laplace()
{
    SpMat AL(n, n);
    SparseVector<double> bc(n);
    VectorXd u(n); //неизв векторы ур ия Лапласа
    VectorXd b(n);
    VectorXd initGuess(n); // начальное значение для солвера
    Construct_guess_L(&initGuess);
    Construct_matrix_Laplace(&AL);
    Construct_BC_Laplace(&bc);
    Construct_load_Laplace(&b, &bc);

    u.fill(0);
    // Решаем уравнение Лапласа
    BiCGSTAB<SparseMatrix<double, RowMajor> > solverL;
    high_resolution_clock::time_point tL1 = high_resolution_clock::now();

// relative residual error: |Ax-b|/|b|
    solverL.setTolerance(error);

    solverL.compute(AL);
// устанавливаем начальное приближение
    solverL.solveWithGuess(b, initGuess);
//запускаем солвер
    //--------------
    u = solverL.solve(b);
    //--------------
    high_resolution_clock::time_point tL2 = high_resolution_clock::now();
//считаем время решения
    auto durationL = duration_cast<seconds>(tL2 - tL1).count();
    std::cout << std::endl <<"Laplace  duration = " << durationL << " || "<<"iterations = " << solverL.iterations()<< std::endl;
// закончили обсчет ур я Лапласа

    return u;

}

//--------------------------------------------------------------------------------------

