//
// Created by gomonov on 11.09.17.
//
#include "functions_Laplace.h"
#include "global.cpp"
#include <fstream>
#include <iomanip>


using namespace std::chrono;
using namespace Eigen;

//--------------------------------------------------------------------------------------
//задаем  начальное раскрытие

void Construct_w0( VectorXd * w0 ) {
    w0->fill(0.0);

    for (int i = 0; i < qx * qy; i++)
    {
        if ( ((i%(qx)-qx/2)*(i%(qx)-qx/2)/A/A + (i/(qx)-qx/2)*(i/(qx)-qx/2)/B/B  ) < 0.9)
        {
            w0->coeffRef(i) =  1.0   ;

        }

    }
}

//--------------------------------------------------------------------------------------


// собираем матрицу СЛАУ для ур ия лапласа
void Construct_matrix_Laplace(SparseMatrix<double> *a)
{
    long row = 0;
    a->reserve(VectorXi::Constant(n, 7));


    for (int k = 1; k <= qz; k++)
    {
        for (int j = 1; j <= qy; j++)
        {
            for (int i = 1; i <= qx; i++)
            {
                if (k > 1)                      a->insert(row, row - qx * qy)   = 1;
                if (j > 1)                      a->insert(row, row - qx)        = 1;
                if (i > 1)                      a->insert(row, row - 1)         = 1;

                                                a->insert(row, row)             = -6;

                if (k < qz && row >= qx * qy)   a->insert(row, row + qx * qy)   = 1;
                if (j < qy)                     a->insert(row, row + qx)        = 1;
                if (i < qx)                     a->insert(row, row + 1)         = 1;
                //implementing Neumann B.C.
                if (k < qz && row < qx * qy)    a->insert(row, row + qx * qy)   = 2;
                row++;

            }
        }

    }

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
void Construct_BC_Laplace(VectorXd *bc, VectorXd * w0)
{
    std::vector<double> E(qz,0), v(qz,0), lamda(qz,0), G(qz,0); // упругие параметры
    // заполняем вектотора упругих параметров в разных слоях соотв функциями
    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);

    bc->fill(0.0);

    for (int i = 1; i < qx * qy -  qx; i++)
    {
        bc->coeffRef(i) = 2*G[0] / (lamda[0]+2*G[0]) * ( w0->coeff(i-1) - 4*w0->coeff(i) +  w0->coeff(i+1) +   w0->coeff(i-qx) + w0->coeff(i+qx)  )  / h / h  ;

    }


}

//--------------------------------------------------------------------------------------


void Construct_load_Laplace(VectorXd *b, VectorXd *bc)
{
    b->fill(0);
    for (int i = 0; i < qx * qy; i++)
    {
        b->coeffRef(i) = 2 * h * bc->coeff(i);
    }


}

//--------------------------------------------------------------------------------------


VectorXd Solve_Laplace()
{
    SpMat AL(n, n);
    VectorXd bc(n);
    VectorXd u(n); //неизв векторы ур ия Лапласа
    VectorXd b(n);
    VectorXd initGuess(n); // начальное значение для солвера
    VectorXd w0(qx*qy); // начальное поле перемещения по Oz

    Construct_w0(&w0);
    Construct_guess_L(&initGuess);
    Construct_matrix_Laplace(&AL);
    Construct_BC_Laplace(&bc, &w0);
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

