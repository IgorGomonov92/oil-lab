//
// Created by gomonov on 11.09.17.
//
#include "functions_Laplace.h"
#include "global.cpp"



using namespace std::chrono;
using namespace Eigen;

//--------------------------------------------------------------------------------------

void Construct_w( VectorXd * w )
{
    for (int i = 0; i < qx*qy; ++i)
    {
        w->coeffRef(i) = 12;
    }
}

//--------------------------------------------------------------------------------------


// собираем матрицу СЛАУ для ур ия лапласа
void Construct_matrix_Laplace(SparseMatrix<double> *a)
{
    int row = 0;
    a->reserve(VectorXi::Constant(n, 7));


    for (int k = 1; k <= qz; k++)
    {
        for (int j = 1; j <= qy; j++)
        {
            for (int i = 1; i <= qx; i++)
            {
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
void Construct_BC_Laplace(SparseVector<double> *bc, VectorXd * w)
{
    std::vector<double> E(qz), v(qz), lamda(qz), G(qz); // упругие параметры
    // заполняем вектотора упругих параметров в разных слоях соотв функциями
    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);

    //задаем граничные условия Неймана
    bc->insert(0) = 2*G[0] / (lamda[0]+2*G[0]) * ( 4*w->coeff(0) + w->coeff(1) + w->coeff(qx) ) / h / h;
    bc->insert(qx * qy) = 2*G[qx*qy] / (lamda[qx*qy]+2*G[qx*qy]) * (w->coeff(qx*qy-1) - 4*w->coeff(qx*qy) + w->coeff(qx*qy-qx)  ) / h / h;


    for (int i = 1; i < qx * qy - 1; i++)
    {
        bc->insert(i) = 2*G[i] / (lamda[i]+2*G[i]) * ( w->coeff(i-1) - 4*w->coeff(i) + w->coeff(i+1) + w->coeff(i-qx) + w->coeff(i+qx) ) / h / h;


    }

}

//--------------------------------------------------------------------------------------


void Construct_load_Laplace(VectorXd *b, SparseVector<double> *bc)
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
    SparseVector<double> bc(n);
    VectorXd u(n); //неизв векторы ур ия Лапласа
    VectorXd b(n);
    VectorXd initGuess(n); // начальное значение для солвера
    VectorXd w(qx*qy); // начальное поле перемещения по Oz

    Construct_w(&w);
    Construct_guess_L(&initGuess);
    Construct_matrix_Laplace(&AL);
    Construct_BC_Laplace(&bc, &w);
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

