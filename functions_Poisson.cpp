//
// Created by igor on 21.09.17.
//
#include "global.cpp"
#include "functions_Poisson.h"
#include <fstream>
#include <iomanip>


using namespace std::chrono;
using namespace Eigen;

//--------------------------------------------------------------------------------------
void Construct_matrix_Poisson(SparseMatrix<double> *a)
{

    int row = 0;
    a->reserve(VectorXi::Constant(n, 7));


    for (int k = 1; k <= qz; k++)
    {
        for (int j = 1; j <= qy; j++)
        {
            for (int i = 1; i <= qx; i++)
            {
                if (k > 1)       a->insert(row, row - qx * qy) = 1;
                if (j > 1)       a->insert(row, row - qx) = 1;
                if (i > 1)       a->insert(row, row - 1) = 1;
                                 a->insert(row, row) = -6;
                if (k < qz)       a->insert(row, row + qx * qy) = 1;
                if (j < qy)      a->insert(row, row + qx) = 1;
                if (i < qx)      a->insert(row, row + 1) = 1;
                row++;

            }
        }

    }


    a->makeCompressed();

}

//--------------------------------------------------------------------------------------


void Construct_guess_P(VectorXd *initGuess)
{
    initGuess->fill(0);
    for (int i = 0; i < n; i++)
        initGuess->coeffRef(i) = -1.0 / 6.0;
}

//--------------------------------------------------------------------------------------

void Construct_f(VectorXd *f, VectorXd * uL)
{
    VectorXd E(qz), v(qz), lamda(qz), G(qz); // упругие параметры
    // заполняем вектотора упругих параметров в разных слоях соотв функциями
    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);

    //заполняем вектор правой части ур я пуассона для каждого слоя
    for (int i = 0; i < qz; i++)
    {
        for (int j = 0; j < qx*qy; j++)
        {
            if ( i > 0  && i < (qz-1) )
                 f->coeffRef(i*qz+j) = ( -G[i] - lamda[i] )/G[i]*(uL->coeff(i*qx*qy+j+qx*qy) - uL->coeff(i*qx*qy+j-qx*qy))/2/h;

            else if ( i == 0 )
                        f->coeffRef(i*qz+j) = ( -G[i] - lamda[i] )/G[i]*(uL->coeff(i*qx*qy+j+qx*qy))/2/h;

            else if ( i == (qz-1) )
                                f->coeffRef(i*qz+j) = ( -G[i] - lamda[i] )/G[i]*(-uL->coeff(i*qx*qy+j-qx*qy))/2/h;

        }
    }
}

//--------------------------------------------------------------------------------------


// задаем граничные условия
void Construct_BC_Poisson(VectorXd *bc)
{
    bc->fill(0.0);

    VectorXd E(qz), v(qz), lamda(qz), G(qz) ; // упругие параметры
    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);

    for (int i = 0; i < qx * qy; i++)
    {
        if ( ((i%(qx)-qx/2.0)*(i%(qx)-qx/2.0)/A/A + (i/(qx)-qx/2.0)*(i/(qx)-qx/2.0)/B/B  ) < 0.9)
        {
            bc->coeffRef(i) =   4e7/G[0]/E[0]*(1-v[0])*B*sqrt(1- ((i%(qx)-qx/2.0)*(i%(qx)-qx/2.0)/A/A + (i/(qx)-qx/2.0)*(i/(qx)-qx/2.0)/B/B  ));

            //           bc->coeffRef(i) =  1.0*(1- (i%(qx)-qx/2)*(i%(qx)-qx/2)/A/A - (i/(qx)-qx/2)*(i/(qx)-qx/2)/B/B  ) ;
        }
    }

}

//--------------------------------------------------------------------------------------


void Construct_load_Poisson(VectorXd *bP, VectorXd *bcP, VectorXd *f)
{
    bP->fill(0);
    for (int i = 0; i < qx * qy; i++)
    {
        bP->coeffRef(i) = - bcP->coeff(i) + f->coeff(i) * h * h;
    }
}

//--------------------------------------------------------------------------------------


VectorXd Solve_Poisson(VectorXd * uL)
{
    SpMat AP(n, n);
    VectorXd bcP(n);
    VectorXd uP(n); //неизв векторы ур ия Пуассона
    VectorXd bP(n);
    VectorXd f(n);
    VectorXd initGuess(n); // начальное значение для солвера

    Construct_f(&f, uL);

    Construct_guess_P(&initGuess);
    Construct_matrix_Poisson(&AP);
    Construct_BC_Poisson(&bcP);
    Construct_load_Poisson(&bP, &bcP, &f);


//--------решаем уравнение Пуассона в каждом слое, где упругие параметры = const

    BiCGSTAB<SparseMatrix<double, RowMajor> > solverP;
// relative residual error: |Ax-b|/|b|
    solverP.setTolerance(error);
    high_resolution_clock::time_point tP1 = high_resolution_clock::now();


// раскладываем матрицу
    solverP.compute(AP);
// устанавливаем начальное приближение
    solverP.solveWithGuess(bP, initGuess);
//запускаем солвер
    //--------------
        uP = solverP.solve(bP);
   //--------------
    high_resolution_clock::time_point tP2 = high_resolution_clock::now();
//считаем время решения
    auto durationP = duration_cast<seconds>(tP2 - tP1).count();
    std::cout << std::endl <<"Poisson  duration = " << durationP << " || "<<"iterations = " << solverP.iterations()<< std::endl;

    return uP;

}

//--------------------------------------------------------------------------------------
//считаем производную dw/dz
void Construct_w_Derivative_z(  VectorXd * w_Derivative_z, VectorXd * w )
{
    w_Derivative_z->fill(0.0);
    for (int i = 0; i < qz; ++i)
    {
        for (int j = 0; j <qx*qy ; ++j)
        {
            if ( i > 0  && i < (qz-1) )
                w_Derivative_z->coeffRef(i*qz+j) = (w->coeff(i*(qz+1)+j) - w->coeff(i*(qz-1)+j))/2/h;

            else if ( i == 0 )
                w_Derivative_z->coeffRef(i*qz+j) = (-w->coeff(i*qz+j))/2/h;

            else if ( i == qz-1 )
                w_Derivative_z->coeffRef(i*qz+j) = (w->coeff(i*qz+j))/2/h;

        }
    }

}