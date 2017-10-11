//
// Created by igor on 21.09.17.
//
#include "global.cpp"
#include "functions_Poisson.h"
#include "functions_Laplace.h"
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
                if (k < qz)      a->insert(row, row + qx * qy) = 1;
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

    //заполняем вектор правой части ур я пуассона
    for (int i = 0; i < qz; i++)
    {
        for (int j = 0; j < qx*qy; j++)
        {
            if ( i > 0  && i < (qz-1) )
                f->coeffRef(i*qy*qx+j) = ( -G[i] - lamda[i] )/G[i]*(uL->coeff(i*qx*qy+j+qx*qy) - uL->coeff(i*qx*qy+j-qx*qy))/2.0/h;

            else if ( i == 0 )
                f->coeffRef(i*qy*qx+j) = (uL->coeff(i*qx*qy+j+2*qx*qy) - 4*uL->coeff(i*qx*qy+j+qx*qy) + 3*(uL->coeff(i*qx*qy+j)))/2.0/h;

            else if ( i == (qz-1) )
                f->coeffRef(i*qy*qx+j) = ( -G[i] - lamda[i] )/G[i]*(-uL->coeff(i*qx*qy+j-qx*qy))/2.0/h;

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


    for (int i = 0; i < qx ; i++)
    {
        for (int j = 0; j < qy ; ++j)
        {
            if ( (((((double)i-(double)qx/2.0)*((double)i-(double)qx/2.0)/A/A + ((double)j-qx/2.0)*((double)j-(double)qx/2.0)/B/B  )  ) < 1.0) )
            {
                bc->coeffRef(i*qx+j) =  4.0e6/G[0]/3.14*2*(1.0-v[0])*B*sqrt(1.0-  (((double)i-(double)qx/2.0)*((double)i-(double)qx/2.0)/A/A + ((double)j-qx/2.0)*((double)j-(double)qx/2.0)/B/B  )  );

            }

        }

    }

}

//--------------------------------------------------------------------------------------


void Construct_load_Poisson(VectorXd *bP, VectorXd *bcP, VectorXd *f)
{
    bP->fill(0);
    for (int i = 0; i < qx * qy; i++)
    {
        bP->coeffRef(i) = f->coeff(i) * h * h - bcP->coeff(i) ;
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

VectorXd Solve_Poisson_XX(VectorXd * uL)
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
   // Construct_BC_Poisson_XX(&bcP);
    //Construct_load_Poisson_XX(&bP, &bcP, &f);


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
    VectorXd w0(qx*qy); // начальное поле перемещения по Oz
    Construct_w0(&w0);

    w_Derivative_z->fill(0.0);
    for (int i = 0; i < qz; ++i)
    {
        for (int j = 0; j <qx*qy ; ++j)
        {
            if ( i > 0  && i < (qz-1) )
                w_Derivative_z->coeffRef(i*qy*qx+j) = (w->coeff(i*(qx*qy)+qx*qy+j) - w->coeff(i*(qx*qy)-qx*qy+j))/2.0/h;


            else if ( i == 0 )
                w_Derivative_z->coeffRef(i*qy*qx+j) = ( w->coeff(j+qx*qy) - w0.coeff(j) )/2.0/h;

            else if ( i == qz-1 )
                w_Derivative_z->coeffRef(i*qy*qx+j) = (w->coeff(i*qx*qy+j-qx*qy))/2.0/h;

        }
    }

}

//--------------------------------------------------------------------------------------
//считаем производную dw/dx
void Construct_w_Derivative_x(  VectorXd * w_Derivative_x, VectorXd * w )
{
    VectorXd w0(qx*qy); // начальное поле перемещения по Ox
    Construct_w0(&w0);

    w_Derivative_x->fill(0.0);
    for (int i = 0; i < qz; ++i)
    {
        for (int j = 0; j <qz*qy ; ++j)
        {
            if ( i > 0  && i < (qx-1) )
                w_Derivative_x->coeffRef(i*qy*qx+j) = (w->coeff(i*(qz*qy)+qz*qy+j) - w->coeff(i*(qz*qy)-qz*qy+j))/2.0/h;


            else if ( i == 0 )
                w_Derivative_x->coeffRef(i*qy*qz+j) = ( w->coeff(j+qz*qy) - w0.coeff(j) )/2.0/h;

            else if ( i == qx-1 )
                w_Derivative_x->coeffRef(i*qy*qz+j) = (w->coeff(i*qz*qy+j-qz*qy))/2.0/h;

        }
    }

}

