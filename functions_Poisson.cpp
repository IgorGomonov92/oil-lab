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
    a->reserve(VectorXi::Constant(nP, 7));

//решаем задачу в одном слое
    for (int k = 1; k <= 2; k++)
    {
        for (int j = 1; j <= qy; j++)
        {
            for (int i = 1; i <= qx; i++)
            {
                if (k > 1) a->insert(row, row - qx * qy) = 1;
                if (j > 1) a->insert(row, row - qx) = 1;
                if (i > 1) a->insert(row, row - 1) = 1;
                a->insert(row, row) = -6;
                if (k < 2) a->insert(row, row + qx * qy) = 1;
                if (j < qy) a->insert(row, row + qx) = 1;
                if (i < qx) a->insert(row, row + 1) = 1;
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
    for (int i = 0; i < nP; i++)
        initGuess->coeffRef(i) = -1.0 / 6.0;
}

//--------------------------------------------------------------------------------------

void Construct_f(std::vector<VectorXd> *f, VectorXd * uL)
{
    std::vector<double> E(qz), v(qz), lamda(qz), G(qz); // упругие параметры
    // заполняем вектотора упругих параметров в разных слоях соотв функциями
    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);

    f->at(0).resize(nP);

    //заполняем вектор правой части ур я пуассона для каждого слоя
    for (int i = 0; i < qz; i++)
    {
        f->at(i).resize(nP);

        for (int j = 0; j < qx*qy; j++)
        {
           // std::cout<<f->at(i).coeffRef(j)<< std::endl;
            if ( i > 0  && i < (qz-1) )
                 f->at(i).coeffRef(j) = ( -G[i] - lamda[i] )/G[i]*(uL->coeff(i*qx*qy+j+qx*qy) - uL->coeff(i*qx*qy+j-qx*qy))/2/h;

                    else if ( i == 0 )
                        f->at(i).coeffRef(j) = ( -G[i] - lamda[i] )/G[i]*(uL->coeff(i*qx*qy+j+qx*qy))/2/h;

                            else if ( i == (qz-1) )
                                f->at(i).coeffRef(j) = ( -G[i] - lamda[i] )/G[i]*(-uL->coeff(i*qx*qy+j-qx*qy))/2/h;

        }
    }
}

//--------------------------------------------------------------------------------------


// задаем граничные условия
void Construct_BC_Poisson(VectorXd *bc)
{
    bc->fill(0.0);
    for (int i = 0; i < qx * qy; i++)
    {
        if ( ((i%(qx)-qx/2)*(i%(qx)-qx/2)/100.0 + (i/(qx)-qx/2)*(i/(qx)-qx/2)/400.0  ) <= 1)
        {
            bc->coeffRef(i) = -10.1;

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


std::vector<VectorXd> Solve_Poissons(VectorXd * uL)
{
    SpMat AP(nP, nP);
    VectorXd bcP(nP);
    std::vector<VectorXd> uP(qz + 1); //неизв векторы ур ия Пуассона
    VectorXd bP(nP);
    std::vector<VectorXd> f(qz + 1);
    VectorXd *bc_prom_P = nullptr; //промежуточный вектор
    VectorXd initGuess(nP); // начальное значение для солвера

    Construct_f(&f, uL);

    //bc_prom_P->resize(nP);

    Construct_guess_P(&initGuess);
    Construct_matrix_Poisson(&AP);
    Construct_BC_Poisson(&bcP);
    Construct_load_Poisson(&bP, &bcP, &f[0]);


    std::ofstream outputB ("B.txt");
    //for (int k = 0; k < qz; ++k)

    for (int l = 0; l < qy; ++l)

        for (int m = 0; m < qx; ++m)
        {
            outputB<< std::scientific << std::setprecision(5) <<bcP.coeff(l*qy+m)<< std::setw(10) <<" "<< m+1 << " " << l+1 << " " << 1 << std::endl;
        }
    outputB.close();


    for (int i = 0; i < qz; ++i)
    {
        uP[i].resize(nP);
        uP[i].fill(0);
    }

//--------решаем уравнение Пуассона в каждом слое, где упругие параметры = const

    BiCGSTAB<SparseMatrix<double, RowMajor>, Eigen::IncompleteLUT<double> > solverP;
// relative residual error: |Ax-b|/|b|
    solverP.setTolerance(error);
// считаем предобуславливатель
    high_resolution_clock::time_point tP1 = high_resolution_clock::now();

    solverP.preconditioner().setFillfactor(7);

    solverP.preconditioner().compute(AP);
// раскладываем матрицу
    solverP.compute(AP);
// устанавливаем начальное приближение
    solverP.solveWithGuess(bP, initGuess);
//запускаем солвер
    //--------------
    for (int i = 0; i < qz; i++)
    {
        uP[i] = solverP.solve(bP);
// граничное условие на след шаге по z равно первому слою решения на предыдущем шаге
        bc_prom_P = &uP[i];

        // bc_prom_P->fill(0.0);
     /*   for (int j = 0; j < qx*qy; ++j)
        {
            bc_prom_P->coeffRef(j) = uP[i].coeff(j);
        }

        */
// считаем новую правую часть
        Construct_load_Poisson(&bP, bc_prom_P, &f[i]);
    }





    //--------------
    high_resolution_clock::time_point tP2 = high_resolution_clock::now();
//считаем время решения
    auto durationP = duration_cast<seconds>(tP2 - tP1).count();
    std::cout << std::endl <<"Poisson  duration = " << durationP << " || "<<"iterations = " << solverP.iterations()<< std::endl;

    return uP;

}

//--------------------------------------------------------------------------------------

void Construct_w_Derivative_z(  std::vector<VectorXd> * w_Derivative_z, std::vector<VectorXd> * w )
{
    for (int i = 0; i < qz; ++i)
    {
        w_Derivative_z->at(i).resize(qx*qy);
        w_Derivative_z->at(i).fill(.0);

        for (int j = 0; j <qx*qy ; ++j)
        {
            if ( i > 0  && i < (qz-1) )
                w_Derivative_z->at(i).coeffRef(j) = (w->at(i+1).coeff(j) - w->at(i-1).coeff(j))/2/h;

                     else if ( i == 0 )
                          w_Derivative_z->at(0).coeffRef(j) = (-w->at(1).coeff(j))/2/h;

                            else if ( i == qz-1 )
                                w_Derivative_z->at(qz-1).coeffRef(j) = (w->at(qz-1).coeff(j))/2/h;

        }
    }

}