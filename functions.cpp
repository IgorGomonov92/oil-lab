//
// Created by gomonov on 11.09.17.
//
#include "functions.h"
#include <omp.h>


#include "/home/igor/My_linalg/my_linalg.cpp"
#include "/home/igor/Eigen/Eigen/SparseCore"


using namespace Eigen;





// собираем матрицу СЛАУ для ур ия лапласа
void Construct_matrix_Laplace(SparseMatrix<double > * a )
{

    int row=0;
   // SpMat a(n,n);// матрица СЛАУ
    a->reserve(VectorXi::Constant(n,7));


    for(int k=1; k<=qz; k++)
    {
        for(int j=1; j<=qy; j++)
        {
            for(int i=1; i<=qx; i++)
            {
                if (k>1)        a->insert(row,row-qx*qy) = 1;
                if (j>1)        a->insert(row,row-qx) = 1;
                if (i>1)        a->insert(row,row-1) = 1;
                a->insert(row,row) = -6;
                if (k<qz && row>=qx*qy)       a->insert(row,row+qx*qy) = 1;
                if (j<qy)       a->insert(row,row+qx) = 1;
                if (i<qx)       a->insert(row,row+1) = 1;
                //implementing Neumann B.C.
                if (k<qz && row<qx*qy)    a->insert(row,row+qx*qy) = 2;
                row++;

            }
        }

    }


    a->makeCompressed();
    //return a;
}

void Construct_matrix_Poisson(SparseMatrix<double > * a )
{

    int row=0;
    // SpMat a(n,n);// матрица СЛАУ
    a->reserve(VectorXi::Constant(n,7));


    for(int k=1; k<=qz; k++)
    {
        for(int j=1; j<=qy; j++)
        {
            for(int i=1; i<=qx; i++)
            {
                if (k>1)        a->insert(row,row-qx*qy) = 1;
                if (j>1)        a->insert(row,row-qx) = 1;
                if (i>1)        a->insert(row,row-1) = 1;
                a->insert(row,row) = -6;
                if (k<qz && row>=qx*qy)       a->insert(row,row+qx*qy) = 1;
                if (j<qy)       a->insert(row,row+qx) = 1;
                if (i<qx)       a->insert(row,row+1) = 1;
                row++;

            }
        }

    }


    a->makeCompressed();
    //return a;
}

void Construct_guess(VectorXd * initGuess)
{
    for(int i=0;i<n;i++)
        initGuess->coeffRef(i) = -1.0/6.0;
}

/*
// задаем граничный условия
*/
void Construct_BC_Laplace(SparseVector<double> * bc)
{

    for(int i=0; i<qx*qy; i++)
    {
        bc->insert(i) = .2;
    }

}


// задаем граничные условия

void Construct_BC_Poisson(SparseVector<double> * bc)
{
    for(int i=0; i<qz; i++)
    {
        bc->insert(i,0) = .1;
    }

}

void Construct_load_Laplace(double h, VectorXd * b, SparseVector<double> * bc)
{
    b->fill(0);
    for(int i=0 ; i<qx*qy; i++)
    {
        b->coeffRef(i) = 2*h*bc->coeff(i);
    }


}


void Construct_load_Poisson(double h, VectorXd * b, SparseVector<double> * bc, VectorXd * f)
{
    b->fill(0);
    for(unsigned long i=0 ; i<qx*qy; i++)
    {
        b->coeffRef(i) = bc->coeff(i) + f->coeffRef(i)*2.0*h;

    }


}



