//
// Created by gomonov on 11.09.17.
//
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"

#include "/home/igor/Eigen/Eigen/SparseCore"


using namespace Eigen;



// реализация метода BiCGSTAB

SpMat BiCGSTAB(SpMat a, SpMat b) {
        SpMat r, r0, u0, rt, p, p0, v, v0, prom1, prom2, s, t, u;
        double alpha0 = 1, alpha, beta = 0, ro0 = 1, ro, w, w0 = 1;
/*
        p.resize(a.size1(), true);
        u0.resize(a.size1(), true);
        v0.resize(a.size1(), true);
        p0.resize(a.size1(), true);
        r0.resize(a.size1(), true);
        u.resize(a.size1(), true);

        r0 = b - prod(a, u);
        rt = r0;
    int i=0;
        do {
            i++;
                ro = inner_prod(rt, r0);
                beta = ro / ro0 * alpha0 / w0;
                prom1 = w0 * v0;
                prom2 = p0 - prom1;
                p = r0 + beta * prom2;
                v = prod(a, p);
                alpha = ro0 / inner_prod(rt, v);
                s = r0 - alpha * v;
                t = prod(a, s);
                w = inner_prod(t, s) / inner_prod(t, t);
                u = u0 + w * s + alpha * p;
                r = s - w * t;


                ro0 = ro;
                alpha0 = alpha;
                w0 = w;
                r0 = r;
                p0 = p;
                v0 = v;
                u0 = u;

        } while (norm_2(r) / norm_2(b) > 1.0e-5);
        std::cout <<std::endl << i << std::endl;
        return u;
*/
}



// собираем матрицу СЛАУ для ур ия лапласа
SpMat Construct_matrix_Laplace()
{
    std::vector<T> tripletList;
    int row=0;
    SpMat a(n, n);// матрица СЛАУ

    for(int k=1; k<=qz; k++)
    {
        for(int j=1; j<=qy; j++)
        {
            for(int i=1; i<=qx; i++)
            {
                if (k>1)        tripletList.push_back(T(row,row-qx*qy, 1));
                if (j>1)        tripletList.push_back(T(row,row-qx, 1));
                if (i>1)        tripletList.push_back(T(row,row-1, 1));
                tripletList.push_back(T(row,row, -6));
                if (k<qz)       tripletList.push_back(T(row,row+qx*qy, 1));
                if (j<qy)       tripletList.push_back(T(row,row+qx, 1));
                if (i<qx)       tripletList.push_back(T(row,row+1, 1));
                //implementing Neumann B.C.
                if (k<qz && row<qx*qy)        tripletList.push_back(T(row,row+qx*qy, 1));
                row++;
            }
        }

    }
    a.setFromTriplets(tripletList.begin(), tripletList.end());
    a.makeCompressed();
    return a;
}

SpMat Construct_matrix_Poisson()
{
    std::vector<T> tripletList;
    int row=0;
    SpMat a(n, n);// матрица СЛАУ

    for(int k=1; k<=qz; k++)
    {
        for(int j=1; j<=qy; j++)
        {
            for(int i=1; i<=qx; i++)
            {
                if (k>1)        tripletList.push_back(T(row,row-qx*qy, 1));
                if (j>1)        tripletList.push_back(T(row,row-qx, 1));
                if (i>1)        tripletList.push_back(T(row,row-1, 1));
                tripletList.push_back(T(row,row, -6));
                if (k<qz)       tripletList.push_back(T(row,row+qx*qy, 1));
                if (j<qy)       tripletList.push_back(T(row,row+qx, 1));
                if (i<qx)       tripletList.push_back(T(row,row+1, 1));
                row++;
            }
        }

    }

    a.setFromTriplets(tripletList.begin(), tripletList.end());
    a.makeCompressed();
    return a;
}

/*
// вывод матрицы
void Print_matrix(SpMat a)
{
    std::vector<T> tripletList;
    std::cout << std::endl;
    for(int i=0; i<a.size(); i++)
    {
        for(int j=0; j<a.size(); j++)
        {

            //  std::cout << tripletList.at(i,j) << ' ';

        }
        std::cout << std::endl;
    }

}


// создаем нагрузку

SpMat Construct_load_Laplace(double h, SpMat bc)
{
    std::vector<T> tripletList;
    SpMat b(n, .0);
    for(int i=0 ; i<n; i++)
    {
        if (i < qx*qy)    b(i) = 2*h*bc(i);
    }

    return b;
}

// задаем граничный условия

SpMat Construct_BC_Laplace()
{
    SpMat bc(n, .0);
    for(int i=0; i<n; i++)
    {
        bc(i) = .2;
    }

    return bc;
}


// задаем граничные условия

SpMat Construct_BC_Poisson()
{
    SpMat bc1(n, .0);
    for(int i=0; i<n; i++)
    {
        bc1(i) = .1;
    }

    return bc1;
}

// print unknowns
void Print_vectors(SpMat u)
{
    std::cout << std::endl;
    for(int i=0; i<u.size(); i++)
    {
        std::cout << u[i]<< " ";

        if((i + 1) % (qx) == 0)  std::cout << std::endl;
        if((i + 1) % (qx*qy) == 0) std::cout << std::endl;

    }
    std::cout << std::endl;
}



SpMat Construct_load_Poisson(double h, SpMat bc, SpMat f)
{
    SpMat b(n, .0);
    for(unsigned long i=0 ; i<n; i++)
    {
        if (i < qx*qy)    b(i) = bc(i) + f(i)*2.0*h;

    }

    return b;
}

*/

