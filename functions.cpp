//
// Created by gomonov on 11.09.17.
//
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"

inline const int qx=3, qy=4, qz=5;
inline  const int n = qx * qy * qz;

using namespace boost::numeric::ublas;



// реализация метода BiCGSTAB

vector<double> BiCGSTAB(matrix<double> a, vector<double> b) {
        vector<double> r, r0, u0, rt, p, p0, v, v0, prom1, prom2, s, t, u;
        double alpha0 = 1, alpha, beta = 0, ro0 = 1, ro, w, w0 = 1;
        u.clear();
        u0.clear();
        v0.clear();
        p0.clear();
        r0.clear();
        p.clear();
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

}



// вывод матрицы
void Print_matrix(matrix<double> a)
{
    std::cout << std::endl;
    for(int i=0; i<a.size1(); i++)
    {
        for(int j=0; j<a.size2(); j++)
        {

            std::cout << a(i,j) << ' ';

        }
        std::cout << std::endl;
    }

}

// собираем матрицу СЛАУ для ур ия лапласа
matrix<double> Construct_matrix_Laplace()
{
    int row=0;
    matrix<double> a(n, n, .0);// матрица СЛАУ
    a.clear();
    for(int k=1; k<=qz; k++)
    {
        for(int j=1; j<=qy; j++)
        {
            for(int i=1; i<=qx; i++)
            {
                if (k>1)        a(row,row-qx*qy) = 1;
                if (j>1)        a(row,row-qx)   = 1;
                if (i>1)        a(row,row-1)   = 1;
                a(row,row) = -6;
                if (k<qz)        a(row,row+qx*qy) = 1;
                if (j<qy)        a(row,row+qx)   = 1;
                if (i<qx)        a(row,row+1)   = 1;

                //implementing Neumann B.C.
                if (k<qz && row<qx*qy)        a(row,row+qx*qy) = 2;
                row++;
            }
        }

    }
    return a;
}



// создаем нагрузку

vector<double> Construct_load_Laplace(double h, vector<double> bc)
{
    vector<double> b(n, .0);
    for(int i=0 ; i<n; i++)
    {
        if (i < qx*qy)    b(i) = 2*h*bc(i);
    }

    return b;
}

// задаем граничный условия

vector<double> Construct_BC_Laplace()
{
    vector<double> bc(n, .0);
    for(int i=0; i<n; i++)
    {
        bc(i) = .2;
    }

    return bc;
}


// задаем граничные условия

vector<double> Construct_BC_Poisson()
{
    vector<double> bc1(n, .0);
    for(int i=0; i<n; i++)
    {
        bc1(i) = .1;
    }

    return bc1;
}

// print unknowns
void Print_vectors(vector<double> u)
{
    int j=0, k=0;
    std::cout << std::endl;
    for(int i=0; i<u.size(); i++)
    {
        std::cout << u[i]<< " ";

        if((i + 1) % (u.size() / qx / qy) == 0) std::cout << std::endl;
        if((i + 1) % (u.size() / qz ) == 0) std::cout << std::endl;

    }
    std::cout << std::endl;
}

matrix<double> Construct_matrix_Poisson()
{
    int row=0;
    matrix<double> a(n, n, .0);// матрица СЛАУ
    a.clear();
    for(int k=1; k<=qz; k++)
    {
        for(int j=1; j<=qy; j++)
        {
            for(int i=1; i<=qx; i++) {

                if (k > 1) a(row, row - qx * qy) = 1;
                if (j > 1) a(row, row - qx) = 1;
                if (i > 1) a(row, row - 1) = 1;
                a(row, row) = -6;
                if (k < qz) a(row, row + qx * qy) = 1;
                if (j < qy) a(row, row + qx) = 1;
                if (i < qx) a(row, row + 1) = 1;
                row++;
            }
        }

    }
    return a;
}


vector<double> Construct_load_Poisson(double h, vector<double> bc, vector<double> f)
{
    vector<double> b(n, .0);
    for(unsigned long i=0 ; i<n; i++)
    {
        if (i < qx*qy)    b(i) = bc(i) + f(i)*2.0*h;

    }

    return b;
}



