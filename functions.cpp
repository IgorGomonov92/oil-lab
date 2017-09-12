//
// Created by gomonov on 11.09.17.
//
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <omp.h>

using namespace boost::numeric::ublas;

//Установить количество потоков

// реализация метода BiCGSTAB

vector<double> BiCGSTAB(matrix<double> a, vector<double> b) {
    vector<double> r, r0, u0, rt, p, p0, v, v0, prom1, prom2, s, t, u;
    double alpha0 = 1, alpha, beta = 0, ro0 = 1, ro, w, w0 = 1;

    u.clear();
    u0.clear();
    v0.clear();
    p0.clear();
    p.clear();
    p.resize(b.size(), true);
    u0.resize(b.size(), true);
    v0.resize(b.size(), true);
    p0.resize(b.size(), true);
    u.resize(b.size(), true);

    r0 = b - prod(a, u);
    rt = r0;

    do {
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

    return u;
}

double y=0.0;
// вывод матрицы
void Print_matrix(matrix<double> a)
{
    std::cout << std::endl;
    for(unsigned long i=0; i<a.size1(); i++)
    {
        for(unsigned long j=0; j<a.size2(); j++)
        {

            std::cout << a(i,j) << ' ';
            y = a(i,j);
        }
        std::cout << std::endl;
    }

}

// собираем матрицу СЛАУ
matrix<double> Construct_matrix_Laplace(int q)
{
    unsigned long row=0;
    unsigned long n;
    n = q * q * q;
    matrix<double> a(n, n);// матрица СЛАУ
    a.clear();
    for(unsigned long k=1; k<=q; k++)
    {
        for(unsigned long j=1; j<=q; j++)
        {
            for(unsigned long i=1; i<=q; i++)
            {
                if (k>1)        a(row,row-q*q) = 1;
                if (j>1)        a(row,row-q)   = 1;
                if (i>1)        a(row,row-1)   = 1;
                a(row,row) = -6;
                if (k<q)        a(row,row+q*q) = 1;
                if (j<q)        a(row,row+q)   = 1;
                if (i<q)        a(row,row+1)   = 1;

                //implementing Neumann B.C.
                if (k<q && row<q*q)        a(row,row+q*q) = 2;
                row++;
            }
        }

    }
    return a;
}



// создаем нагрузку

vector<double> Construct_load(int q, double h, vector<double> bc)
{
    unsigned long n=q*q*q;
    vector<double> b(n);
    unsigned long j=0;
    for(unsigned long i=0 ; i<n; i++)
    {
        if (i < q)    b(i) = 2*h*bc(j);
        j++;
    }

    return b;
}

// задаем граничный условия

vector<double> Construct_BC(int q)
{
    unsigned long n;
    n = q * q * q;
    vector<double> bc(n);
    for(unsigned long i=0; i<n; i++)
    {
        bc(i) = .2;
    }

    return bc;
}

// print unknowns
void Print_vectors(int q, vector<double> u)
{
    unsigned long j=0, k=0;
    std::cout << std::endl;
    for(unsigned long i=0; i<u.size(); i++)
    {
        std::cout << u[i]<< " ";

        if((i + 1) % (u.size() / q / q) == 0) std::cout << std::endl;
        if((i + 1) % (u.size() / q ) == 0) std::cout << std::endl;

    }
    std::cout << std::endl;
}





