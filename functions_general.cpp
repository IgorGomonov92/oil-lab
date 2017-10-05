//
// Created by igor on 22.09.17.
//
#include "functions_general.h"
#include "global.cpp"


using namespace std::chrono;
using namespace Eigen;


//--------------------------------------------------------------------------------------

void Construct_E( VectorXd * E)
{
    for (int i = 0; i < qz; ++i)
    {
        E->coeffRef(i) = 1.0e+10;
    }
}

//--------------------------------------------------------------------------------------

void Construct_v(VectorXd * v )
{
    for (int i = 0; i < qz; ++i)
    {
        v->coeffRef(i) = 0.3;
    }

}

//--------------------------------------------------------------------------------------

void Construct_lamda(VectorXd * lamda ,VectorXd * E,VectorXd * v)
{
    for (int i = 0; i < qz; ++i)
    {
        lamda->coeffRef(i) = E->coeff(i)*v->coeff(i) / (1.0+v->coeff(i)) / (1.0-2.0*v->coeff(i));
    }


}

//--------------------------------------------------------------------------------------

void Construct_G( VectorXd * G ,VectorXd * E,VectorXd * v  )
{
    for (int i = 0; i < qz; ++i)
    {
        G->coeffRef(i) = E->coeff(i) / 2.0 / (1.0+v->coeff(i));
    }

}

//--------------------------------------------------------------------------------------

