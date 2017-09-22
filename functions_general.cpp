//
// Created by igor on 22.09.17.
//
#include "functions_general.h"
#include "global.cpp"


using namespace std::chrono;
using namespace Eigen;


//--------------------------------------------------------------------------------------

void Construct_E( std::vector<double> * E)
{
    for (int i = 0; i < qz; ++i)
    {
        E->at(i) = i;
    }
}

//--------------------------------------------------------------------------------------

void Construct_v(std::vector<double> * v )
{
    for (int i = 0; i < qz; ++i)
    {
        v->at(i) = i;
    }

}

//--------------------------------------------------------------------------------------

void Construct_lamda(std::vector<double> * lamda ,std::vector<double> * E,std::vector<double> * v)
{
    for (int i = 0; i < qz; ++i)
    {
        lamda->at(i) = E->at(i)*v->at(i) / (1+v->at(i)) / (1-2*v->at(i));
    }

}

//--------------------------------------------------------------------------------------

void Construct_G( std::vector<double> * G ,std::vector<double> * E,std::vector<double> * v  )
{
    for (int i = 0; i < qz; ++i)
    {
        G->at(i) = E->at(i) / 2 / (1+v->at(i));
    }

}

//--------------------------------------------------------------------------------------

