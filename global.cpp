//
// Created by igor on 14.09.17.
//
#ifndef FDM_GLOBAL_CPP
#define FDM_GLOBAL_CPP

#include "/home/igor/Eigen/Eigen/SparseCore"
//размеры сетки - глоб константы

const long qx=10, qy=10, qz=4;
const long n = qx * qy * qz;
const long nP = qx * qy * 2;
const double h=1.0/qx, error=1.e-20; // шаг сетки и требуемая точность

typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;

#endif
