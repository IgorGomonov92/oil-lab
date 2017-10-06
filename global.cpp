//
// Created by igor on 14.09.17.
//
#ifndef FDM_GLOBAL_CPP
#define FDM_GLOBAL_CPP

#include "/home/igor/Eigen/Eigen/SparseCore"
//размеры сетки - глоб константы

const long qx=100, qy=100, qz=100;
const long n = qx * qy * qz;
const double h=1.0/(qx-1.0), error=1.e-5; // шаг сетки и невязка
const double A = 10.0, B = 10.0;// полуоси эллипса раскрытия

typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;

#endif
