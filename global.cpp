//
// Created by igor on 14.09.17.
//
#ifndef FDM_GLOBAL_CPP
#define FDM_GLOBAL_CPP

#include "/home/igor/Eigen/Eigen/SparseCore"
//размеры сетки - глоб константы

const long qx=300, qy=300, qz=50;
const long n = qx * qy * qz;
const double h = (double)qx/(qx-1.0), error=1.e-5; // шаг сетки и невязка
const double A = 10.0, B = 10.0;// полуоси эллипса раскрытия

typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;

#endif
