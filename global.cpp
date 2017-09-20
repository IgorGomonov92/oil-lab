//
// Created by igor on 14.09.17.
//
#include "/home/igor/Eigen/Eigen/SparseCore"

//размеры сетки - глоб константы

const long qx=4, qy=4, qz=2;
const long n = qx * qy * qz;
const double h=1.0; // шаг сетки

typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;
