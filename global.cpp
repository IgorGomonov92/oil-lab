//
// Created by igor on 14.09.17.
//
#include "/home/igor/Eigen/Eigen/SparseCore"

//размеры сетки - глоб константы

const long qx=10, qy=10, qz=10;
const long n = qx * qy * qz;

typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;
