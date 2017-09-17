//
// Created by igor on 14.09.17.
//
#include "/home/igor/Eigen/Eigen/SparseCore"

//размеры сетки - глоб константы

const long qx=20, qy=20, qz=20;
const long n = qx * qy * qz;

typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<double> T;
