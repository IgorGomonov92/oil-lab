//
// Created by igor on 14.09.17.
//
#include "/home/igor/Eigen/Eigen/SparseCore"

//размеры сетки - глоб константы

const long qx=100, qy=100, qz=100;
const long n = qx * qy * qz;
const long nP = qx * qy * 2;
const double h=1.0, error=1.e-10; // шаг сетки и требуемая точность

typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;
