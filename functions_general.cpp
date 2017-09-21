//
// Created by igor on 21.09.17.
//

#include "functions_Laplace.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>

#include "/home/igor/My_linalg/my_linalg.cpp"
#include "/home/igor/Eigen/Eigen/SparseCore"
#include <chrono>

using namespace std::chrono;
using namespace Eigen;

void Construct_f(std::vector<VectorXd> *f) {
    f->at(0).resize(nP);

    for (int i = 0; i < nP; i++) {
        f->at(0).coeffRef(i) = 1.0 / (i + 1);
    }

    for (int i = 1; i <= qz; i++) {
        f->at(i).resize(nP);

        for (int j = 0; j < nP; j++) {
            f->at(i).coeffRef(j) = 1.0 / (j + 1);
        }
    }

}