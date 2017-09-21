
#include <iostream>
#include "global.cpp"
#include "functions_Laplace.h"
#include "functions_Poisson.h"
#include "functions_general.h"

#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include "../../My_linalg/my_linalg.h"
#include <chrono>


using namespace Eigen;
using namespace std::chrono;

int main(int argc, char **argv) {
    omp_set_num_threads(omp_get_max_threads());
    Eigen::setNbThreads(omp_get_max_threads());

    VectorXd uL;
    std::vector<VectorXd> uP; // вектора решений

    uL = Solve_Laplace();
    uP = Solve_Poissons();

    // std::cout<< u <<std::endl << "----------" <<std::endl << u1;
    return 0;
}
