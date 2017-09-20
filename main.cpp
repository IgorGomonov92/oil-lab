
#include <iostream>
#include "functions.h"
#include <omp.h>

#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include "../../My_linalg/my_linalg.h"
#include <chrono>



using namespace Eigen;
using namespace std::chrono;

int main(int argc, char **argv)
{
    omp_set_num_threads(omp_get_max_threads());
    Eigen::setNbThreads(omp_get_max_threads());

    VectorXd u, u1; // вектора решений

    u = Solve_Laplace();
    u1 = Solve_Poissons();

   // std::cout<< u <<std::endl << "----------" <<std::endl << u1;
        return 0;
}
