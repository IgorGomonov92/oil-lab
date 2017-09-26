
#include <iostream>
#include "functions_Laplace.h"
#include "functions_Poisson.h"
#include "functions_general.h"
#include "global.cpp"
#include <omp.h>
#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include <chrono>
#include <fstream>
#include <iomanip>

using namespace Eigen;
using namespace std::chrono;

int main(int argc, char **argv) {
    omp_set_num_threads(omp_get_max_threads());
    Eigen::setNbThreads(omp_get_max_threads());
    std::ofstream outputW ("W.txt");


    VectorXd uL;// вектора решений
    std::vector<VectorXd> uP(qz+1), uPseparated(qz+1); // Вектор решения одного слоя уравнения пуассона содержит два слоя из=за необх учитывать ГУ Дирихле

    uL = Solve_Laplace();
    uP = Solve_Poissons( &uL );

    for (int i = 0; i < qz; ++i)
    {
        uPseparated[i].resize(qx*qy);
        uPseparated[i].fill(0);
        for (int j = 0; j <qx*qy ; ++j)
        {
            uPseparated[i].coeffRef(j) = uP[i].coeff(j);
        }
    }

     //std::cout<< uL <<std::endl << "----------";

    Construct_w_Derivative_z( &uPseparated);

   // for (int k = 0; k < qz; ++k)

        for (int l = 0; l < qy; ++l)

            for (int m = 0; m < qx; ++m)
            {
                outputW<< std::scientific << std::setprecision(5) <<uPseparated[10].coeff(l*qy+m)<< std::setw(10) <<" "<< m+1 << " " << l+1 << " " << 1 << std::endl;

            }



    return 0;

}
