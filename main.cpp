
#include <iostream>
#include "functions_Laplace.h"
#include "functions_Poisson.h"
#include "functions_general.h"
#include "global.cpp"
#include <omp.h>
#include "/home/igor/Eigen/Eigen/SparseCore"
#include </home/igor/Eigen/Eigen/IterativeLinearSolvers>
#include <chrono>


using namespace Eigen;
using namespace std::chrono;

int main(int argc, char **argv) {
    omp_set_num_threads(omp_get_max_threads());
    Eigen::setNbThreads(omp_get_max_threads());

    std::vector<double> E(qz), v(qz), lamda(qz), G(qz); // упругие параметры
    // заполняем вектотора упругих параметров в разных слоях соотв функциями
    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);

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

    // std::cout<< u <<std::endl << "----------" <<std::endl << u1;
    return 0;
}
