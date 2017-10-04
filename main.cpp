
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
   // omp_set_num_threads(omp_get_max_threads());
   // Eigen::setNbThreads(omp_get_max_threads());


    VectorXd uL(n);// вектор  решений уравнения лапласа
    VectorXd uP(n);

    uL = Solve_Laplace();
    uP = Solve_Poisson( &uL );


    VectorXd w0(qx*qy), bcP(n);
    VectorXd E(qz), v(qz), lamda(qz), G(qz) ; // упругие параметры
    VectorXd w_Derivative_z(n);

    Construct_E(&E);
    Construct_v(&v);
    Construct_lamda(&lamda, &E, &v);
    Construct_G(&G, &E, &v);
    Construct_w_Derivative_z( &w_Derivative_z, &uP);
    Construct_w0(&w0);

    std::ofstream output ("output.txt");
    for (int k = 0; k < qz; ++k)
      for (int l = 0; l < qy; ++l)
        for (int m = 0; m < qx; ++m)
        {
            output << std::scientific << std::setprecision(5) << uP(k*qx*qy+l*qy+m) <<"  "<<lamda[0]* uL.coeff(k*qx*qy+l*qy+m) <<"  "<< w_Derivative_z(k*qx*qy+l*qy+m) <<"  " <<  lamda[0]*uL.coeff(k*qx*qy+l*qy+m) + 2*G[0]*  w_Derivative_z.coeff(k*qx*qy+l*qy+m)  << "   "<< m+1 << " " << l+1 << " " <<  k+1 << std::endl;
        }
    output.close();


    return 0;

}
