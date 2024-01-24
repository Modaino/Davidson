/*                                                  *
 *       C\C++ library for the Davidson method      *
 *                                                  *
 *   Last modified by Aron Vizkeleti on 2021-05-19  *
 *                                                  *
 *                                                  */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <cstdlib>
#include "DavidsonSolver.hpp"

using namespace DavidsonSolvers;

int main(){
    std::cout << "Start of davidson" << std::endl;
    arma::cx_mat A = * new arma::cx_mat(10, 10, arma::fill::randu);
    A = A + A.t();
    A = A * 0.01 + 100 * arma::diagmat(A); //Make it diagonally dominant 
    
    // std::ofstream myfile;
    // myfile.open("matrix.txt");

    // A.print(myfile);
    // myfile.close();

    DavidsonSolver mySolver = *new DavidsonSolver(A, 0.001, 2000);


    std::cout << "smallest eig_val:" << mySolver.get_smallest_eigen_value() << std::endl;
    
    delete &mySolver;
    std::exit(EXIT_SUCCESS);
}
