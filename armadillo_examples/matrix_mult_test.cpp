#include <iostream>
#include <armadillo>

/*
 *Armadillo example written by Áron Vízkeleti
 */

int main(){
    std::cout << "Start of davidson" << std::endl;

    arma::mat A(1, 5, arma::fill::randu);
    arma::vec b(5, arma::fill::randu);
    
    A.print();
    std::cout << ".............." << std::endl;
    b.print();
    std::cout << ".............." << std::endl;
    (A*b).print();
    return 0;
}