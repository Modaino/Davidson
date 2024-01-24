#include <iostream>
#include <armadillo>


int main(){
    std::cout << "Start Armadillo test" << std::endl;

    arma::cx_mat A = * new arma::cx_mat(3, 3, arma::fill::randu);
    arma::cx_mat mI(arma::size(A), arma::fill::eye); 
    arma::cx_mat X = ( 1.5 * mI - arma::diagmat(A) ); 

    X.print();
    std::cout << ".............." << std::endl;
    arma::cx_mat Y = arma::diagmat( mI/X );
    (Y).print();
    std::cout << ".............." << std::endl;

    arma::cx_vec b(5, arma::fill::randn);



    return 0;
}