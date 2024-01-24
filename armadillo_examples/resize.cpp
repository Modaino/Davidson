#include <iostream>
#include <armadillo>

int main(){
    arma::mat A = * new arma::mat(5, 5, arma::fill::randu);
    A = A + A.t();
    A.print();
    std::cout << "----------" << std::endl;
    arma::vec b(5, arma::fill::ones);

    //A.resize(5, 6);
    //A.col(6-1) = b;
    
    A.resize(6, 5);
    A.row(6-1) = b.t();
    A.print();
}