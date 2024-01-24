#include <iostream>
#include <armadillo>


int main(){
    std::cout << "Start Armadillo test" << std::endl;

    arma::cx_mat X(3, 6, arma::fill::randu);

    arma::cx_mat Q;
    arma::cx_mat R;

    arma::qr(Q, R, X.cols(0,3));

    X.print();
    std::cout << ".............." << std::endl;
    Q.print();
    std::cout << ".............." << std::endl;
    R.print();
    std::cout << ".............." << std::endl;
    (Q*R).print();
    return 0;
}
