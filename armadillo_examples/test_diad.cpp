#include <iostream>
#include <armadillo>


int main(){
    std::cout << "Start Armadillo test" << std::endl;

    arma::mat A(5, 5, arma::fill::randu);
    arma::mat V(5, 2, arma::fill::randu);
    //arma::vec b(5, arma::fill::randu);
    
    V.print();
    std::cout << "Results .............." << std::endl;

    try{    (V.t()*A*V).print();
    } catch(std::exception e){        std::cout << e.what() << std::endl;    }
    std::cout << ".............." << std::endl;
    try{    (V*A*V.t()).print();
    } catch(std::exception e){        std::cout << e.what() << std::endl;    }
    return 0;
}