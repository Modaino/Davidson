#include <iostream>
#include <armadillo>

int main(){
    //To find an orthogonal vector to a set of vectors A = [v_1, v_2, ...], solve A^T x = 0 for x (under determined, a constraint is required)
    arma::mat A = * new arma::mat(5, 5, arma::fill::randu);
    A.resize(5, 4);
    A.print();
    std::cout << "----------" << std::endl;
    arma::vec b(A.n_cols, arma::fill::zeros);
    b(0) = 1; //constraint
    b = arma::shuffle(b);
    b.print();
    std::cout << "----------" << std::endl;
    

    arma::vec x = arma::solve(A.t(), b, arma::solve_opts::fast);
    x.print();
}