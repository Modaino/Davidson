/* Block davidson algorithm implemented by Áron Vízkeleti, 
 *         Wigner Research institute 2021-05-06
 * Description:
 * 
 */
#ifndef BLOCKDAVIDSON_H
#define BLOCKDAVIDSON_H

#include <armadillo>
#include <vector>

namespace DavidsonSolvers
{

    class BlockDavidsonSolver{
        private:
            arma::cx_vec eigen_values;
            arma::cx_mat eigen_vectors;
            size_t number_of_eigenvalues_;
        public:
            BlockDavidsonSolver(const arma::cx_mat & A,
                                size_t number_of_eigenvalues, 
                                double tolerance,
                                size_t maximum_number_of_iteration);
            BlockDavidsonSolver();
            std::vector<double> get_eigen_values();
    };

    class DavidsonSolver{
        private:
            double tolerance_;
            size_t dimension_;
            size_t maximum_number_of_iteration_;
            double smallest_eigen_value_;
            arma::cx_vec smallest_eigen_vector_;
            arma::cx_vec project (const arma::cx_vec & u, const arma::cx_vec & v);
            void get_smallest_eigenpair(const arma::cx_mat & A, arma::cx_double & lambda, arma::cx_vec & v);
            void solve(const arma::cx_mat & A);
            void gram_schmidt_orthogonalize_in_place(arma::cx_mat & A);
        public:
            DavidsonSolver(const arma::cx_mat & A,
                           double tolerance,
                           size_t maximum_number_of_iteration);
            double get_smallest_eigen_value();
            arma::cx_vec get_smallest_eigen_vector();
            bool test_MGS(arma::cx_mat & A);
    };
}

#endif