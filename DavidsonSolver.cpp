#include "DavidsonSolver.hpp"

using namespace DavidsonSolvers;


DavidsonSolver::DavidsonSolver(const arma::cx_mat & A, double tolerance, size_t maximum_number_of_iteration){
    if (tolerance <= 0.0){ tolerance = arma::datum::eps; }
    tolerance_ = tolerance;
    dimension_ = A.n_rows;
    if (maximum_number_of_iteration < dimension_){
        maximum_number_of_iteration_ = maximum_number_of_iteration;
    }else{
        maximum_number_of_iteration_ = dimension_;
    }
    solve(A);
}

arma::cx_vec DavidsonSolver::project (const arma::cx_vec & u, const arma::cx_vec & v){
    arma::cx_vec *result = new arma::cx_vec(u);
    arma::cx_double scale = arma::dot(u, v) / arma::dot(u, v);          //moved outside of lambda for safty reasons:
    result->for_each([&](arma::cx_double & X) { return X * scale; } );  //(O:) can const reference be captured by reference in lambda?
    return *result;                                                     //dereferencing return value for ease of use
}

void DavidsonSolver::gram_schmidt_orthogonalize_in_place(arma::cx_mat &A){
    //Modified GS algorithm
    arma::cx_mat *P = new arma::cx_mat(arma::size(A), arma::fill::zeros);
    size_t iter = 0;
    arma::cx_vec previous_p = A.col(0);
    P->each_col( [&](arma::cx_vec & p){ 
        p = project(previous_p, A.col(iter));
        previous_p = p;
        iter++;
    } );
    A = A - *P;
}

bool DavidsonSolver::test_MGS(arma::cx_mat &A){
    gram_schmidt_orthogonalize_in_place(A);

    arma::cx_double zero_cp = 0;
    bool is_diagonal = true;
    arma::cx_vec previous = A.tail_rows(1);
    A.each_col( [&](arma::cx_vec &vec){
        if ( std::abs(arma::dot(vec, previous)) >  tolerance_ ){
            is_diagonal = false;
        }
    } );

    return is_diagonal;
}

void DavidsonSolver::get_smallest_eigenpair(const arma::cx_mat & A, arma::cx_double & lambda, arma::cx_vec & v){
    arma::cx_vec eig_vals;
    arma::cx_mat eig_vecs;
    arma::eig_gen(eig_vals, eig_vecs, A); //generating eigen pairs
    // size_t position = 0;
    // size_t iter = 0;
    // eig_vals.for_each( [&] (const arma::cx_vec::elem_type& value) { //finding minimal value & position
    //                                                                 if (std::abs(value) < std::abs(lambda)){
    //                                                                     position = iter;
    //                                                                     lambda = value;
    //                                                                     v = eig_vecs.row(position);
    //                                                                 } 
    //                                                                 iter++;
    //                                                               } );
    lambda = eig_vals.min();
    v = arma::normalise( eig_vecs.col(eig_vals.index_min()) );
}

void DavidsonSolver::solve(const arma::cx_mat & A){
    arma::cx_vec v_1(dimension_, arma::fill::zeros);              //initial vector
    v_1(0) = 1.0;                                                 //unit basis
    v_1 = arma::shuffle(v_1);                                     //randomise
    
    arma::cx_mat I(arma::size(A), arma::fill::eye);               //Identity for later use
    arma::cx_mat V(dimension_, 1, arma::fill::zeros);             //Guess vector(s)
    V.col(0) = v_1;

    arma::cx_double previous_eigen_value = 10;
    arma::cx_vec smallest_eigen_vec;
    arma::cx_double smallest_eigen_val;

    for (size_t k = 0; k < maximum_number_of_iteration_; k++)
    {
        arma::cx_mat H_k = V.t() * A * V; //Generating the Rayleigh matrix
        get_smallest_eigenpair(H_k, smallest_eigen_val, smallest_eigen_vec);
        //arma::cx_vec Ritz_vector = V.col(k) * smallest_eigen_vec;
        //arma::cx_vec residual = (smallest_eigen_val * I - A) * Ritz_vector;

        auto residual_value = std::abs(smallest_eigen_val - previous_eigen_value);
        if (residual_value <= tolerance_){
            smallest_eigen_value_ = std::abs(smallest_eigen_val);
            smallest_eigen_vector_ = smallest_eigen_vec;
            std::cout << "Algorithm converged after " << k << " iterations" << std::endl;
            return;
        }
        else{
            //Eigen value must be real, complex part is kept in bounds / forced to be zero
            //to increase stability (against complex vibrations)
            if (std::imag(smallest_eigen_val) >= tolerance_){
                smallest_eigen_val.imag(std::imag(smallest_eigen_val) / 2); //TODO: look at stability
            }

            arma::cx_vec Ritz_vector = V * smallest_eigen_vec;
            arma::cx_vec residual = (smallest_eigen_val * I - A) * Ritz_vector;
            V.resize(dimension_, k+2);
            // Method #1
            arma::cx_mat correction_matrix = ( smallest_eigen_val * I - arma::diagmat(A) );    //this is a diagonal matrix
            //correction_matrix = 0.5 * (correction_matrix + correction_matrix.t());           //forcing it to be real -> unstable
            arma::cx_mat inverse_correction_matrix = arma::diagmat(I / correction_matrix);
            arma::cx_vec new_direction = inverse_correction_matrix*residual;
            
            //Check if we accidentaly got an eigen value (not the smallest)
            new_direction.transform( [](arma::cx_double val){ return (std::isnan(val.real()) ? 0 : val); } );
            
            V.col(k+1) = new_direction;
            //Orthogonalize vector in V



            if ( std::abs(smallest_eigen_val) <= std::abs(previous_eigen_value) ){
                previous_eigen_value = smallest_eigen_val;
            }
        }
    }
    std::cout << "Algorithm did not converge" << std::endl;
}

double DavidsonSolver::get_smallest_eigen_value(){
    return smallest_eigen_value_;
}

arma::cx_vec DavidsonSolver::get_smallest_eigen_vector(){
    return smallest_eigen_vector_;
}

BlockDavidsonSolver::BlockDavidsonSolver(const arma::cx_mat & A,
                                        size_t number_of_eigenvalues,
                                        double tolerance,
                                        size_t maximum_number_of_iteration){
    number_of_eigenvalues_ = number_of_eigenvalues;
    auto l = 2*number_of_eigenvalues_;                       // number of initial guess vectors
    arma::cx_mat V(arma::size(A.n_rows, l), arma::fill::eye);// set of k unit vectors as guess
    arma::cx_mat I(arma::size(A), arma::fill::eye);          // identity matrix same dimen as A

    arma::cx_vec theta_old(number_of_eigenvalues, arma::fill::ones);
    arma::cx_vec theta;

    //Begin block Davidson routine
    for (size_t k = 0; k < maximum_number_of_iteration; k++)
    {
        arma::cx_mat W_k = A*V;
        arma::cx_mat Rayleigh_matrix = V.t()*W_k;
        arma::cx_vec eig_vals;
        arma::cx_mat eig_vecs;
        arma::eig_gen(eig_vals, eig_vecs, A);
        arma::cx_mat Ritz_vectors = V*eig_vecs;
        //arma::cx_mat residuals(arma::size(A.n_rows, l), arma::fill::zeros);
        //for (size_t i = 0; i < l; i++)
        //{
        //    residuals.row(i) = eigen_values(i) * Ritz_vectors.row(i) - W_k * eig_vecs.row(i);
        //}
        arma::cx_mat residuals = eig_vals * Ritz_vectors - W_k * eig_vecs;
        /*
        arma::cx_mat Q;
        arma::cx_mat R;
        arma::qr(Q, R, V.cols(0, m));
        V.cols(0, m) = Q;
        arma::cx_mat T = V.cols(0, m).t() * A * V.cols(0, m);
        arma::cx_vec THETA;
        arma::cx_mat S;
        arma::eig_gen(THETA, S, T, "balance");
        theta = arma::sort(THETA);
        arma::cx_mat s = arma::sort(S);
        for (size_t j = 0; j < k; j++)
        {
            arma::cx_mat w = (A - theta(j)*I) * (V.cols(0, m) * s.col(j));
            arma::cx_mat q = w / (theta(j)-A(j,j));
            V.cols(0, m+j) = q;
        }
        auto norm = arma::norm(theta.head(number_of_eigenvalues) - theta_old);
        if (norm < tolerance){
            break;
        }*/
    }
    
    //eigen_values = theta_old.head(number_of_eigenvalues);
    //eig_vecs = V.cols(0, number_of_eigenvalues);
}

BlockDavidsonSolver::BlockDavidsonSolver(){

}

std::vector<double> BlockDavidsonSolver::get_eigen_values(){
    std::vector<double> result;
    for (size_t i = 0; i < number_of_eigenvalues_; i++){
        result.push_back( arma::as_scalar( eigen_values(i).real() ));
    }
    return result;
}