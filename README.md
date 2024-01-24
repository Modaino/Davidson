# Block davidson method implementation using the Armadillo library.

This method is an iterative (not exact) eigenvalue solver designed to discover the lowest eigenvector of a matrix. In (solid-state) physics and quantum chemistry, a critical challenge lies in determining the ground state (lowest eigen value) of a system (Hamiltonian). As the dimension of quantum systems grows exponentially with the number of degrees of freedom, the practical size of the Hamiltonian becomes so extensive that employing exact methods (or often even storing the Hamiltonian in memory) is infeasible. Iterative methods prove advantageous, conserving valuable computational resources and thereby offering a more effective approach to addressing this formidable problem. The method is closely related to the Arnoldi method (https://en.wikipedia.org/wiki/Arnoldi_iteration).

The algorithms is based on the article: https://joshuagoings.com/2013/08/23/davidsons-method/ implemented in C++ and extendinf it to complex matrices.

Code snipets used as a base to write the c++ code (python from the article) can be found in the python directory.

## Example output of main.exe (compiled from main.cpp)

"
Start of davidson
Algorithm converged after 3 iterations
smallest eig_val:5.16376e-12

"
