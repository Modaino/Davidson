#!/bin/python
from __future__ import division
from __future__ import print_function
import math
import numpy as np
import time

''' Block Davidson, Joshua Goings (2013)

    Block Davidson method for finding the first few
	lowest eigenvalues of a large, diagonally dominant,
    sparse Hermitian matrix (e.g. Hamiltonian)

    Revised by Áron Vízkeleti to complex
'''

n = 120					# Dimension of matrix
tol = 1e-8				# Convergence tolerance
mmax = n//2				# Maximum number of iterations	

''' Create sparse, diagonally dominant matrix A with 
	diagonal containing 1,2,3,...n. The eigenvalues
    should be very close to these values. You can 
    change the sparsity. A smaller number for sparsity
    increases the diagonal dominance. Larger values
    (e.g. sparsity = 1) create a dense matrix
'''

sparsity = 0.0001
A = np.zeros((n,n), dtype=complex)
for i in range(0,n):
    A[i,i] = i + 1 + 1j
A = A + sparsity*np.random.randn(n,n) 
A = (np.conjugate(A.T) + A)/2 


k = 10					# number of initial guess vectors 
eig = 6			# number of eignvalues to solve 
t = np.eye(n,k, dtype=complex)			# set of k unit vectors as guess
V = np.zeros((n,n), dtype=complex)		# array of zeros to hold guess vec
I = np.eye(n, dtype=complex)			# identity matrix same dimen as A

# Begin block Davidson routine

start_davidson = time.time()

for m in range(k,mmax,k):
    if m <= k:
        for j in range(0,k):
            V[:,j] = t[:,j]/np.linalg.norm(t[:,j])
        theta_old = 1 
    elif m > k:
        theta_old = THETA[:eig]
    V[:,:m],R = np.linalg.qr(V[:,:m])
    T = np.dot(np.conjugate(V[:,:m].T),np.dot(A,V[:,:m]))
    THETA,S = np.linalg.eig(T)
    idx = THETA.argsort()
    for j in range(0,k):
        w = (A - THETA[j]*I)@(V[:,:m]@S[:,j]) 
        q = w/(THETA[j]-A[j,j])
        V[:,(m+j)] = q
    norm = np.linalg.norm(THETA[:eig] - theta_old)
    if norm < tol:
        break

end_davidson = time.time()

# End of block Davidson. Print results.

print("davidson = ", THETA[:eig],";",
    end_davidson - start_davidson, "seconds")

# Begin Numpy diagonalization of A

start_numpy = time.time()

E,Vec = np.linalg.eig(A)
E = np.sort(E)

end_numpy = time.time()

# End of Numpy diagonalization. Print results.

print("numpy = ", E[:eig],";",
     end_numpy - start_numpy, "seconds") 