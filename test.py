#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes.
#
# Last thing tested:
#   Compare H(s) from run_magnus.py dictionary and saved H(s) in
#   Potentials/vsrg_macos.
#
#------------------------------------------------------------------------------


import numpy as np
from run_magnus import run_magnus
from Potentials.vsrg_macos import load_save_potentials as lp

def vector2matrix(B,anti_hermitian=True):
    '''Takes the vector of a top right matrix and returns the full matrix. If 
    anti_hermitian = True, assumes B corresponds to a anti-hermitian matrix. If 
    false, assumes B corresponds to a hermitian matrix.'''
        
    # Dimension of the vectorized matrix
    n = len(B)
    # Dimension of matrix (given by solving N*(N+1)/2 = n for N where n = len(B))
    N = int( (-1 + np.sqrt(1+8*n) ) / 2 )
    # Initialize matrix
    A = np.zeros((N,N))
    
    # Build upper half of A with diagonal

    first = 0
    last = N

    for i in range(N):

        A[i,i:] = B[first:last]
        first = last
        last += N-i-1

    # Reflect upper half to lower half to build full matrix
    # [np.transpose(A)-np.diag(np.diag(A))] is the lower half of A excluding 
    # the diagonal
    if anti_hermitian:
        return A-(np.transpose(A)-np.diag(np.diag(A)))
    else: # Hermitian
        return A+(np.transpose(A)-np.diag(np.diag(A)))

# Set-up
kvnn = 901
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120

generator = 'Wegner'
lamb = 2.8
k_magnus = 6


d = run_magnus(kvnn, channel, kmax, kmid, ntot, generator, np.array([lamb]),
               k_magnus, ds=1e-5, save=False)
H_matrix_1 = d['hamiltonian'][lamb] * 41.47

H_matrix_2 = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'magnus', 
                                 generator, lamb, k_magnus=k_magnus, ds=1e-5)

print(H_matrix_1[:6, :6])
print(H_matrix_2[:6, :6])