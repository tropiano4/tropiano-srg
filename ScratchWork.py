# Scratch Work

import numpy as np
from os import chdir,getcwd
from scipy.integrate import odeint
import time


class srg(object):
    
    def __init__(self,vnn,k_array,k_weights,generator,lamb):
        
        # Dimension of Hamiltonian
        self.N = len(vnn)
        
        # Build kinetic energy matrix
        # Coupled channel??
        ksq_array = np.concatenate((k_array**2,k_array**2))
        T_matrix = np.diag(ksq_array)
        
        # Convert vnn to fm^-2
        row,col = np.meshgrid(k_array*np.sqrt(k_weights),k_array*np.sqrt(k_weights))
        V0_matrix = vnn*2/np.pi*row*col
        
        # Save Hamiltonian and kinetic energy
        self.H0_matrix = T_matrix+V0_matrix
        self.T_matrix = T_matrix
        
    def commutator(self,A,B):
        
        return A@B-B@A
        
    def matrix2vector(self,A):
    
        # Dimension of matrix
        N = self.N
        # Dimension of vectorized matrix
        n = int(N*(N+1)/2)
        # Initialize vectorized matrix
        B = np.zeros(n)
    
        first = 0
        last = N
    
        for i in range(N):
        
            B[first:last] = A[i][i:]
            first = last
            last += N-i-1

        return B
 
    def vector2matrix(self,B):
        
        # Dimension of the vectorized matrix
        n = len(B)
        # Dimension of matrix (given by solving N*(N+1)/2 = n)
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
        return A+(A.T-np.diag(np.diag(A)))
    
    # Derivs
    
    
    
    # evolve_hamiltonian
    
    # main
        