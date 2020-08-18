#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a seperate file with an
# extension _testv#.py where v# corresponds to the version number. For example,
# momentum_projection_operator_testv1.py. Use the revision history below to
# document when and why these files are created.
#
# Revision history:
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py in Old_codes.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Compare times on Magnus v. SRG evolution for a couple potentials.


import numpy as np
from scipy.integrate import solve_ivp
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import vnn
from run_magnus import run_magnus


class SRG(object):
    
    
    def __init__(self, H_initial):
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_m = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / hbar_sq_over_m
        
        # Save length of matrix
        self.N = len(H_initial)

    
    def commutator(self, A, B):
        
        return A @ B - B @ A
    
    
    def matrix2vector(self, A):
    
        # Length of matrix
        N = self.N
        # Length of vectorized matrix
        n = int( N * (N+1) / 2 )
        
        # Initialize vectorized matrix
        B = np.zeros(n)
    
        a = 0
        b = N
    
        for i in range(N):
        
            B[a:b] = A[i][i:]
            a = b
            b += N-i-1

        return B

 
    def vector2matrix(self, B):
        
        # Length of matrix (given by solving N*(N+1)/2 = n)
        N = self.N
    
        # Initialize matrix
        A = np.zeros( (N, N) )
    
        # Build upper half of A with diagonal

        a = 0
        b = N

        for i in range(N):

            A[i, i:] = B[a:b]
            a = b
            b += N-i-1

        # Reflect upper half to lower half to build full matrix
        # [np.transpose(A)-np.diag(np.diag(A))] is the lower half of A 
        # excluding the diagonal
        return A + ( np.transpose(A) - np.diag( np.diag(A) ) )
    
    
    def derivative(self, s, H_evolved):
        
        # Matrix form of the evolving Hamiltonian
        H_matrix = self.vector2matrix(H_evolved)

        # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
        # evolving Hamiltonian)
        eta = self.commutator( np.diag( np.diag(H_matrix) ), H_matrix)
            
        # RHS of flow equation in matrix form
        dH_matrix = self.commutator(eta, H_matrix)
        
        # Returns vector form of RHS of flow equation
        dH_vector = self.matrix2vector(dH_matrix)
        
        return dH_vector


    def evolve_hamiltonian(self, lambda_array):
        
        N = self.N

        # Set-up ODE
        
        # Initial Hamiltonian as a vector and dictionary
        H_initial = self.matrix2vector(self.H_initial)
        d = {}
        
        s_array = np.append( np.array( [0.0] ), 1.0 / lambda_array**4.0)
        
        sol = solve_ivp(self.derivative, s_array, H_initial, method='BDF')

        # # Store transformation matrix in dictionary for each lambda
        # i = 1 # Start at 1 to skip s = 0 transformation
        # for lamb in lambda_array:
            
        #     Hs_vector = sol[i]
        #     d[lamb] = np.reshape(Hs_vector, (N, N))
        #     i += 1

        # return d
        return sol
    
    
def run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array):

    # Load initial Hamiltonian, kinetic energy, momentum, and weights
    H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)

    # Initialize SRG class
    evolve = SRG(H_initial)
        
    # Time the evolution and return dictionary d of evolved Hamiltonians where
    # the keys are lambda values
    t0 = time.time() # Start time
    d = evolve.evolve_hamiltonian(lambda_array)
    t1 = time.time() # End time
    
    # Print details
    mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
    print('_'*85)
    print( 'H(s) done evolving to final lambda = %.2f fm^-1 after %f minutes'
          % (lambda_array[-1], mins) )
    print('_'*85)
    print('\nSpecifications:\n')
    print( 'kvnn = %d, channel = %s' % (kvnn, channel) )
    print( 'kmax = %.1f, kmid = %.1f, ntot = %d' % (kmax, kmid, ntot) )
    print( 'method = srg, generator = %s' % generator )

    # Otherwise, only return the dictionary d
    return d


if __name__ == '__main__':
    
    #kvnns = (111, 900)
    kvnns = [901]
    channel = '3S1'
    ntot = 120
    generator = 'Wegner'
    lambda_array = np.array( [6.0, 3.0, 2.0, 1.2] )
    methods = ['SRG', 'Magnus']
    k_magnus_values = (2, 10)
    
    for kvnn in kvnns:
        
        if kvnn == 111:
            
            kmax = 10.0
            kmid = 2.0
            
        else:
            
            kmax = 30.0
            kmid = 4.0
            
        for method in methods:
            
            if method == 'SRG':
                
                d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, 
                            lambda_array)
                
            else:
                
                for k_m in k_magnus_values:
                
                    d = run_magnus(kvnn, channel, kmax, kmid, ntot, generator, 
                                   lambda_array, k_magnus=k_m, ds=1e-6, 
                                   save=False)