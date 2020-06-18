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
#   Test error tolerances of SRG.


import numpy as np
import numpy.linalg as la
from scipy.integrate import ode
import time
from Potentials.vsrg_macos import load_save_potentials as lsp
  

class SRG(object):
    
    
    def __init__(self, H_initial):
        """
        Saves the initial Hamiltonian in units fm^-2 and the length of the 
        matrix.
        
        Parameters
        ----------
        H_initial : 2-D ndarray
            Initial Hamiltonian matrix [MeV].
        
        """
        
        # h-bar^2 / M [MeV fm^2]
        self.hbar_sq_over_M = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / self.hbar_sq_over_M
        
        # Save length of matrix
        self.N = len(H_initial)

    
    def commutator(self, A, B):
        """
        Commutator of A and B, [A,B] where A and B are square matrices.
        
        Parameters
        ----------
        A : 2-D ndarray
            First input square matrix.
        B : 2-D ndarray
            Second input square matrix.
            
        Returns
        -------
        out : 2-D ndarray
            Commutator of the two input matrices.
            
        """
        
        return A @ B - B @ A
    
    
    def matrix2vector(self, A):
        """
        Takes the upper triangle of the matrix A (including the diagonal) 
        and reshapes it into a vector B of length N*(N+1)/2.
        
        Parameters
        ----------
        A : 2-D ndarray
            Input matrix.
            
        Returns
        -------
        B : 1-D ndarray
            Output vector.
            
        """
    
        # Length of matrix
        N = self.N
        # Length of vectorized matrix
        n = int(N*(N+1)/2)
        
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
        """
        Takes the vector of a upper triangle matrix and returns the full 
        matrix. Use only for symmetric matrices.
        
        Parameters
        ----------
        B : 1-D ndarray
            Input vector.
        
        Returns
        -------
        A : 2-D ndarray
            Output matrix.
            
        """
        
        # Length of matrix (given by solving N*(N+1)/2 = n)
        N = self.N
    
        # Initialize matrix
        A = np.zeros((N, N))
    
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
    
    
    def derivative(self, lamb, H_evolved):
        """
        Right-hand side of the SRG flow equation using the Wegner generator.
        
        Parameters
        ----------
        lamb : float
            Evolution parameter lambda [fm^-1].
        H_evolved : 1-D ndarray
            Evolving Hamiltonian which is a vector and function of lambda 
            [fm^-2].
        
        Returns
        -------
        dH_vector : 1-D ndarray
            Derivative with respect to lambda of the evolving Hamiltonian which 
            is a vector [fm^-2].

        """
        
        # Matrix form of the evolving Hamiltonian
        H_matrix = self.vector2matrix(H_evolved)

        # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
        # evolving Hamiltonian)
        eta = self.commutator( np.diag( np.diag(H_matrix) ), H_matrix)
            
        # RHS of flow equation in matrix form
        dH_matrix = -4.0 / lamb**5 * self.commutator(eta, H_matrix)
        
        # Returns vector form of RHS of flow equation
        dH_vector = self.matrix2vector(dH_matrix)
        
        return dH_vector


    def evolve_hamiltonian(self, lamb, abs_tol, rel_tol):
        """
        Evolved Hamiltonian at each value of lambda in lambda_array.
        
        Parameters
        ----------
        lambda_array : 1-D ndarray
            Lambda evolution values [fm^-1].
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each lambda value (e.g. d[1.5] returns the evolved
            Hamiltonian at lambda = 1.5 fm^-1).
            
        """
        
        # Initial Hamiltonian as a vector and dictionary
        H_initial = self.matrix2vector(self.H_initial)

        # Set-up ODE
        
        # Use SciPy's ode function to solve flow equation
        solver = ode(self.derivative)
        # Following the example in Hergert:2016iju with modifications to nsteps
        # and error tolerances
        solver.set_integrator('vode', method='bdf', order=5, atol=abs_tol,
                              rtol=rel_tol, nsteps=5000000)
        # Set initial value of Hamiltonian at lambda = lambda_initial
        solver.set_initial_value(H_initial, 12.0)
    
        # Loop over lambda values in lambda_array
        # Solve ODE up to lamb and store in dictionary
        while solver.successful() and solver.t > lamb:
            
            # Select step-size depending on extent of evolution
            if solver.t >= 6.0:
                dlamb = 1.0
            elif solver.t < 6.0 and solver.t >= 2.5:
                dlamb = 0.5
            elif solver.t < 2.5 and solver.t >= lamb:
                dlamb = 0.1
                
            # This if statement prevents the solver from over-shooting 
            # lambda and takes a step in lambda equal to the exact amount
            # necessary to reach the specified lambda value
            if solver.t - dlamb < lamb:
                
                dlamb = solver.t - lamb
                
            # Integrate to next step in lambda
            H_vector = solver.integrate(solver.t - dlamb)
                
        H_matrix = self.vector2matrix(H_vector) * self.hbar_sq_over_M

        # Return in MeV
        return H_matrix
    
    
# --- Main program --- #


# Potential specifications
kvnn = 901
#kvnn = 111
channel = '3S1'

# Initial Hamiltonian, momentum, and weights
H_initial = lsp.load_hamiltonian(kvnn, channel)

# Calculate initial eigenvalues
if kvnn == 901:
    deuteron_index = 1
else:
    deuteron_index = 0
eigenvalues_init, eigenvectors_init = la.eig(H_initial)
deuteron_exact = np.sort(eigenvalues_init)[deuteron_index]
print('Deuteron exact = %.3f MeV\n' % deuteron_exact)

# SRG specifications
generator = 'Wegner'
lamb = 1.2
srg = SRG(H_initial)

for tol in [1e-3, 1e-6, 1e-9, 1e-12]:

    t0 = time.time() # Start time
    H_evolved = srg.evolve_hamiltonian(lamb, tol, tol)
    t1 = time.time() # End time
    
    # Print details
    secs = round(t1 - t0, 3) # Seconds elapsed evolving H(s)
    print('_'*85)
    print( 'H(s) done evolving to final lambda = %.1f fm^-1 after %.3f seconds'
          % (lamb, secs) )
    print('_'*85)
    print( 'Error tolerances = %.1e' % tol )
    
    eigenvalues_evol, eigenvectors_evol = la.eig(H_evolved)
    deuteron_lamb = np.sort(eigenvalues_evol)[deuteron_index]
    
    deuteron_error = abs((deuteron_lamb-deuteron_exact)/deuteron_exact)
    print('Deuteron error = %.5e'%deuteron_error)
    print('_'*85)