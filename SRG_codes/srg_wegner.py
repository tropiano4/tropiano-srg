#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: srg_wegner.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# Evolves Hamiltonian to band-diagonal, decoupled form with flow parameter 
# lambda [fm^-1] using the Wegner generator.
#
# Revision history:
#   05/27/19 --- Solve flow equation with respect to parameter lambda and use 
#                SciPy's ode function.
#   09/18/19 --- Added odeint as an option for the ODE solver instead of ode.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.integrate import ode
  

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
        hbar_sq_over_M = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / hbar_sq_over_M
        
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


    def evolve_hamiltonian(self, lambda_initial, lambda_array, method='ode'):
        """
        Evolved Hamiltonian at each value of lambda in lambda_array.
        
        Parameters
        ----------
        lambda_initial : float
            Initial value of lambda [fm^-1].
        lambda_array : 1-D ndarray
            Lambda evolution values [fm^-1].
        method : str, optional
            ODE solver to use: SciPy 'ode' or 'odeint'. 'ode' is the default.
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each lambda value (e.g. d[1.5] returns the evolved
            Hamiltonian at lambda = 1.5 fm^-1).
            
        """

        # Set-up ODE
        
        # Initial Hamiltonian as a vector and dictionary
        H_initial = self.matrix2vector(self.H_initial)
        d = {}
        
        # Use SciPy's ode function to solve flow equation
        solver = ode(self.derivative)
        # Following the example in Hergert:2016iju with modifications to nsteps
        # and error tolerances
        solver.set_integrator('vode', method='bdf', order=5, atol=1e-6,
                              rtol=1e-6, nsteps=5000000)
        # Set initial value of Hamiltonian at lambda = lambda_initial
        solver.set_initial_value(H_initial, lambda_initial)
    
        # Loop over lambda values in lambda_array
        for lamb in lambda_array:
            
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
                H_evolved = solver.integrate(solver.t - dlamb)
                
            # Store evolved Hamiltonian matrix in dictionary
            d[lamb] = self.vector2matrix(H_evolved)

        return d