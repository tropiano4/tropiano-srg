#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: srg_block_diagonal.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# Revision history:
#   May 28, 2019 --- Solve flow equation with respect to parameter lambda and 
#                    use SciPy's ode function.
# 
# Evolves Hamiltonian to block-diagonal, decoupled form with flow parameter 
# lambda [fm^-1] using the block-diagonal generator.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.integrate import ode
  

class SRG(object):
    
    
    def __init__(self, H_initial, lambda_bd, k_array, coupled_channel):
        """
        Saves the initial Hamiltonian in units fm^-2, the length of the 
        matrix. 
        
        Parameters
        ----------
        H_initial : 2-D ndarray
            Initial Hamiltonian matrix in units MeV.
        lambda_bd : float
            Cutoff for block-diagonal decoupling.
        k_array : 1-D ndarray
            Momentum array (used for creating the projection matrices).
        coupled_channel : bool
            True if the channel is coupled channel and false otherwise.
        
        """
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_M = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / hbar_sq_over_M
        
        # Save length of matrix
        self.N = len(H_initial)
        
        # Array of True and False values
        # In index notation: bool_array[i] = k_array[i] < lambda_bd[i]
        bool_array = k_array < lambda_bd
        
        # Length of momentum array
        n = len(k_array)
        
        # Build projection operators using bool_array
        
        # Matrix of ones along the diagonal up to k > lambda_bd
        p = np.diag( np.ones(n) * bool_array )
        
        # Opposite of p
        q = np.identity(n) - p
        
        # Projection operators for coupled channel potentials
        if coupled_channel:
            o = np.zeros((n,n))
            self.P = np.vstack( ( np.hstack( (p, o) ), np.hstack( (o, p) ) ) )
            self.Q = np.vstack( ( np.hstack( (q, o) ), np.hstack( (o, q) ) ) )

        else:
            self.P = p
            self.Q = q

    
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
    
    
    def derivs(self, lamb, H_evolved):
        """
        Right-hand side of the SRG flow equation using the block-diagonal 
        generator.
        
        Parameters
        ----------
        lamb : float
            Evolution parameter lambda in units fm^-1.
        H_evolved : 2-D ndarray
            Evolving Hamiltonian which is a vector and function of lambda. 
            Units are fm^-2.
        
        Returns
        -------
        dH_vector : 1-D ndarray
            Derivative with respect to lambda of the evolving Hamiltonian which 
            is a vector. Units are fm^-2.

        """
        
        # Matrix form of the evolving Hamiltonian
        H_matrix = self.vector2matrix(H_evolved)

        # Block-diagonal SRG generator, eta = [G,H] where G = H_BD 
        H_bd = self.P @ H_matrix @ self.P + self.Q @ H_matrix @ self.Q
        eta = self.commutator( H_bd, H_matrix)
            
        # RHS of flow equation in matrix form
        dH_matrix = -4.0 / lamb**5 * self.commutator(eta, H_matrix)
        
        # Returns vector form of RHS of flow equation
        dH_vector = self.matrix2vector(dH_matrix)
        
        return dH_vector


    def evolve_hamiltonian(self, lambda_initial, lambda_array):
        """
        Evolved Hamiltonian H_matrix at each value of lambda in lambda_array.
        
        Parameters
        ----------
        lambda_initial : float
            Initial value of lambda in units fm^-1. (Most potentials in
            Potentials/vsrg_macos are generated at lambda = 12 fm^-1 but check
            run_generate_vsrg_vlowk.pl to be sure.)
        lambda_array : 1-D ndarray
            Lambda evolution values in units fm^-1.
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each lambda value (e.g. d[1.5] returns the evolved
            Hamiltonian at lambda = 1.5 fm^-1).
            
        """
        
        # Set-up ODE
        
        # Initial Hamiltonian as a vector
        H_initial = self.matrix2vector(self.H_initial)

        # Use SciPy's ode function to solve flow equation
        solver = ode(self.derivative)
        # Following the example in Hergert:2016iju with modifications to nsteps
        # and error tolerances
        solver.set_integrator('vode', method='bdf', order=5, nsteps=100000, 
                              atol=1e-10, rtol=1e-10)
        # Set initial value of Hamiltonian at lambda = lambda_initial
        solver.set_initial_value(H_initial, lambda_initial)
    
        # Initialize dictionary
        d = {}
    
        # Loop over lambda values in lambda_array
        for lamb in lambda_array:
            
            # Solve ode up to lamb and store in dictionary
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