#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: srg_wegner_alt.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     January 27, 2020
# 
# Evolves Hamiltonian to band-diagonal, decoupled form with flow parameter 
# lambda [fm^-1] using the Wegner generator by solving the evolution equation
# for U(s).
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.integrate import odeint
  

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
        hbar_sq_over_m = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / hbar_sq_over_m
        
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
    
    
    def derivative(self, U_vector, s):
        """
        Right-hand side of the SRG flow equation using the Wegner generator.
        
        Parameters
        ----------
        U_vector : 1-D ndarray
            Unitary transformation which is a vector and function of s.
        s : float
            Evolution parameter s [fm^4].
        
        Returns
        -------
        dU_vector : 1-D ndarray
            Derivative with respect to s of the unitary transformation which 
            is a vector [fm^-2].

        """
        
        # Dimension of matrix
        N = self.N
        
        # Matrix form of the unitary transformation
        U_matrix = np.reshape(U_vector, (N, N))
        
        # Evolve the Hamiltonian to compute \eta(s)
        H_matrix = U_matrix @ self.H_initial @ U_matrix.T

        # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
        # evolving Hamiltonian)
        eta = self.commutator( np.diag( np.diag(H_matrix) ), H_matrix)
            
        # RHS of equation in matrix form
        dU_matrix = eta @ U_matrix
        
        # Returns vector form of RHS
        dU_vector = np.reshape(dU_matrix, -1)
        
        return dU_vector


    def solve(self, lambda_array):
        """
        Unitary transformation at each value of lambda in lambda_array.
        
        Parameters
        ----------
        lambda_array : 1-D ndarray
            Lambda evolution values [fm^-1].
            
        Returns
        -------
        d : dict
            Dictionary storing each transformation with keys (floats)
            corresponding to each lambda value (e.g. d[1.5] returns the 
            transformation at lambda = 1.5 fm^-1).
            
        """

        # Set-up ODE
        
        # Initialize dictionary
        d = {}
        
        # Dimension of matrices
        N = self.N
        # Initial transformation and convert to vector
        U_initial = np.eye(N).reshape(-1)
        
        # Convert from lambda values to s including initial value of 0
        s_array = np.append( np.array( [0.0] ), 1.0 / lambda_array**4.0)
        
        # Use SciPy's odeint function to solve flow equation
        # Each index of sol (corresponding to an index of s_array) gives the
        # vectorized solution of the ODE, that is, U_vector
        sol = odeint(self.derivative, U_initial, s_array)
        
        # Store transformation matrix in dictionary for each lambda
        i = 1 # Start at 1 to skip s = 0 transformation
        for lamb in lambda_array:
            
            U_vector = sol[i]
            d[lamb] = np.reshape(U_vector, (N, N))
            i += 1

        return d