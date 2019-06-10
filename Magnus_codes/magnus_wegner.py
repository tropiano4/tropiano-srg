#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: magnus_wegner.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     June 3, 2019
# 
# Evolves Hamiltonian to band-diagonal, decoupled form with flow parameter s 
# using the Wegner generator and the Magnus implementation.
#
# Revision history:
#   June 5, 2019 --- Added an Euler method function.
#   June 7, 2019 --- Adding SciPy's expm function as an option instead of BCH
#
# Notes:
#   * Tried to do lambda differential equation similar to SRG codes but kept
#     getting infinity errors in computing omega matrix. Thus, we use the flow
#     parameter s, which has worked before.
#
#------------------------------------------------------------------------------


from math import factorial
import numpy as np
from scipy.linalg import expm
from sympy import bernoulli


class Magnus(object):
    
    
    def __init__(self, H_initial, k_magnus, ds=1e-5):
        """
        Saves the initial Hamiltonian in units fm^-2 and the length of the 
        matrix. Also saves the number of terms in the Magnus omega derivative
        equation, k_magnus, and the Euler method step-size ds. Initializes
        arrays for factorials and Bernoulli numbers.
        
        Parameters
        ----------
        H_initial : 2-D ndarray
            Initial Hamiltonian matrix in units MeV.
        k_magnus : int
            Number of terms to include in Magnus sum (that is,
            dOmega / ds ~ \sum_0^k_magnus ... )
        ds : float, optional
            Step-size in the flow parameter s.
        
        """
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_M = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / hbar_sq_over_M
        
        # Save length of matrix
        self.N = len(H_initial)
        
        # Save k_magnus
        self.k_magnus = k_magnus
        
        # Save step-size ds
        self.ds = ds
        
        # Initialize factorial and Bernoulli number arrays for summations up to
        # 30 terms
        self.factorial_array = np.zeros(31)
        self.magnus_factors = np.zeros(31)
        for i in range(31):
            self.factorial_array[i] = factorial(i)
            self.magnus_factors[i] = bernoulli(i) / factorial(i)

    
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
    
    
    def bch_formula(self, M_initial, O_evolved, k_truncate):
        """
        Evolved operator given an initial operator M and evolved Magnus omega
        matrix.

        Parameters
        ----------
        M_initial : 2-D ndarray
            Initial operator which is a matrix.
        O_evolved : 2-D ndarray
            Evolved omega matrix.
        k_truncate : int
            Truncation in BCH sum (this function includes the term k_truncate).
        
        Returns
        -------
        M_evolved : 2-D ndarray
            Magnus evolved operator as a matrix.
            
        """
        
        # Initial nested commutator ad_Omega^0(M)
        ad = M_initial
        
        # Load factorials
        factorial_array = self.factorial_array
        
        # k = 0 term
        M_evolved = ad / factorial_array[0]
        
        # Sum from k = 1 to k = k_truncate
        for k in range(1, k_truncate + 1):
            
            ad = self.commutator(O_evolved, ad)
            M_evolved += ad / factorial_array[k]
            
        return M_evolved

    
    def derivative(self, O_evolved):
        """
        Right-hand side of the Magnus derivative omega equation using the 
        Wegner generator.
        
        Parameters
        ----------
        O_evolved : 2-D ndarray
            Evolving omega matrix which is a function of s.
        
        Returns
        -------
        dO_matrix : 2-D ndarray
            Derivative with respect to s of the evolving omega matrix. 

        """
        
        # Compute the evolving Hamiltonian with the BCH formula
        #H_evolved = self.bch_formula(self.H_initial, O_evolved, 25)
        # Use expm instead
        H_evolved = expm(O_evolved) @ self.H_initial @ expm(-O_evolved)

        # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
        # evolving Hamiltonian)
        eta = self.commutator( np.diag( np.diag(H_evolved) ), H_evolved)
        
        # Load Magnus factors in sum
        magnus_factors = self.magnus_factors
        
        # Initial nested commutator ad_Omega^0(eta)
        ad = eta
        
        # k = 0 term
        dO_matrix = magnus_factors[0] * ad
        
        # Sum from k = 1 to k = k_magnus
        for k in range(1, self.k_magnus + 1):
            
            ad = self.commutator(O_evolved, ad)
            dO_matrix += magnus_factors[k] * ad
    
        return dO_matrix
    
    
    def euler_method(self, O_initial, s_init, s_final):
        """
        Solves the Magnus derivative omega equation using the first-order Euler
        method.
        
        Parameters
        ----------
        O_initial : 2-D ndarray
            Initial omega matrix which is a function of s.
        s_init : float
            Initial s value.
        s_final : float
            Final s value.
            
        Returns
        -------
        O_evolved : 2-D ndarray
            Evolved omega matrix.
            
        """
        
        # Load Euler method step-size
        ds = self.ds
        
        # Set initial s value and initial value of omega
        s = s_init
        O_evolved = O_initial
        
        # Step in s until s_final is reached
        while s <= s_final:

            # Next step in s
            O_evolved += self.derivative(O_evolved) * ds
                
            # Step to next s value
            s += ds
                
        # To ensure omega stops at s = s_final step-size is fixed to go from 
        # last value of s < s_final to s_final where 
        # ds_exact = s_final - (s - ds)
        ds_exact = s_final - s + ds
        O_evolved += self.derivative(O_evolved) * ds_exact
        
        return O_evolved
    

    def evolve_hamiltonian(self, lambda_array):
        """
        Magnus evolved Hamiltonian and omega matrix at each value of lambda in
        lambda_array using the first-order Euler method.
        
        Parameters
        ----------
        lambda_array : 1-D ndarray
            Lambda evolution values in units fm^-1.
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian and omega matrix with
            keys (floats) corresponding to each lambda value (e.g. 
            d['hamiltonian'][1.5] and d['omega'][1.5] returns the evolved 
            Hamiltonian and omega matrix at lambda = 1.5 fm^-1, respectively).
            
        """

        # Load initial Hamiltonian and length
        H_initial = self.H_initial
        N = self.N

        # Set-up ODE
        
        # Initial Magnus omega matrix
        O_initial = np.zeros( (N, N) )
    
        # Initialize dictionary
        d = {}
        # Key for evolved Hamiltonians
        d['hamiltonian'] = {}
        # Key for evolved omega matrices
        d['omega'] = {}
        
        # Initial s value, step-size, and omega matrix
        s_init = 0.0
        O_evolved = O_initial
    
        # Loop over lambda values in lambda_array
        for lamb in lambda_array:
            
            # Convert lamb to s value
            s_val = 1.0 / lamb**4.0
            
            # Solve ODE up to s_val using the Euler method and store in 
            # dictionary
            O_evolved = self.euler_method(O_evolved, s_init, s_val)
                
            # Store evolved omega matrix in dictionary
            d['omega'][lamb] = O_evolved
                
            # Evaluate the evolved Hamiltonian matrix using the BCH formula
            #H_evolved = self.bch_formula(H_initial, O_evolved, 25)
            # Use expm instead
            H_evolved = expm(O_evolved) @ H_initial @ expm(-O_evolved)
            
            # Store evolved Hamiltonian matrix in dictionary
            d['hamiltonian'][lamb] = H_evolved
            
            # Reset initial s value to last s_val and continue the lambda for
            # loop
            s_init = s_val
                
        return d