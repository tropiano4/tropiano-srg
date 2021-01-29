#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: magnus_split.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     June 3, 2019
# 
# Evolves Hamiltonian to band-diagonal, decoupled form with flow parameter s 
# using the Wegner generator and the Magnus "split thing" implementation. The
# split thing method involves solving the Euler method equation for the 
# difference in omega matrices,
# 
#   delta_omega = Omega(s_{i+1}) - Omega(s_i) = derivative(Omega, s) * ds.
# 
# We take the matrix delta_omega and use the BCH formula to update the evolving
# Hamiltonian at each step in s,
#
#   H(s_{i+1}) = exp^(delta_omega) * H(s_i) * exp^(-delta_omega).
#
# Since Omega(s_{i+1}) - Omega(s_i) is directly proportional to ds, we can tune
# ds to a small enough value to prevent delta_omega from diverging to infinity.
# This method prevents the Magnus evolution from running into NaN or infinity
# errors. However, we cannot save the full omega matrix since we only deal with
# differences in the omega matrix.
#
# Notes:
#   * For kvnn = 902 or 6 where this method is necessary, use ds = 10^-7. This
#     will take a considerable amount of time to evolve but eventually works 
#     (t ~ 8 days on a laptop).
#   * Currently this script is not set-up in run_magnus.py. For split thing
#     runs, must use this script manually.
#   * For very small step-sizes, ds ~ 10^-7, the Wegner SRG generator matches
#     the relative kinetic energy generator. For this reason, we only need
#     this Wegner version of the split thing method.
#
#------------------------------------------------------------------------------


from math import factorial
import numpy as np
from sympy import bernoulli


class Magnus(object):
    
    
    def __init__(self, H_initial, k_magnus, ds):
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
        ds : float
            Step-size in the flow parameter s
        
        """
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_m = 41.47
        
        # Save matrices in scattering units [fm^-2]
        self.H_initial = H_initial / hbar_sq_over_m
        
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
        matrix. Can also evolve by one step in s by using this function as
        follows: H(s_{i+1}) = bch_formula(H(s_i), delta_omega, k_truncate)
        where delta_omega = Omega(s_{i+1}) - Omega(s_i).

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
        
        # k = 0 term
        ad = M_initial
        
        # Load factorials
        factorial_array = self.factorial_array
        
        # Zeroth term
        M_evolved = ad / factorial_array[0]
        
        # Sum from k = 1 to k = k_truncate
        for k in range(1, k_truncate + 1):
            
            ad = self.commutator(O_evolved, ad)
            M_evolved += ad / factorial_array[k]
            
        return M_evolved

    
    def derivative(self, delta_omega, H_evolved):
        """
        Right-hand side of the Magnus derivative omega equation using the 
        Wegner generator and "split thing" approach. In the split thing
        approach, the RHS is Omega(s_{i+1}) - Omega(s_i).
        
        Parameters
        ----------
        delta_omega : 2-D ndarray
            Difference in evolving omega matrices from flow parameter s_{i+1}
            to s_i.
        H_evolved : 2-D ndarray
            Evolving Hamiltonian which is a matrix and function of s. Units are
            fm^-2.
        
        Returns
        -------
        dO_matrix : 2-D ndarray
            Derivative with respect to s of the evolving omega matrix. 

        """

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
            
            ad = self.commutator(delta_omega, ad)
            dO_matrix += magnus_factors[k] * ad
    
        return dO_matrix


    def evolve_hamiltonian(self, lambda_array):
        """
        Magnus evolved Hamiltonian at each value of lambda in lambda_array
        using the first-order Euler method and "split thing" approach.
        
        Parameters
        ----------
        lambda_array : 1-D ndarray
            Lambda evolution values in units fm^-1.
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each lambda value (e.g. d[1.5]).
            
        """

        # Load initial Hamiltonian and length
        H_matrix = self.H_initial
        N = self.N

        # Set-up ODE
        
        # Initial Magnus omega matrix
        O_initial = np.zeros( (N, N) )
    
        # Initialize dictionary
        d = {}
        
        # Initial s value, step-size, and delta omega matrix
        s = 0.0
        ds = self.ds
        delta_omega = O_initial
    
        # Loop over lambda values in lambda_array
        for lamb in lambda_array:
            
            # Convert lamb to s value
            s_val = 1.0 / lamb**4.0
            
            # Solve ODE up to s_val using the Euler method and store in 
            # dictionary
            while s <= s_val:

                # Next step in lambda using the Euler method
                delta_omega = self.derivative(delta_omega, H_matrix) * ds
                
                # Step to the next value of evolving Hamiltonian with the BCH
                # formula
                H_matrix = self.bch_formula(H_matrix, delta_omega, 25)
                
                # Check for NaN's and infinities
                # If true, stop evolving
                if np.isnan(H_matrix).any() or np.isinf(H_matrix).any():
                    
                    print('_'*85)
                    error = 'Infinities or NaNs encountered in Hamiltonian.'
                    print(error)
                    print('s = %.5e' % s)
                    suggestion = 'Try setting ds to a smaller value.'
                    print(suggestion)
                    
                    return d
                
                # Step to next s value
                s += ds
                
            # To ensure omega stops at s = s_val step-size is fixed to go from 
            # last value of s < s_max to s_max where ds_exact = s_max - (s-ds)
            ds_exact = s_val - s + ds
            delta_omega = self.derivative(delta_omega, H_matrix) * ds_exact
                
            # Step to the next value of evolving Hamiltonian with the BCH 
            # formula
            H_matrix = self.bch_formula(H_matrix, delta_omega, 25)
            
            # Store evolved Hamiltonian matrix in dictionary
            d[lamb] = H_matrix
            
            # Reset initial s value to last s_val and continue the lambda for
            # loop
            s = s_val
                
        return d