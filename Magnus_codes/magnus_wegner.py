# Created 05/06/19 by A.T. (tropiano.4@osu.edu)

# Magnus code: Evolves Hamiltonian to band-diagonal, decoupled form with 
# parameter s using the Wegner generator.


from math import factorial
import numpy as np
from sympy import bernoulli


class Magnus(object):
    
    
    def __init__(self, H0_matrix, k_magnus):
        '''Saves the initial Hamiltonian in units fm^-2 and dimension of the 
        matrix. Also saves the number of terms to be summed in the Magnus
        Omega(s) derivative sum.'''
        
        # Arguments
        
        # H0_matrix (2-D NumPy array): Hamiltonian matrix in units MeV
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_M = 41.47
        
        # Save matrices in units fm^-2
        self.H0_matrix = H0_matrix / hbar_sq_over_M
        # Save dimension of matrix
        self.N = len(H0_matrix)
        
        # Save truncation in Magnus sum
        self.k_magnus = k_magnus
        
        # Initialize factorial and Bernoulli number arrays for summations up to
        # 30 numbers
        self.k_factorials = np.array([factorial(i) for i in range(31)])
        self.magnus_factors = np.array([bernoulli(j)/factorial(j) for j in range(31)])
        
        
    def commutator(self, A, B):
        '''Returns commutator of A and B, [A,B] where A and B are square matrix
        NumPy arrays.'''
        
        return A @ B - B @ A
    

    def bch_formula(self, M0_matrix, Os_matrix, k_truncate):
        '''Returns an evolved operator given an initial operator M0_matrix and 
        Omega(s).'''
        
        # Arguments
        
        # M0_matrix (2-D NumPy array): Initial operator (which is an N x N matrix)
        # Os_matrix (2-D NumPy array): Evolved Omega matrix
        # k_truncate (integer): Number of terms to include in the sum
        
        # k = 0 ad_Omega^k(H) term
        ad = M0_matrix
        
        # Load factorials
        factorial_array = self.factorial_array
        
        # Zeroth term
        Ms_matrix = ad/factorial_array[0]
        
        # Sum from k = 0 to k = k_truncate
        for k in range(1, k_truncate):
            
            ad = self.commutator(Os_matrix, ad)
            Ms_matrix += ad/factorial_array[k]
            
        return Ms_matrix
    
    
    def derivs(self, Os_matrix, H0_matrix):
        '''Returns the RHS of the Magnus Omega(s) equation where Os_matrix is
        the solution matrix.'''
        
        # Arguments
        
        # Os_matrix (2-D NumPy array): Evolved Omega matrix
        # H0_matrix (2-D NumPy array): Initial Hamiltonian
        
        # Load initial Hamiltonian
        H0_matrix = self.H0_matrix
        
        # Obtain evolved Hamiltonian with BCH formula summing through 25 terms
        Hs_matrix = self.bch_formula(H0_matrix, Os_matrix, 25)
        
        # Wegner SRG generator, eta = [G,H] where G = H_D
        G = np.diag( np.diag(Hs_matrix) )
        
        # SRG generator [G, H(s)]
        eta = self.commutator(G, Hs_matrix)  
        
        # Load Magnus factors in sum
        magnus_factors = self.magnus_factors
        
        # Initial nested commutator ad(eta)
        # k = 0 ad_Omega^k(eta)
        ad = eta
    
        dO_matrix = ad*magnus_factors[0] 
    
        # Sum through 0 to k_magnus
        for k in range(1,self.k_magnus+1):
        
            ad = self.commutator(Os_matrix, ad)
            dO_matrix += ad*magnus_factors[k]
    
        return dO_matrix
    
    
    def euler_method(self, O0_matrix, ds, s_max):
        '''Use first-order Euler method with fixed step-size ds to solve Magnus
        Omega(s) equation.'''
        
        # Arguments
        
        # O0_matrix (2-D NumPy array): Initial Omega matrix
        # ds (float): Step-size in the flow parameter s
        # s_max (float): Maximum value of s
        
        # Initialize Omega(s) matrix
        Os_matrix = O0_matrix
        
        # Step through s values until fully evolved
        s = 0.0
        
        while s <= s_max:
            
            Os_matrix += self.derivs(Os_matrix)*ds
            s += ds
            
        # To ensure Omega(s) stops at s = s_max step size is fixed to go from 
        # last value of s < s_max to s_max where final_ds = s_max - (s-ds)
        Os_matrix += self.derivs(Os_matrix)*(s_max-s+ds)

        return Os_matrix
        
        
    def evolve_hamiltonian(self, lambda_array):
        '''Returns evolved Hamiltonian Hs_matrix at several values of lambda 
        for given lambda array.'''
        
        # Convert 
        
        return None