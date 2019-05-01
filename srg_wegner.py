# Created 05/01/19 by A.T. (tropiano.4@osu.edu)

# SRG code: Evolves Hamiltonian to band-diagonal, decoupled form with parameter 
# s using the Wegner generator.


import numpy as np
from scipy.integrate import odeint
  

class SRG(object):
    
    
    def __init__(self, H0_matrix):
        '''Saves the initial Hamiltonian in units fm^-2 and dimension of the 
        matrix.'''
        
        # Arguments
        
        # H0_matrix (NumPy array): Hamiltonian matrix in units MeV
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_M = 41.47
        
        # Save matrices in units fm^-2
        self.H0_matrix = H0_matrix / hbar_sq_over_M
        # Save dimension of matrix
        self.N = len(H0_matrix)

    
    def commutator(self, A, B):
        '''Returns commutator of A and B, [A,B] where A and B are square matrix
        NumPy arrays.'''
        
        return A @ B - B @ A

    
    def vector2matrix(self, B):
        '''Takes the vector of a upper triangle matrix and returns the full 
        matrix - i.e. the reverse of using np.triu_indices(). Use only for
        hermitian matrices.'''
        
        # Arguments
        
        # B (NumPy array): Vectorized top right piece of a square matrix
        
        # Dimension of matrix (given by solving N*(N+1)/2 = n)
        N = self.N
    
        # Indices of the upper right triangle of A
        i,j = np.triu_indices(N)
    
        # Initialize matrix
        A = np.zeros((N,N))
    
        # Set upper right triangle
        A[i,j] = B
        # Using hermiticity reflect for the lower triangle
        A[j,i] = B
    
        return A
    
    
    def derivs(self, Hs_vector, s):
        '''Returns RHS of SRG flow equation using the Wegner generator.'''
        
        # Arguments
        
        # Hs_vector (NumPy array): Solution vector (which is a function of s)
        # s (float): SRG flow parameter
        
        # Dimension of matrix
        N = self.N
        
        # Matrix of the solution vector
        Hs_matrix = self.vector2matrix(Hs_vector)

        # Wegner SRG generator, eta = [G,H] where G = H_D(s) 
        G = np.diag( np.diag(Hs_matrix) )
        
        # SRG generator [G, H(s)]
        eta = self.commutator(G, Hs_matrix)
            
        # RHS of flow equation in matrix form
        dH_matrix = self.commutator(eta, Hs_matrix)
        
        # Returns vector form of RHS of flow equation using NumPy's 
        # triu_indices function (which returns the indices of the upper 
        # triangle)
        return dH_matrix[np.triu_indices(N)]


    def evolve_hamiltonian(self, lambda_array):
        '''Returns evolved Hamiltonian Hs_matrix at several values of lambda 
        for given lambda array.'''
    
        # Arguments
        
        # lambda_array (NumPy array): Array of lambda values to be evolved to

        # Set-up ODE
        
        # Dimension of matrix
        N = self.N
        # Reshape initial hamiltonian to a vector
        H0_vector = self.H0_matrix[np.triu_indices(N)]
        
        # Evaluate H(s) at the following values of lambda (or s)
        s_array = np.zeros(len(lambda_array)+1)
        s_array[1:] = 1.0/lambda_array**4.0 # This array includes s = 0

        # Solve the flow equations
        # sol returns a vectorized H(s) at several values of s
        sol = odeint(self.derivs, H0_vector, s_array, atol=1e-06, rtol=1e-06, \
                     mxstep=5000000)

        # Return a dictionary of H(s) at the values of s (i.e., d[1.2] returns 
        # H(lambda=1.2)) which is a matrix
        d = {}
        i = 1
        for lamb in lambda_array:
            
            d[lamb] = self.vector2matrix(sol[i])
            i += 1
        
        return d