# Created 05/01/19 by A.T. (tropiano.4@osu.edu)

# SRG code: Evolves Hamiltonian to band-diagonal, decoupled form with parameter 
# s using the relative kinetic energy generator.


import numpy as np
from scipy.integrate import odeint
  

class SRG(object):
    
    
    def __init__(self, H0_matrix, T0_matrix):
        '''Saves the initial Hamiltonian and relative kinetic energy in units
        fm^-2 and dimension of the matrices.'''
        
        # Arguments
        
        # H0_matrix (NumPy array): Hamiltonian matrix in units MeV
        # T0_matrix (NumPy array): Relative kinetic energy matrix in units MeV
        
        # h-bar^2 / M [MeV fm^2]
        hbar_sq_over_M = 41.47
        
        # Save matrices in units fm^-2
        self.H0_matrix = H0_matrix / hbar_sq_over_M
        self.T0_matrix = T0_matrix/ hbar_sq_over_M
        # Save dimension of matrix
        self.N = len(H0_matrix)

    
    def commutator(self, A, B):
        '''Returns commutator of A and B, [A,B] where A and B are square matrix
        NumPy arrays.'''
        
        return A @ B - B @ A
    
    
    def matrix2vector(self, A):
        '''Takes the upper triangle of the matrix A (including the diagonal) 
        and reshapes it into a vector B of dimension N*(N+1)/2.'''
    
        # Dimension of matrix
        N = self.N
        # Dimension of vectorized matrix
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
        '''Takes the vector of a upper triangle matrix and returns the full 
        matrix. Use only for hermitian matrices.'''
        
        # Dimension of matrix (given by solving N*(N+1)/2 = n)
        N = self.N
    
        # Initialize matrix
        A = np.zeros((N,N))
    
        # Build upper half of A with diagonal

        a = 0
        b = N

        for i in range(N):

            A[i,i:] = B[a:b]
            a = b
            b += N-i-1

        # Reflect upper half to lower half to build full matrix
        # [np.transpose(A)-np.diag(np.diag(A))] is the lower half of A 
        # excluding the diagonal
        return A+(np.transpose(A)-np.diag(np.diag(A)))
    
    
    def derivs(self, Hs_vector, s):
        '''Returns RHS of SRG flow equation using the Wegner generator.'''
        
        # Arguments
        
        # Hs_vector (NumPy array): Solution vector (which is a function of s)
        # s (float): SRG flow parameter
        
        # Matrix of the solution vector
        Hs_matrix = self.vector2matrix(Hs_vector)

        # Relative kinetic energy SRG generator, eta = [G,H] where G = T_rel 
        G = self.T0_matrix
        
        # SRG generator [G, H(s)]
        eta = self.commutator(G, Hs_matrix)
            
        # RHS of flow equation in matrix form
        dH_matrix = self.commutator(eta, Hs_matrix)
        
        # Returns vector form of RHS of flow equation
        dH_vector = self.matrix2vector(dH_matrix)
        
        return dH_vector


    def evolve_hamiltonian(self, lambda_array):
        '''Returns evolved Hamiltonian Hs_matrix at several values of lambda 
        for given lambda array.'''
    
        # Arguments
        
        # lambda_array (NumPy array): Array of lambda values to be evolved to

        # Set-up ODE
        
        # Reshape initial hamiltonian to a vector
        H0_vector = self.matrix2vector(self.H0_matrix)
        
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