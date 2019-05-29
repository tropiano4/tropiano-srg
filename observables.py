# Created 05/03/19 by A.T. (tropiano.4@osu.edu)

# Deuteron code: Calculates observables for deuteron given an NN Hamiltonian 
# (can be SRG evolved). The class contains functions that obtain the deuteron 
# wave function in momentum-space, functions that calculate observable 
# quantities, and functions that return operators in momentum-space.

# Note to self: should probably generalize operator code to another script
# independent of this. This code should use functions from the operator script
# like r^2.


import numpy as np
import numpy.linalg as la
from scipy.special import spherical_jn
from scipy.interpolate import CubicSpline


# Generalize this to any NN state (give energy of state and select closest
# value)? Make epsilon = -2.22 MeV default. Or split this into returning the
# deuteron wave function and observables using NN operators from operator.py.
# Rename this observables.py?


class Deuteron(object):
    
    
    def __init__(self, H_matrix, k_array, k_weights, r_array=np.empty(0), dr=0.0):
        '''Saves the initial Hamiltonian in units MeV, dimension of the matrix,
        the momentum array, and momentum weights. Option to enter in a 
        coordinates array and step-size dr for coordinate-space wave functions
        or operators.''' 
        
        # Arguments
        
        # H_matrix (2-D NumPy array): Hamiltonian matrix in units MeV
        # k_array (1-D NumPy array): Momentum array
        # k_weights (1-D NumPy array): Momentum weights
        # r_array (1-D NumPy array): Coordinates array (default is empty array)
        # dr (float): Step-size for spacing between radial coordinates (default
        # is zero)
        
        # Save Hamiltonian for usage in other functions
        self.H_matrix = H_matrix
        # Save dimension of Hamiltonian
        self.N = len(H_matrix)
        
        # Save momentum and weights
        self.k_array = k_array
        self.k_weights = k_weights
        # Dimension of k_array
        self.m = len(k_array)
        
        # Save coordinates and weight if given
        if dr: # This is false if dr = 0
            self.r_array = r_array
            self.dr = dr
            # Dimension of r_array
            self.n = len(r_array)
    
    
    def wave_func(self, U=np.empty(0)):
        '''Diagonalizes the Hamiltonian and returns the deuteron wave function 
        as u and w corresponding to the 3S1 and 3D1 channels. The wave function
        is unitless, that is, the momenta and weights are factored in such that
        \sum_i { u(k_i)^2 + w(k_i)^2 } = 1. For an evolved wave function, enter
        in a unitary transformation U.'''
        
        # Arguments
        
        # U (2-D NumPy array): Option to apply an SRG or Magnus unitary
        # transformation to the wave function
        
        # Load Hamiltonian
        H_matrix = self.H_matrix 

        # Diagonalize Hamiltonian
        eigenvalues,eigenvectors = la.eig(H_matrix)
    
        # Find the deuteron wave function
        bool_array = (eigenvalues<0.0)*(eigenvalues>-3.0)
        # This is the integer index where eigenvalue[index] = eps_d
        deuteron_index = list(bool_array).index(True) 
        # Full wave function (unitless)
        psi = eigenvectors[:,deuteron_index] 

        # Evolve wave function by applying unitary transformation U
        if U.any():
            psi = U @ psi
        
        u = psi[:120] # 3S1 part 
        w = psi[120:] # 3D1 part
        
        # Check normalization in momentum-space
        #normalization = np.sum((u**2+w**2))
        #print('Normalization = %.4f (k-space)'%normalization)
            
        return u, w
    
    
    def momentum_distribution(self, u, w):
        '''Returns an array of the deuteron momentum distribution given the
        deuteron wave function in components u and w.'''
        
        # Arguments
        
        # u (1-D NumPy array): 3S1 part of the deuteron wave function (unitless)
        # w (1-D NumPy array): 3D1 part of the deuteron wave function (unitless)
        
        # The momentum distribution is given by u^2 + w^2
        psi_squared = u**2+w**2 # Unitless
        
        # Load momentum and weights
        k_array = self.k_array
        k_weights = self.k_weights
        
        # Divide out momenta and weights
        psi_squared_units = psi_squared / (k_array**2 * k_weights) # Units fm^3
        
        return psi_squared_units
    
    
    def momentum_proj_operator(self, q, U=np.empty(0)):
        '''Returns the a_q^dagger a_q projection operator in momentum-space.
        For an evolved operator, enter in a unitary transformation U. Unevolved 
        this matrix should be zero everywhere except where k, k' = q. Note, we 
        return this operator with momenta and weights factored in (the operator 
        is unitless). For presentation, one should divide out the momenta and 
        weights by dividing by k_i * k_j * Sqrt( w_i * w_j ).'''
        
        # Arguments
        
        # q (float): Momentum value in fm^-1
        # U (2-D NumPy array): Option to apply an SRG or Magnus unitary
        # transformation to the operator
        
        # Load momentum
        k_array = self.k_array
        k_weights = self.k_weights
        
        # Dimension of k_array
        m = self.m
        # Matrix of zeros (m x m)
        o = np.zeros((m,m))
        
        # q needs to be a value in the momentum array
        if q in k_array:
        
            # Find index of q in k_array
            q_index = list(k_array).index(q)
            # Weight for q value
            q_weight = k_weights[q_index]
            # Build projection operator 
            proj_operator = np.zeros((m,m))
            proj_operator[q_index,q_index] = 1 / ( q**2 * q_weight )
            # Return coupled channel operator
            proj_operator = np.vstack((np.hstack((proj_operator,o)),\
                                       np.hstack((o,proj_operator))))
            
            # Evolve operator by applying unitary transformation U
            if U.any():
                proj_operator = U @ proj_operator @ U.T

            return proj_operator
        
        # q is not in the momentum mesh
        else:
            print('You need to specify a q value in k_array')
            return np.zeros((2*m,2*m))
    
    
    def hankel_transformation(self, channel):
        '''Returns the <r|k;channel> matrix for given partial wave channel. If 
        len(r_array) = n and len(k_array) = m, then this function returns an 
        n x m matrix. Note, one must specify an r_array and step-size dr.'''
        
        # Arguments
        
        # channel (string): The partial wave channel ('1S0', '3S1', etc.)
        
        # Grids of k (col), and r (row) values   
        k_cols, r_rows = np.meshgrid(self.k_array, self.r_array)
        
        # L = 0 (0th spherical Bessel function)
        if channel == '3S1':
            M = np.sqrt(2/np.pi) * k_cols**2 * r_rows * spherical_jn(0, k_cols*r_rows)
            
        # L = 2 (2nd spherical Bessel function)
        elif channel == '3D1':
            M = np.sqrt(2/np.pi) * k_cols**2 * r_rows * spherical_jn(2, k_cols*r_rows)

        return M    
    
    
    def r2_operator(self, U=np.empty(0)):
        '''Returns the r^2 operator in momentum-space. For an evolved operator,
        enter in a unitary transformation U. If len(k_array) = m, 
        then this function returns an m x m matrix.'''
        
        # Arguments:
        
        # U (2-D NumPy array): Option to apply an SRG or Magnus unitary
        # transformation to the operator
        
        # Load r_array
        r_array = self.r_array
        
        # Initialize r^2 in coordinate-space first where r^2 is a diagonal matrix
        r2_coordinate_space = np.diag(r_array**2)
        
        # Matrix of zeros (m x m)
        m = self.m
        o = np.zeros((m,m))
        
        # Transform to momentum-space and build coupled channel matrix
        s_trans = self.hankel_transformation('3S1') # n x m matrices
        d_trans = self.hankel_transformation('3D1')
        # Each variable here corresponds to a sub-block of the coupled channel matrix
        ss = s_trans.T @ r2_coordinate_space @ s_trans
        dd = d_trans.T @ r2_coordinate_space @ d_trans
        
        # Full coupled channel matrix
        r2_momentum_space = np.vstack((np.hstack((ss,o)),np.hstack((o,dd))))
    
        # Evolve operator by applying unitary transformation U
        if U.any():
            r2_momentum_space = U @ r2_momentum_space @ U.T
        
        # Factor of dr for one integration over dr (the other dr' integration
        # is killed by delta function) ???
        return r2_momentum_space*self.dr
        #return r2_momentum_space
    
    
    def r2_integrand(self, u, w, r2_operator):
        '''Returns the integrand of <r^2> given the wave function in components
        u(k) and w(k), and the r^2 operator in momentum-space. If len(k_array)
        = m, then this function returns an m x m matrix.'''
        
        # Arguments
        
        # u (1-D NumPy array): 3S1 part of the deuteron wave function (unitless)
        # w (1-D NumPy array): 3D1 part of the deuteron wave function (unitless)
        # r2_operator (2-D NumPy array): r^2 operator in momentum-space
        
        # Build full wave function 
        psi = np.concatenate((u,w))
        
        psi_row, psi_col = np.meshgrid(psi, psi)
        # The i,j-th component of the integrand is given by
        # psi(k_i)*r2_operator(k_i,k_j)*psi(k_j)
        # Here we return the entire matrix using np.meshgrid instead of two for
        # loops - this is not matrix multiplication!
        r2_integrand = psi_row * r2_operator * psi_col
        
        return r2_integrand
    
    
    def rms_radius_from_rspace(self, u, w, r2_operator):
        '''Returns the RMS half-radius of deuteron by evaluating 
        0.5 * Sqrt( <psi|r^2|psi> ) in momentum-space where the r^2 operator is
        transformed from coordinate- to momentum-space. Experimental value is 
        ~2.14 fm.'''
        
        return None
    
    
    def rms_radius_from_kspace(self, U=np.empty(0)):
        '''Same as the above function but calculates the r^2 operator in 
        momentum-space explicitly, that is, r^2 depends on derivates d/dk. One
        can check any violation of unitarity by entering in a unitary 
        transformation U.'''
        
        return None
    
    
    def quadrupole_moment_from_rspace(self):
        '''Description.'''
        
        return None
    
    
    def quadrupole_moment_from_kspace(self, U=np.empty(0)):
        '''Same as above function but calculates the Q operator in momentum-
        space explicitly, that is, Q depends on derivatives d/dk. One can check
        any violation of unitarity by entering in a unitary transformation 
        U.'''
        
        return None