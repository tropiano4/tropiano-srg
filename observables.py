#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: observables.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 3, 2019
# 
# Contains several functions for calculation of NN observables given a
# Hamiltonian. These functions can be used to look at wave functions or to 
# calculate observable quantities with functions from operators.py
#
# Revision history:
#   05/29/19 --- This file was renamed from deuteron.py to observables.py. The 
#                idea was to generalize this code to observables for any state
#                given an NN potential.
#   06/06/19 --- Added a function that returns the energies of a given
#                Hamiltonian.
#   08/22/19 --- Changed find_eps_index function to use np.fabs() and .argmin()
#   08/29/19 --- Merged phase_shifts.py code to this script.
#   09/06/19 --- Split phase_shifts function into coupled-channel and normal
#                functions, coupled_channel_phase_shifts and phase_shifts.
#
# Notes:
#   * Some functions here only work for the 3S1 - 3D1 coupled channel. This 
#     code still needs to be further generalized.
#
#------------------------------------------------------------------------------


import numpy as np
import numpy.linalg as la
from scipy.interpolate import RectBivariateSpline
from scipy.special import spherical_jn
#from scipy.interpolate import CubicSpline
# Scripts made by A.T.


def find_eps_index(eps, e_array):
    """
    Finds the index of the energy value nearest to eps in the given energy
    array. For instance, say eps = -2.22 MeV and energies does not contain
    -2.22 exactly. Then this function would return an index corresponding to 
    the energy value nearest to -2.22 in e_array (e.g. 
    e_array[eps_index] = -2.25 MeV).
    
    Parameters
    ----------
    eps : float
        Energy value in units MeV.
    e_array : 1-D ndarray
        Energies array.
        
    Returns
    -------
    eps_index : int
        Index of eps (or nearest e) in e_array.
        
    """
    
    e_difference_array = np.fabs(e_array - eps)
    eps_index = e_difference_array.argmin()
    
    return eps_index


def wave_function(H_matrix, eps=-2.22, U=np.empty(0)):
    """
    Diagonalizes the Hamiltonian and returns the wave function of the state
    nearest energy = eps. The wave function is unitless, that is, the momenta 
    and weights are factored in such that \sum_i { psi(k_i)^2 } = 1. For an 
    evolved wave function, enter in a unitary transformation U.
    
    Parameters
    ----------
    H_matrix : 2-D ndarray
        Hamiltonian matrix in units MeV.
    eps : float, optional
        Energy of the desired state in MeV. Default is the deuteron state.
    U : 2-D ndarray, optional
        Unitary transformation matrix. If no unitary transformation is 
        provided, the function will skip the line where it evolves the wave
        function.
    
    Returns
    -------
    psi : 1-D ndarray
        Wave function for the specified state (unitless).
        
    """
    
    # Diagonalize Hamiltonian
    eigenvalues, eigenvectors = la.eig(H_matrix)
    
    # Index of the wave function
    eps_index = find_eps_index(eps, eigenvalues)
    
    # Full wave function (unitless)
    psi = eigenvectors[:, eps_index] 

    # Evolve wave function by applying unitary transformation U
    if U.any():
        psi = U @ psi
        
    # Check normalization in momentum-space
    #normalization = np.sum( psi**2 )
    #print('Normalization = %.4f (k-space)'%normalization)
            
    return psi


def energies(H_matrix, bound_states_only=True):
    """
    Energies of a given Hamiltonian in units MeV. Option to return only bound
    state energies.
    
    Parameters
    ----------
    H_matrix : 2-D ndarray
        Hamiltonian matrix in units MeV.
    bound_states_only : bool, optional
        If true, returns only bound state energies.
        
    Returns
    -------
    output : 1-D ndarray
        Array of energies in units MeV.
        
    """
        
    # Diagonalize Hamiltonian and obtain eigenvalues
    eigenvalues = np.sort( la.eig(H_matrix)[0] )
        
    # Return only the bound state energies
    if bound_states_only:
            
        return eigenvalues[ eigenvalues < 0.0 ]
        
    # Otherwise return all energies
    else:
            
        return eigenvalues


def phase_shifts(e_array, V_matrix, k_array, k_weights):
    """
    Calculates NN phase shifts as a function of lab energy for a given 
    potential. Note, this function will not work for a coupled-channel
    potential. For details on the calculation, see Phase_Shift_Notes.pdf or
    PHY989_Project1.pdf in the Notes folder
    
    Parameters
    ----------
    e_array : 1-D ndarray
        Array of lab energies in MeV.
    V_matrix : 2-D ndarray
        Potential matrix in units fm.
    k_array : 1-D ndarray
        Momentum array.
    k_weights : 1-D ndarray
        Momentum weights.
    
    Returns
    -------
    phase_shifts : 1-D ndarray
        Array of phase shifts in degrees for each lab energy in e_array.
    
    Notes
    -----
    * Difficult to understand what causes shifts in pi. At the moment, manually
      correcting these shifts is a sloppy fix. Is there a more elegant
      solution?
    
    """
    
    # Set-up
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # Length of the energy array
    M = len(e_array)
    
    # Maximum momentum value in fm^-1
    k_max = max(k_array)
        
    # Length of the momentum array
    N = len(k_array)
        
    # Interpolate potential with RectBivariateSpline
    V_func = RectBivariateSpline(k_array, k_array, V_matrix)
    
    # Initialize array for phase shifts
    phase_shifts = np.zeros(M)


    # Loop over each lab energy
    for i in range(M):
        
        # Lab energy
        e = e_array[i]
        
        # Momentum corresponding to center of mass energy E_lab / 2 where the 
        # factor of 41.47 converts from MeV to fm^-1
        k0 = np.sqrt( e / 2.0 / hbar_sq_over_M )
        
        # Build D_vector
        
        # First N elements of D_vector
        D_vector = 2.0/np.pi * ( k_weights * k_array**2 ) / \
                   ( k_array**2 - k0**2 )
        # N+1 element of D_vector
        D_last = -2.0/np.pi * k0**2 *( np.sum( k_weights /
                 ( k_array**2 - k0**2 ) ) + np.log( ( k_max + k0 ) / \
                 ( k_max - k0 ) ) / ( 2.0*k0 ) )
        # Append N+1 element to D_vector
        D_vector = np.append(D_vector, D_last) # Length is now N+1
        
        # Append k0 to k_array
        k_full = np.append(k_array, k0)
        
        # Create meshes for interpolation
        col, row = np.meshgrid(k_full, k_full)
        
        # Append k0 points by using the interpolated potential
        V_matrix = V_func.ev(row, col)
            
        # Build F matrix, N+1 x N+1, unitless where F_ij = delta_ij + D_j V_ij
        F_matrix = np.identity(N+1) + np.tile( D_vector, (N+1, 1) ) * V_matrix

        # Calculate R matrix and define extremes of R_matrix
        R_matrix = la.solve(F_matrix, V_matrix) # Units fm

        phase_shifts[i] = np.arctan( -k0 * R_matrix[N, N] )

    # Return phase shifts in degrees
    return np.degrees(phase_shifts)


def coupled_channel_phase_shifts(e_array, V_matrix, k_array, k_weights,
                                 convention='Stapp'):
    """
    Calculates NN phase shifts as a function of lab energy for a given coupled-
    channel potential. For details on the calculation, see
    Phase_Shift_Notes.pdf or PHY989_Project1.pdf in the Notes folder.
    
    Parameters
    ----------
    e_array : 1-D ndarray
        Array of lab energies in MeV.
    V_matrix : 2-D ndarray
        Potential matrix in units fm.
    k_array : 1-D ndarray
        Momentum array.
    k_weights : 1-D ndarray
        Momentum weights.
    convention : str, optional
        Phase shift calculation convention 'Stapp' or 'Blatt'.
    
    Returns
    -------
    phase_shifts : 2-D ndarray
        Array of phase shifts delta_a, delta_b, and epsilon for each lab energy
        in e_array. For example, phase_shifts[i, 0] returns delta_a at the ith
        lab energy in e_array. For an entire array of one type of phase shift,
        take phase_shifts[:, i] where i = 0, 1, or 2.
    
    Notes
    -----
    * Difficult to understand what causes shifts in pi. At the moment, manually
      correcting these shifts is a sloppy fix. Is there a more elegant
      solution?
    
    """
    
    # Set-up
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # Length of the energy array
    M = len(e_array)
    
    # Maximum momentum value in fm^-1
    k_max = max(k_array)
        
    # Length of the momentum array
    N = len(k_array)
        
    # Interpolate potential with RectBivariateSpline doing each sub-block
    # separately for the coupled-channel potential
    V11_func = RectBivariateSpline(k_array, k_array, V_matrix[:N, :N])
    V12_func = RectBivariateSpline(k_array, k_array, V_matrix[:N, N:2*N])
    V21_func = RectBivariateSpline(k_array, k_array, V_matrix[N:2*N, :N])
    V22_func = RectBivariateSpline(k_array, k_array, V_matrix[N:2*N, N:2*N])
    
    # Initialize array for phase shifts
    phase_shifts = np.zeros( (M, 3) )


    # Loop over each lab energy
    for i in range(M):
        
        # Lab energy
        e = e_array[i]
        
        # Momentum corresponding to center of mass energy E_lab / 2 where the 
        # factor of 41.47 converts from MeV to fm^-1
        k0 = np.sqrt( e / 2.0 / hbar_sq_over_M )
        
        # Build D_vector
        
        # First N elements of D_vector
        D_vector = 2.0/np.pi * ( k_weights * k_array**2 ) / \
                   ( k_array**2 - k0**2 )
        # N+1 element of D_vector
        D_last = -2.0/np.pi * k0**2 *( np.sum( k_weights /
                 ( k_array**2 - k0**2 ) ) + np.log( ( k_max + k0 ) /
                 ( k_max - k0 ) ) / ( 2.0*k0 ) )
        # Append N+1 element to D_vector
        D_vector = np.append(D_vector, D_last) # Length is now N+1
        
        # Append k0 to k_array
        k_full = np.append(k_array, k0)
        
        # Create meshes for interpolation
        col, row = np.meshgrid(k_full, k_full)
        
        # Append k0 points by using the interpolated potential
        V11_matrix = V11_func.ev(row, col)
        V12_matrix = V12_func.ev(row, col)
        V21_matrix = V21_func.ev(row, col)
        V22_matrix = V22_func.ev(row, col)
            
        # Build coupled channel potential with k0 points included
        V_matrix = np.vstack( ( np.hstack( (V11_matrix, V12_matrix) ), 
                                np.hstack( (V21_matrix, V22_matrix) ) ) )
            
        # Build F matrix, 2*(N+1) x 2*(N+1), unitless where 
        # F_ij = delta_ij + D_j V_ij
        F_matrix = np.identity( 2*(N+1) ) + np.tile( D_vector, (2*(N+1), 2) ) \
                   * V_matrix

        # Calculate R matrix and define extremes of R_matrix
        R_matrix = la.solve(F_matrix, V_matrix) # Units fm

        # These are scalars!
        R11 = R_matrix[N ,N]
        R12 = R_matrix[N, 2*N+1]
        # R21 = R12
        R22 = R_matrix[2*N+1, 2*N+1]

        # Coupled-channel variables
        eps = 0.5 * np.arctan( 2.0 * R12 / ( R11 - R22 ) )
        R_eps = ( R11 - R22 ) / np.cos( 2.0*eps )
        delta_a = -np.arctan( 0.5 * k0 * ( R11 + R22 + R_eps ) )
        delta_b = -np.arctan( 0.5 * k0 * ( R11 + R22 - R_eps ) )
            
        # Restrict values on phases (not sure how this works!)
        while delta_a - delta_b <= 0:
            delta_a += np.pi
        while delta_a - delta_b > np.pi/2.0:
            delta_b += np.pi
        
        # Blatt convention
        if convention == 'Blatt':
            
            # Manually fix +/- shifts in pi
            if delta_b > 0.0: 
                delta_b -= np.pi
                
            phases = np.array( (delta_a, delta_b, eps) )
        
        # Stapp convention
        else:
        
            eps_bar = 0.5 * np.arcsin( np.sin( 2.0*eps ) * \
                      np.sin( delta_a - delta_b ) )
            delta_bar_a = 0.5 * ( delta_a + delta_b + np.arcsin(
                          np.tan( 2.0*eps_bar ) / np.tan( 2.0*eps ) ) )
            delta_bar_b = 0.5 * ( delta_a + delta_b - np.arcsin(
                          np.tan( 2.0*eps_bar ) / np.tan( 2.0*eps ) ) )
            
            # Manually fix +/- shifts in pi
            if delta_b > 0.0: 
                delta_bar_b -= np.pi
                eps_bar *= -1.0
                
            while delta_bar_a < -100.0 * np.pi / 180.0:
            #while delta_bar_a < 0.0:
                delta_bar_a += np.pi
            if e > 120.0:
                ang = 80.0 * np.pi / 180.0
            else:
                ang = np.pi
            while delta_bar_a > ang:
                delta_bar_a -= np.pi
                
            phases = np.array( (delta_bar_a, delta_bar_b, eps_bar) )
            
        # Append phases to phase_shifts in degrees
        phase_shifts[i, :] = np.degrees(phases)
            
    # End of the loop and return the phases
    # This is an M x 3 dimensional array where M is the length of e_array
    return phase_shifts


def phase_corrector(phase_array):
    """
    Description.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    
    return None


def rms_radius_from_rspace(psi, r2_operator):
    """
    Description.
    
    Parameters
    ----------
    
    Returns
    -------
    
    """
    
    # psi should be unitless
    r2 = psi.T @ r2_operator @ psi
    
    return 0.5 * np.sqrt(r2)


def rms_radius_from_kspace(psi):
    """
    Description.
    
    Parameters
    ----------
    
    Returns
    -------
    
    """
    
    return None


def quadrupole_moment_from_kspace(psi):
    """
    Description.
    
    Parameters
    ----------
    
    Returns
    -------
    
    """
    
    return None

# Old deuteron observables code
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