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
#   12/13/19 --- Finished quadrupole_moment_from_kspace function. 
#   03/16/20 --- Updated functions involving observables and their associated 
#                operators to follow the conventions in the notes
#                "NN operator conventions".
#
# Notes:
#   * Some functions here only work for the 3S1 - 3D1 coupled channel. This 
#     code still needs to be further generalized.
#
#------------------------------------------------------------------------------


import numpy as np
import numpy.linalg as la
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import CubicSpline


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
        Hamiltonian matrix [MeV].
    eps : float, optional
        Energy of the desired state [MeV]. Default is the deuteron state.
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
        take phase_shifts[:, j] where j = 0, 1, or 2.
    
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
    Calculates the RMS radius of deuteron using momentum-space wave functions
    and a Fourier transformed r^2 operator (does not involve derivatives in k!)
    Note, set r_max > 25 fm in defining the r^2 operator for mesh-independent
    result.
    
    Parameters
    ----------
    psi : 1-D ndarray
        Full deuteron wave function in momentum-space (unitless). Should be 
        twice the length of k_array and k_weights since it includes the 3S1 
        and 3D1 components.
    r2_operator : 2-D ndarray
        r^2 operator in momentum-space with momenta and weights factored in 
        [fm^2].
      
    Returns
    -------
    output : float
        RMS radius of deuteron [fm].
    
    """
    
    
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


def quadrupole_moment_from_kspace(psi, k_array, k_weights):
    """
    Calculates the quadrupole moment of deuteron fully in momentum-space.
    
    Parameters
    ----------
    psi : 1-D ndarray
        Full deuteron wave function in momentum-space (unitless). Should be 
        twice the length of k_array and k_weights since it includes the 3S1 
        and 3D1 components.
    k_array : 1-D ndarray
        Momentum array.
    k_weights: 1-D ndarray
        Momentum weights.
    
    Returns
    -------
    output : float
        Quadrupole moment of deuteron in units fm^2.
    
    """
        
    # Split psi into 3S1 and 3D1 components
    n = int( len(psi) / 2 )
    u_unitless = psi[:n]
    w_unitless = psi[n:]
    
    # Divide out momenta and weights
    u = u_unitless / ( k_array * np.sqrt(k_weights) ) # fm^3/2
    w = w_unitless / ( k_array * np.sqrt(k_weights) ) # fm^3/2
    
    # Interpolate with CubicSpline (these are functions of k)
    u_func = CubicSpline(k_array, u)
    w_func = CubicSpline(k_array, w)
        
    # Calculate derivatives of u(k) and w(k) using original k_array
    # This is a function
    u_deriv_func = u_func.derivative()
    # This is an array with as many points as k_array and k_weights
    u_deriv = u_deriv_func(k_array)
    w_deriv_func = w_func.derivative()
    w_deriv = w_deriv_func(k_array)
        
    # Quadrupole moment integrand in momentum-space
    integrand = np.sqrt(8) * ( k_array**2 * u_deriv * w_deriv + \
                3 * k_array * w * u_deriv ) + ( k_array * w_deriv )**2 + \
                6 * w**2
        
    # The sum over the integrand (which is weighted by k_weights) gives the 
    # value of <Q>
    return -1/20 * np.sum( k_weights * integrand )