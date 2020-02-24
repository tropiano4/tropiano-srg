#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: operators.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 29, 2019
# 
# Contains several functions for momentum-space operators. These functions can
# be used in conjunction with functions from observables.py to calculate 
# observable quantities such as the RMS radius. Or the functions can be used to
# look at SRG evolved operators by applying unitary transformations to the
# operator matrices.
#
# Revision history:
#   08/07/19 --- Minor revisions to r^2 operator.
#   08/22/19 --- Changed find_q_index function to use np.fabs() and .argmin()
#   09/23/19 --- Generalized momentum_projection_operator to any channel, not
#                just coupled-channels.
#
# Notes:
#   * The operators here only work for the 3S1 - 3D1 coupled channel. This code
#     still needs to be generalized.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp


def find_q_index(q, k_array):
    """
    Finds the index of the k value nearest to q in the given momentum array.
    For instance, say q = 3.0 fm^-1 and k_array does not contain 3.0 exactly.
    Then this function would return an index corresponding to the k value
    nearest to 3.0 in k_array (e.g. k_array[q_index] = 2.98 fm^-1).
    
    Parameters
    ----------
    q : float
        Momentum value in units fm^-1.
    k_array : 1-D ndarray
        Momentum array.
        
    Returns
    -------
    q_index : int
        Index of q (or nearest k) in k_array.
        
    """
    
    k_difference_array = np.fabs(k_array - q)
    q_index = k_difference_array.argmin()
    
    return q_index


def momentum_projection_operator(q, k_array, k_weights, channel,
                                 U=np.empty(0)):
    """
    ( a_q^dagger a_q ) momentum projection operator in momentum-space. When
    applied to a wave function, returns the wave function at momentum value q.
    For an evolved operator, enter in a unitary transformation U. The initial
    operator is zero everywhere except where k, k' = q. For presentation, one
    should divide out the momenta and weights by dividing by k_i * k_j *
    Sqrt( w_i * w_j). This gives a mesh independent result.

    Parameters
    ----------
    q : float
        Momentum value in units fm^-1.
    k_array : 1-D ndarray
        Momentum array.
    k_weights: 1-D ndarray
        Momentum weights.
    channel : str
        The partial wave channel ('1S0', '3S1', etc.) This allows the function
        to distinguish whether the operator should work for coupled-channels or
        not.
    U : 2-D ndarray, optional
        Unitary transformation matrix. If no unitary transformation is
        provided, the function will skip the line where it evolves the
        operator.
        
    Returns
    -------
    coupled_channel_operator : 2-D ndarray
        Momentum projection operator in units fm^3.
        
    """
        
    # Length of k_array
    m = len(k_array)
        
    # Find index of q in k_array
    q_index = find_q_index(q, k_array)
        
    # Weight for q value
    q_weight = k_weights[q_index]
        
    # Build momentum projection operator 
    operator = np.zeros( (m, m) )
    operator[q_index, q_index] = 1 / ( q**2 * q_weight )
    #operator[q_index, q_index] = np.pi / ( 2 * q**2 * q_weight )
    #operator[q_index, q_index] = 1
    #operator[q_index, q_index] = np.pi / ( 2* q**4 * q_weight**2 )
    
    # Build coupled channel operator 
    if lp.coupled_channel(channel):
    
        # Matrix of zeros (m x m) for coupled-channel operator
        o = np.zeros( (m, m) )
    
        # Build coupled channel operator
        operator = np.vstack( ( np.hstack( (operator, o) ),
                                np.hstack( (o, operator) ) ) )
            
    # Evolve operator by applying unitary transformation U
    if U.any():
        operator = U @ operator @ U.T

    return operator


def hankel_transformation(channel, k_array, r_array):
    """
    <r|k;channel> matrix for given partial wave channel. If len(k_array) = m
    and len(r_array) = n, then this function returns an n x m matrix in units
    UNITS.
    
    Parameters
    ----------
    channel : str
        The partial wave channel ('3S1' or '3D1').
    k_array : 1-D ndarray
        Momentum array.
    r_array : 1-D ndarray
        Coordinates array.
        
    Returns
    -------
    M : 2-D ndarray
        Hankel transformation matrix.

    """
    
    # Grids of k (col), and r (row) values   
    k_cols, r_rows = np.meshgrid(k_array, r_array)
        
    # L = 0 (0th spherical Bessel function)
    if channel[1] == 'S':
        
        L = 0
        
    elif channel[1] == 'P':
        
        L = 1
            
    # L = 2 (2nd spherical Bessel function)
    elif channel[1] == 'D':
        
        L = 2
        
    #M = np.sqrt(2/np.pi) * k_cols**2 * r_rows * spherical_jn(L, k_cols*r_rows)
    M = np.sqrt(2/np.pi) * r_rows * spherical_jn(L, k_cols*r_rows)
    #M = np.sqrt(2/np.pi) * k_cols**2 * spherical_jn(L, k_cols*r_rows)

    return M


def r2_operator(k_array, k_weights, r_array, dr, U=np.empty(0)):
    """
    r^2 operator in momentum-space. For an evolved operator, enter in a unitary 
    transformation U. For presentation, one should divide out the momenta and 
    weights by dividing by k_i * k_j * Sqrt( w_i * w_j). This gives a mesh 
    independent result.

    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array.
    k_weights: 1-D ndarray
        Momentum weights.
    r_array : 1-D ndarray
        Coordinates array.
    dr : float
        Coordinates step-size (weight).
    U : 2-D ndarray, optional
        Unitary transformation matrix. If no unitary transformation is
        provided, the function will skip the line where it evolves the
        operator.
        
    Returns
    -------
    r2_momentum_space : 2-D ndarray
        r2 operator in units XXXXX.
        
    """
        
    # Initialize r^2 in coordinate-space first where r^2 is a diagonal matrix
    r2_coordinate_space = np.diag(r_array**2)
        
    # Transform operator to momentum-space
    s_wave_trans = hankel_transformation('3S1', k_array, r_array)
    d_wave_trans = hankel_transformation('3D1', k_array, r_array)
    
    # Each variable here corresponds to a sub-block of the coupled channel 
    # matrix
    ss_block = s_wave_trans.T @ r2_coordinate_space @ s_wave_trans * dr
    dd_block = d_wave_trans.T @ r2_coordinate_space @ d_wave_trans * dr
        
    # Length of k_array
    m = len(k_array)
        
    # Matrix of zeros (m x m) for coupled-channel operator
    o = np.zeros( (m, m) )
        
    # Build coupled channel operator
    r2_momentum_space = np.vstack( ( np.hstack( (ss_block, o) ),
                                     np.hstack( (o, dd_block) ) ) )
    
    # Evolve operator by applying unitary transformation U
    if U.any():
        r2_momentum_space = U @ r2_momentum_space @ U.T
        
    # Factor of dr for one integration over dr (the other dr' integration is 
    # killed by delta function) - not sure what the weights should be???
    #return r2_momentum_space * dr
    return r2_momentum_space
    #return r2_momentum_space / dr