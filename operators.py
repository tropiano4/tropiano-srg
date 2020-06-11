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
#   03/16/20 --- Updated operators to follow the conventions in the notes
#                "NN operator conventions".
#   06/02/20 --- Added option to use smeared delta functions in the momentum
#                projection operator.
#   06/05/20 --- Added option to use exponential regulator function for r^2
#                operator.
#
# Notes:
#   * The operators here only work for the 3S1 - 3D1 coupled channel. This code
#     still needs to be generalized.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lsp


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
                                 U=np.empty(0), smeared=True):
    """
    ( a_q^dagger a_q ) momentum projection operator in momentum-space. When
    applied to a wave function, returns the wave function at momentum value q.
    For an evolved operator, enter in a unitary transformation U. The initial
    operator is zero everywhere except where k, k' = q. For presentation, one
    should divide out the momenta and weights by dividing by 2/pi * k_i * k_j 
    * Sqrt( w_i * w_j) which gives a mesh-independent result.

    Parameters
    ----------
    q : float
        Momentum value [fm^-1].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights: 1-D ndarray
        Momentum weights [fm^-1].
    channel : str
        The partial wave channel (e.g. '1S0'). This allows the function to 
        distinguish whether the operator should work for coupled-channels.
    U : 2-D ndarray, optional
        Unitary transformation matrix with momenta/weights factored in, that 
        is, the matrix is unitless. If no unitary transformation is provided, 
        the function will skip the line where it evolves the operator.
    smeared : bool, optional
        Option on whether the discretized delta function is a smeared delta
        function. Default is smeared version.
        
    Returns
    -------
    operator : 2-D ndarray
        Momentum projection operator [fm^3].
        
    """
        
    # Length of k_array
    m = len(k_array)
        
    # Find index of q in k_array
    q_index = find_q_index(q, k_array)
        
    # Construct delta function
    delta_function_array = np.zeros(m)

    # Assume \delta(k-q) = \delta_{k_i, q} / 2 + \delta_{k_(i-1), q} / 4 + 
    # \delta_{k_(i+1), q} / 4
    # Note, the weighting here is arbitrary
    if smeared:
    
        delta_function_array[q_index] = 1/2
        delta_function_array[q_index-1] = 1/4
        delta_function_array[q_index+1] = 1/4

    else: # Assume \delta(k-q) = \delta_{k_i, q}
        
        delta_function_array[q_index, q_index] = 1
    
    # Divide by momenta/weights (see "NN operator conventions" LaTeX file for
    # more details on this step)
    factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
    delta_function_array /= factor_array
    
    # Build momentum projection operator
    row, col = np.meshgrid(delta_function_array, delta_function_array)
    operator = row * col

    # Build coupled channel operator 
    if lsp.coupled_channel(channel):
    
        # Matrix of zeros (m x m) for coupled-channel operator
        o = np.zeros( (m, m) )
    
        # Build coupled channel operator
        operator = np.vstack( ( np.hstack( (operator, o) ),
                                np.hstack( (o, operator) ) ) )
            
    # Evolve operator by applying unitary transformation U
    if U.any():
        operator = U @ operator @ U.T

    return operator


def hankel_transformation(channel, k_array, r_array, dr):
    """
    <klm|r> matrix for given partial wave channel. If len(r_array) = m
    and len(k_array) = n, then this function returns an n x m matrix.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    r_array : 1-D ndarray
        Coordinates array [fm].
    dr : float
        Coordinates step-size (weight) [fm].
        
    Returns
    -------
    M : 2-D ndarray
        Hankel transformation matrix [fm^3\2].

    """
        
    
    # L = 0 (0th spherical Bessel function)
    if channel[1] == 'S':
        L = 0
    # L = 1
    elif channel[1] == 'P':
        L = 1
    # L = 2
    elif channel[1] == 'D':
        L = 2
        
    # r_array column vectors and k_array row vectors where both grids are
    # n x m matrices
    r_cols, k_rows = np.meshgrid(r_array, k_array)
        
    M = np.sqrt(dr) * r_cols * spherical_jn(L, k_rows * r_cols)

    return M


def r2_operator(k_array, k_weights, r_array, dr, U=np.empty(0), a=1000):
    """
    r^2 operator in momentum-space. For an evolved operator, enter in a unitary 
    transformation U. For presentation, one should divide out the momenta and 
    weights by dividing by 2/pi * k_i * k_j * Sqrt( w_i * w_j ). This gives a
    mesh-independent result. To calculate the expectation value <r^2>, use this
    operator with the unitless wave functions.

    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights: 1-D ndarray
        Momentum weights [fm^-1].
    r_array : 1-D ndarray
        Coordinates array [fm].
    dr : float
        Coordinates step-size (weight) [fm].
    U : 2-D ndarray, optional
        Unitary transformation matrix with momenta/weights factored in, that 
        is, the matrix is unitless. If no unitary transformation is provided, 
        the function will skip the line where it evolves the operator.
    a : float, optional
        Parameter in regulator function. Gives an option to regulate the r^2
        operator in momentum space with an exponential function exp^[-r^2/a^2]
        where a = 5 fm. If a > 100 fm, we assume the exponential function is
        1 (and thus do not regulate the operator).
        
    Returns
    -------
    r2_momentum_space : 2-D ndarray
        r2 operator [fm^2].
        
    * Note, update this function for coupled-channels that aren't 3S1-3D1.
        
    """
    
    # Set up regulator
    if a < 100:
    
        regulator = np.exp( -r_array**2 / a**2 )
        
    else:
        
        regulator = 1
        
    # Initialize r^2 in coordinate-space first where r^2 is a diagonal matrix
    r2_coordinate_space = np.diag(r_array**2) * regulator
     
    # Transform operator to momentum-space
    s_wave_trans = hankel_transformation('3S1', k_array, r_array, dr)
    d_wave_trans = hankel_transformation('3D1', k_array, r_array, dr)

    # Each variable here corresponds to a sub-block of the coupled channel 
    # matrix
    ss_block = s_wave_trans @ r2_coordinate_space @ s_wave_trans.T
    dd_block = d_wave_trans @ r2_coordinate_space @ d_wave_trans.T
    
    # Grids of momenta and weights
    factor_array = np.concatenate( (np.sqrt(k_weights) * k_array, 
                                    np.sqrt(k_weights) * k_array) ) * \
                   np.sqrt(2/np.pi)
    row, col = np.meshgrid(factor_array, factor_array)
        
    # Length of k_array
    n = len(k_array)
        
    # Matrix of zeros (m x m) for coupled-channel operator
    o = np.zeros( (n, n) )
        
    # Build coupled channel operator with momenta/weights
    r2_momentum_space = np.vstack( ( np.hstack( (ss_block, o) ),
                                     np.hstack( (o, dd_block) ) ) ) * row * col
    
    # Evolve operator by applying unitary transformation U
    if U.any():
        r2_momentum_space = U @ r2_momentum_space @ U.T
        
    return r2_momentum_space