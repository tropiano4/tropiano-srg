#!/usr/bin/env python3

"""
File: momentum_projection_operator.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 29, 2019

Operator for projecting onto the relative momentum state |q>. This is
sometimes referred to as the momentum occuptation operator. In second
quantization, it is given by a^{\dagger}_q a_q.

Last update: July 8, 2022

"""

# Python imports
import numpy as np

# Imports from A.T. codes
from .integration import unattach_weights_from_vector
from .tools import build_coupled_channel_matrix, find_index


def momentum_projection_operator(
        q, k_array, k_weights, coupled=False, U_matrix=None, smeared=True):
    """
    Momentum projection operator in momentum space. When applied to a wave
    function, returns the wave function at momentum value q. For an evolved
    operator, enter an SRG transformation U_matrix. For presentation, one
    should divide out the integration measure by dividing by
    2/pi * k_i * k_j * Sqrt(w_i*w_j) which gives a mesh-independent result.

    Parameters
    ----------
    q : float
        Momentum value [fm^-1].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights: 1-D ndarray
        Momentum weights [fm^-1].
    coupled : bool, optional
        Whether the operator should work for coupled-channels or not.
    U_matrix : 2-D ndarray, optional
        SRG transformation matrix [unitless]. If no transformation is provided,
        the function will not evolve the operator.
    smeared : bool, optional
        Option on whether the discretized delta function is a smeared delta
        function. Default is smeared version.
        
    Returns
    -------
    operator : 2-D ndarray
        Momentum projection operator [fm^3].
        
    """
    
    ntot = len(k_array)
        
    # Find index of nearest q in k_array
    q_index = find_index(q, k_array)
        
    # Construct delta function
    delta_function_array = np.zeros(ntot)

    # Smeared version given by:
    #   \delta(k-q) = \delta_{k_i, q}/2
    #               + \delta_{k_(i-1), q}/4
    #               + \delta_{k_(i+1), q}/4
    if smeared:
    
        delta_function_array[q_index] = 1/2
        delta_function_array[q_index-1] = 1/4
        delta_function_array[q_index+1] = 1/4

    # Assume \delta(k-q) = \delta_{k_i, q}
    else:
        
        delta_function_array[q_index] = 1
    
    # Divide out integration measure
    delta_function_array = unattach_weights_from_vector(k_array, k_weights,
                                                        delta_function_array)
    
    # Build momentum projection operator
    row, col = np.meshgrid(delta_function_array, delta_function_array,
                           indexing='ij')
    operator = row*col

    # Make coupled-channel operator?
    if coupled:
    
        # Matrix of zeros (ntot x ntot) for off-diagonal blocks
        zeros = np.zeros((ntot, ntot))
    
        # Build coupled channel operator
        operator = build_coupled_channel_matrix(operator, zeros, zeros,
                                                operator)
            
    # Evolve operator?
    if U_matrix is not None:
        operator = U_matrix @ operator @ U_matrix.T

    return operator