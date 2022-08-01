#!/usr/bin/env python3

"""
File: src_scaling_factor.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 17, 2022

Function that computes the SRC scaling factor a_2. This is the quantity
associated with the plateau of inclusive cross section ratios A/d of electron-
scattering from nucleus A to deuteron d targets at 1.45 < x < 1.9.

Last update: July 25, 2022

"""

# Python imports
import numpy as np


def compute_a2(A, n_p_array, n_n_array, n_d_array, q_array, q_weights):
    """
    SRC scaling factor a_2 evaluated by integrating single-nucleon momentum
    distributions over high momentum.
    
        a2 ~ 2/A*[\int dq q^2 (n_p^A(q) + n_n^A(q))] / [\int dq q^2 n_d(q)]
        
    where q ranges from q_min > k_F to some high value.
    
    Parameters
    ----------
    A : int
        Nuclear mass number (A = Z + N).
    n_p_array : 1-D ndarray
        Proton momentum distribution [fm^3].
    n_n_array : 1-D ndarray
        Neutron momentum distribution [fm^3].
    n_d_array : 1-D ndarray
        Deuteron momentum distribution [fm^3].
    q_array : 1-D ndarray
        High momentum values to integrate over [fm^-1].
    q_weights : 1-D ndarray
        Momentum weights [fm^-1].
    
    Returns
    -------
    output : float
        SRC scaling factor [unitless].
        
    """

    # Calculate deuteron denominator
    denominator = np.sum(q_weights * q_array**2 * n_d_array)
    
    # Calculate probability distribution of the numerator
    p_a_array = q_array**2 / A * (n_p_array+n_n_array)
        
    # Calculate numerator
    numerator = np.sum(q_weights*p_a_array)
        
    # Compute a_2
    return numerator / denominator