#!/usr/bin/env python3

"""
File: fourier_transform.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 31, 2021

Functions for Fourier transforming from momentum space to coordinate space or
vice versa.

Last update: March 17, 2022.

"""

# To-do: Try to understand difference between function here and one in
# spectroscopic_factor.ipynb.

# Python imports
import numpy as np
from scipy.special import spherical_jn


def hankel_transformation_k2r(L, k_array, k_weights, r_array):
    """
    Hankel transformation matrix < r | k L M_L > for given partial wave
    channel. If len(r_array) = m and len(k_array) = n, then this function
    returns an m x n matrix.
    
    Parameters
    ----------
    L : int
        Orbital angular momentum.
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    r_array : 1-D ndarray
        Coordinates array [fm].
        
    Returns
    -------
    M : 2-D ndarray
        Hankel transformation matrix [fm^-3].
        
    Notes
    -----
    The L > 0 transformations may require factors of i or -1.
    Check conventions.
    
    """
        
    # Create meshgrid of k and r
    k_cols, r_rows = np.meshgrid(k_array, r_array)
    _, k_weights_cols = np.meshgrid(k_weights, r_array, indexing='ij')
        
    M = 2/np.pi * k_cols**2 * k_weights_cols * spherical_jn(L, k_cols*r_rows)

    return M


def hankel_transformation_r2k(L, k_array, r_array, dr):
    """
    Hankel transformation matrix < k L M_L | r > for given partial wave
    channel. If len(r_array) = m and len(k_array) = n, then this function
    returns an n x m matrix.
    
    Parameters
    ----------
    L : int
        Orbital angular momentum.
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

    # Create meshgrid of k and r
    r_cols, k_rows = np.meshgrid(r_array, k_array)
        
    M = np.sqrt(dr) * r_cols * spherical_jn(L, k_rows*r_cols)

    return M