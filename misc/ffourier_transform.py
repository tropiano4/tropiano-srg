#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: fourier_transform.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 31, 2021
# 
# Functions for Fourier transforming from momentum space to coordinate space
# or the other way around. (Note, running from within this subdirectory will
# cause an error since we're importing vnn.py from another subdirectory.)
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from potentials.vsrg_macos import vnn


def hankel_transformation_k2r(channel, k_array, k_weights, r_array):
    """
    <r|klm> matrix for given partial wave channel. If len(r_array) = m and
    len(k_array) = n, then this function returns an m x n matrix.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
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
        
    # Get L value
    L = vnn.channel_L_value(channel)
        
    # r_array row vectors and k_array column vectors where both grids are
    # n x m matrices
    k_cols, r_rows = np.meshgrid(k_array, r_array)
    k_weights_cols, _ = np.meshgrid(k_weights, r_array)
        
    M = 2/np.pi * k_cols**2 * k_weights_cols * spherical_jn(L, k_cols * r_rows)

    return M


def hankel_transformation_r2k(channel, k_array, r_array, dr):
    """
    <klm|r> transformation matrix for given partial wave channel. If
    len(r_array) = m and len(k_array) = n, then this function returns an 
    n x m matrix.
    
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
        
    # Get L value
    L = vnn.channel_L_value(channel)
        
    # r_array column vectors and k_array row vectors where both grids are
    # n x m matrices
    r_cols, k_rows = np.meshgrid(r_array, k_array)
        
    M = np.sqrt(dr) * r_cols * spherical_jn(L, k_rows * r_cols)

    return M