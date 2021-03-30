#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: integration.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 26, 2021
# 
# Functions for creating integration points and weights under Gaussian
# quadrature.
#
# Revision history:
#   xx/xx/xx --- ...
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss


def gaussian_quadrature_mesh(xmax, ntot, xmin=0, xmid=0, nmod=0):
    """
    Creates an integration mesh from [xmin, xmax] under Gaussian quadrature.
    Option to split into two connected meshes by specifying the mid-point and
    number of points in the interval [xmin, xmid].

    Parameters
    ----------
    xmax : float
        Maximum value in the mesh.
    ntot : int
        Number of points in mesh.
    xmin : float, optional
        Minimum value in the mesh.
    xmid : float, optional
        Mid-point value in the mesh.
    nmod : int, optional
        Number of points in the interval [xmin, xmid]. If nmod=0, the
        function will not split up the mesh into two pieces.

    Returns
    -------
    x_array : 1-D ndarray
        Integration points from [xmin, xmax].
    x_weights: 1-D ndarray
        Integration weights.
        
    Notes
    -----
    * For relative momentum k, it suggested to set xmax = 10 fm^-1, xmid = 2
      fm^-1, and ntot = 120.
    * For total momentum K, it suggested to set xmax = 3 fm^-1 and ntot = 15.
    * In doing integration over angles, we can take a fairly low number of
      points (e.g., ntot = 10).

    """
    
    # Standard (no mid-point)
    if nmod == 0 and xmid == 0:
        
        y_array, y_weights = leggauss(ntot) # Interval [-1, 1]
        
        # Convert from interval [-1, 1] to [a, b] (meaning y_array -> x_array)
        x_array = 0.5 * (y_array + 1) * (xmax - xmin) + xmin
        x_weights = (xmax - xmin) / 2 * y_weights
        
    # Connect low- and high-x meshes    
    else:
        
        y_array_1, y_weights_1 = leggauss(nmod) # Interval [-1, 1]
        y_array_2, y_weights_2 = leggauss(ntot - nmod)
        
        # Convert from interval [-1, 1] to [a, b] (meaning y_array -> x_array)
        x_array_1 = 0.5 * (y_array_1 + 1) * (xmid - xmin) + xmin
        x_weights_1 = (xmid - xmin) / 2 * y_weights_1
    
        x_array_2 = 0.5 * (y_array_2 + 1) * (xmax - xmid) + xmid
        x_weights_2 = (xmax - xmid) / 2 * y_weights_2
    
        x_array = np.concatenate( (x_array_1, x_array_2) )
        x_weights = np.concatenate( (x_weights_1, x_weights_2) )
    
    return x_array, x_weights