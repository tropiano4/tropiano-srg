#!/usr/bin/env python3

"""
File: integration.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 26, 2021

Functions for creating integration points and weights under Gaussian
quadrature.

Last update: May 2, 2022

"""

# Python imports
import numpy as np
from numpy.polynomial.chebyshev import chebgauss
from numpy.polynomial.legendre import leggauss


def gaussian_quadrature_mesh(
        xmax, ntot, xmin=0, xmid=None, nmod=None, method='legendre'):
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
    method : str, optional
        Method to compute sample points and weights. Default is 'legendre' for
        Gauss-Legendre quadrature. Specify 'chebyshev' for Gauss-Chebyshev
        quadrature.

    Returns
    -------
    x_array : 1-D ndarray
        Integration points from [xmin, xmax].
    x_weights: 1-D ndarray
        Integration weights.

    """
    
    # Standard (no mid-point)
    if xmid == None:
        
        # Interval [-1, 1]
        if method == 'legendre':
            y_array, y_weights = leggauss(ntot)
        elif method == 'chebyshev':
            y_array, y_weights = chebgauss(ntot)
        
        # Convert from interval [-1, 1] to [a, b]
        x_array = 0.5 * (y_array+1) * (xmax-xmin) + xmin
        x_weights = (xmax-xmin) / 2 * y_weights
        
    # Split into two meshes
    else:
        
        # In case nmod wasn't specified, take nmod = ntot/2 rounding up
        if nmod == None:
            nmod = round(ntot/2, 0)
        
        # Interval [-1, 1]
        if method == 'legendre':
            y_array_1, y_weights_1 = leggauss(nmod)
            y_array_2, y_weights_2 = leggauss(ntot-nmod)
        elif method == 'chebyshev':
            y_array_1, y_weights_1 = chebgauss(nmod)
            y_array_2, y_weights_2 = chebgauss(ntot-nmod)
        
        # Convert from interval [-1, 1] to [a, b] for both sides of the mesh
        x_array_1 = 0.5 * (y_array_1+1) * (xmid-xmin) + xmin
        x_weights_1 = (xmid-xmin) / 2 * y_weights_1
        x_array_2 = 0.5 * (y_array_2+1) * (xmax-xmid) + xmid
        x_weights_2 = (xmax-xmid) / 2 * y_weights_2
    
        # Combine both sides of the mesh
        x_array = np.concatenate((x_array_1, x_array_2))
        x_weights = np.concatenate((x_weights_1, x_weights_2))
    
    return x_array, x_weights


def momentum_mesh(kmax, kmid, ntot):
    """
    Generate a momentum mesh like the ones in data/potentials.
    
    Parameters
    ----------
    kmax : float
        Maximum value in the momentum mesh [fm^-1].
    kmid : float
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int
        Number of momentum points in mesh.
    
    Returns
    -------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    
    """
    
    k_array, k_weights = gaussian_quadrature_mesh(kmax, ntot, xmid=kmid)
    
    return k_array, k_weights


def get_factor_array(k_array, k_weights, coupled_channel=False):
    """
    Gets the square root of the integration measure as an array, that is,
        \sqrt(2/\pi * w_i) * k_i,
    where i specifies each point in the mesh.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
        
    Returns
    -------
    factor_array: 1-D ndarray
        Square root of the integration measure [fm^-3/2].
    
    """
    
    factor_array = np.sqrt(2/np.pi * k_weights) * k_array
    
    if coupled_channel:
        factor_array = np.concatenate((factor_array, factor_array))
    
    return factor_array


def get_factor_meshgrid(k_array, k_weights, coupled_channel):
    """
    Gets the square root of the integration measure as meshgrids, that is,
        \sqrt(2/\pi) * w_i) * k_i and \sqrt(2/\pi) * w_j) * k_j,
    where the indices i and j run over an operator.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
    
    Returns
    -------
    row : 2-D ndarray
        Rows of square root integration measure meshgrid [fm^-3/2].
    col : 2-D ndarray
        Columns of square root integration measure meshgrid [fm^-3/2].
    
    """
    
    factor_array = get_factor_array(k_array, k_weights, coupled_channel)
    
    row, col = np.meshgrid(factor_array, factor_array, indexing='ij')
    
    return row, col


def attach_weights_to_vector(k_array, k_weights, v, coupled_channel=False):
    """
    Attaches weights to a vector.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    v : 1-D ndarray
        Input vector.
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
    
    Returns
    -------
    output : 1-D ndarray
        Output vector with weights attached. Units are [v] * fm^-3/2.

    """
    
    factor_array = get_factor_array(k_array, k_weights, coupled_channel)

    return v * factor_array


def unattach_weights_from_vector(k_array, k_weights, v, coupled_channel=False):
    """
    Unattaches weights from a vector.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    v : 1-D ndarray
        Input vector.
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
    
    Returns
    -------
    output : 1-D ndarray
        Output vector without weights. Units are [v] * fm^3/2.

    """

    factor_array = get_factor_array(k_array, k_weights, coupled_channel)

    return v / factor_array


def attach_weights_to_matrix(k_array, k_weights, M, coupled_channel=False):
    """
    Attaches weights to matrix.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    M : 2-D ndarray
        Input matrix.
    coupled_channel : bool, optional
        True if the operator is coupled-channel.
    
    Returns
    -------
    output : 2-D ndarray
        Output matrix with weights attached. Units are [M] * fm^-3.

    """
    
    row, col = get_factor_meshgrid(k_array, k_weights, coupled_channel)

    return M * row * col


def unattach_weights_from_matrix(k_array, k_weights, M, coupled_channel=False):
    """
    Unattaches weights from a matrix.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    M : 2-D ndarray
        Input matrix.
    coupled_channel : bool, optional
        True if the operator is coupled-channel.
    
    Returns
    -------
    output : 2-D ndarray
        Output matrix without weights. Units are [M] * fm^3.

    """

    row, col = get_factor_meshgrid(k_array, k_weights, coupled_channel)

    return M / row / col