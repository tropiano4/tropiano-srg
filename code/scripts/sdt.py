#!/usr/bin/env python3

"""
File: sdt.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: January 28, 2020

Functions for comparing two Hamiltonians or potentials based off spectral
distribution theory (SDT). See https://arxiv.org/pdf/1308.5963.pdf and 
https://arxiv.org/pdf/nucl-th/0703076.pdf for more details.

Last update: March 17, 2022

"""

# Python imports
import numpy as np

# Imports from A.T. codes
from .integration import get_factor_meshgrid


def expectation_value(H):
    """Calculates the expectation value of H with dimensionality N."""

    return np.trace(H) / len(H)


def sigma(H, k_array, k_weights, coupled_channel=False):
    """
    Calculates the positive square root of the variance of H.
    
    Parameters
    ----------
    H : 2-D ndarray
        Input Hamiltonian [fm].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
    
    Returns
    -------
    output : float
        Square root of the variance of H [fm].
    
    """

    # Get the integration measure (weights)
    col, row = get_factor_meshgrid(k_array, k_weights, coupled_channel)

    # Compute squared Hamiltonian with factors
    H2_with_factors = (H * col) @ (H * row)

    # Convert H^2 back to mesh-independent quantity (units fm)
    H2 = H2_with_factors / row / col

    return np.sqrt(expectation_value(H2) - expectation_value(H) ** 2)


def inner_product(H, Hp, k_array, k_weights, coupled_channel=False):
    """
    Calculates the inner product of two Hamiltonians in SDT.
    
    Parameters
    ----------
    H : 2-D ndarray
        First input Hamiltonian [fm].
    Hp : 2-D ndarray
        Second input Hamiltonian [fm].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
    
    Returns
    -------
    output : float
        Inner product of H and Hp [fm^2].
    
    """

    # Get the integration measure (weights)
    col, row = get_factor_meshgrid(k_array, k_weights, coupled_channel)

    # HHp is the matrix product of H^{\dagger} H (we're assuming H is real)
    HHp_with_factors = (H.T * col) @ (Hp * row)

    # Convert HHp back to mesh-independent quantity (units fm)
    HHp = HHp_with_factors / row / col

    return (expectation_value(HHp)
            - expectation_value(H.T) * expectation_value(Hp))


def correlation_coefficient(H, Hp, k_array, k_weights, coupled_channel=False):
    """
    Calculates the correlation coefficient between two Hamiltonians in SDT.
    
    Parameters
    ----------
    H : 2-D ndarray
        First input Hamiltonian [fm].
    Hp : 2-D ndarray
        Second input Hamiltonian [fm].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    coupled_channel : bool, optional
        True if the channel is a coupled-channel.
    
    Returns
    -------
    output: float
        Correlation coefficient of H and Hp [unitless].
    
    """

    numerator = inner_product(H, Hp, k_array, k_weights, coupled_channel)
    denominator = (sigma(H, k_array, k_weights, coupled_channel)
                   * sigma(Hp, k_array, k_weights, coupled_channel))

    return numerator / denominator
