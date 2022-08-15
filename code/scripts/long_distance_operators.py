#!/usr/bin/env python3

"""
File: long_distance_operators.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 3, 2019

Contains several functions for calculation of long-distance NN observables.

Last update: May 2, 2022

"""

# Todo: Some functions here only work for the 3S1-3D1. This code still needs
#  to be further generalized.

# Python imports
import numpy as np
from scipy.interpolate import CubicSpline

# Imports from A.T. codes
from .fourier_transform import hankel_transformation_r2k
from .integration import attach_weights_to_matrix, unattach_weights_from_vector
from .tools import build_coupled_channel_matrix


def r2_operator(k_array, k_weights, r_array, dr, a=None, U_matrix=None):
    """
    r^2 operator in momentum space. To calculate the expectation value <r^2>,
    use this operator with the unitless wave functions. For an evolved
    operator, enter an SRG transformation U_matrix. For presentation, one
    should divide out the integration measure by dividing by
    2/pi * k_i * k_j * Sqrt(w_i*w_j) which gives a mesh-independent result.

    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights: 1-D ndarray
        Momentum weights [fm^-1].
    r_array : 1-D ndarray
        Coordinates array [fm].
    dr : float
        Coordinates step-size (weight) [fm].
    a : float, optional
        Parameter in a regulator exponential function exp^(-r^2/a^2). If a
        value is specified, the function will apply the regulator.
    U_matrix : 2-D ndarray, optional
        SRG transformation matrix [unitless]. If no transformation is provided,
        the function will not evolve the operator.

    Returns
    -------
    r2_momentum_space : 2-D ndarray
        r2 operator [fm^2].
    
    Notes
    -----
    This function is specific to the 3S1-3D1 coupled-channel.
        
    """

    # Set up regulator function
    if a is not None:

        regulator = np.exp(-r_array ** 2 / a ** 2)

    else:

        regulator = 1

    # Initialize r^2 in coordinate space first where r^2 is a diagonal matrix
    r2_coordinate_space = np.diag(r_array ** 2) * regulator

    # Transform operator to momentum space
    s_transformation = hankel_transformation_r2k(0, k_array, r_array, dr)
    d_transformation = hankel_transformation_r2k(2, k_array, r_array, dr)

    # Each variable here corresponds to a sub-block of the coupled channel 
    # matrix
    ss_block = s_transformation @ r2_coordinate_space @ s_transformation.T
    dd_block = d_transformation @ r2_coordinate_space @ d_transformation.T

    # Matrix of zeros (m x m) for coupled-channel operator
    ntot = len(k_array)
    zeros = np.zeros((ntot, ntot))

    # Build coupled channel operator with integration measure built-in
    r2_momentum_space_no_weights = build_coupled_channel_matrix(
        ss_block, zeros, zeros, dd_block)
    r2_momentum_space = attach_weights_to_matrix(
        k_array, k_weights, r2_momentum_space_no_weights, coupled_channel=True)

    # Evolve operator?
    if U_matrix is not None:
        r2_momentum_space = U_matrix @ r2_momentum_space @ U_matrix.T

    return r2_momentum_space


def rms_radius_from_rspace(psi, r2_operator):
    """
    Calculates the RMS radius of deuteron using momentum-space wave functions
    and a Fourier transformed r^2 operator (does not involve derivatives in k!)
    Note, set r_max > 25 fm in defining the r^2 operator for effective IR
    cutoff.
    
    Parameters
    ----------
    psi : 1-D ndarray
        Wave function in momentum space [unitless].
    r2_operator : 2-D ndarray
        r^2 operator in momentum space with integration measure [fm^2].
      
    Returns
    -------
    output : float
        RMS radius [fm].
    
    """

    # Calculate the expectation value <\psi|r^2|\psi>
    r2 = psi.T @ r2_operator @ psi

    return 0.5 * np.sqrt(r2)


def rms_radius_from_kspace(psi, k_array, k_weights):
    """
    Same as above function but calculates fully in momentum space where the 
    r^2 operator has derivative terms with respect to k. The expression to
    evaluate is:
        
    <r> = 1/2 * \sqrt[ \int[ dk * ( (k*du/dk)^2 + (k*dw/dk)^2 + 6*w^2 ) ] ]

    Parameters
    ----------
    psi : 1-D ndarray
        Wave function in momentum space [unitless].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
        
    Returns
    -------
    output : float
        RMS radius [fm].
        
    Notes
    -----
    This function is specific to the 3S1-3D1 coupled-channel.
        
    """

    ntot = len(k_array)

    # Load wave function in momentum space
    u_array = unattach_weights_from_vector(k_array, k_weights, psi[:ntot])
    w_array = unattach_weights_from_vector(k_array, k_weights, psi[ntot:])

    # Interpolate with CubicSpline
    # This is a function of k (float or array)
    u_func = CubicSpline(k_array, u_array)
    w_func = CubicSpline(k_array, w_array)

    # Calculate derivatives of u(k) and w(k) using original k_array
    u_deriv_func = u_func.derivative()  # This is a function
    # This is an array with as many points as k_array and k_weights
    u_deriv_array = u_deriv_func(k_array)

    # Same but for D-wave
    w_deriv_func = w_func.derivative()
    w_deriv_array = w_deriv_func(k_array)

    # RMS radius integrand in momentum-space
    integrand = ((k_array * u_deriv_array) ** 2 + (k_array * w_deriv_array) ** 2
                 + 6 * w_array ** 2)

    # The sum over the integrand (which is weighted by k_weights) gives the 
    # value of <r^2>
    r2 = np.sum(k_weights * integrand)

    # The RMS radius of deuteron is given by 1/2 * sqrt(r^2)
    return 0.5 * np.sqrt(r2)


def quadrupole_moment_from_kspace(psi, k_array, k_weights):
    """
    Calculates the quadrupole moment of deuteron fully in momentum space.
    
    Parameters
    ----------
    psi : 1-D ndarray
        Wave function in momentum space [unitless].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    
    Returns
    -------
    output : float
        Quadrupole moment of deuteron [fm^2].
    
    """

    ntot = len(k_array)

    # Load wave function in momentum space
    u_array = unattach_weights_from_vector(k_array, k_weights, psi[:ntot])
    w_array = unattach_weights_from_vector(k_array, k_weights, psi[ntot:])

    # Interpolate with CubicSpline
    # This is a function of k (float or array)
    u_func = CubicSpline(k_array, u_array)
    w_func = CubicSpline(k_array, w_array)

    # Calculate derivatives of u(k) and w(k) using original k_array
    u_deriv_func = u_func.derivative()  # This is a function
    # This is an array with as many points as k_array and k_weights
    u_deriv_array = u_deriv_func(k_array)

    # Same but for D-wave
    w_deriv_func = w_func.derivative()
    w_deriv_array = w_deriv_func(k_array)

    # Quadrupole moment integrand in momentum-space
    integrand = np.sqrt(8) * (
            (k_array ** 2 * u_deriv_array * w_deriv_array
             + 3 * k_array * w_array * u_deriv_array)
            + (k_array * w_deriv_array) ** 2 + 6 * w_array ** 2
    )

    # The sum over the integrand (which is weighted by k_weights) gives the 
    # value of <Q>
    return -1 / 20 * np.sum(k_weights * integrand)
