#!/usr/bin/env python3

"""
File: long_distance_operators.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 3, 2019

Contains several functions for calculation of long-distance NN observables.

Last update: March 17, 2022

"""

# To-do: Some functions here only work for the 3S1-3D1. This code still needs
# to be further generalized.

# Python imports
import numpy as np
from scipy.interpolate import CubicSpline

# Imports from A.T. codes
from .fourier_transform import hankel_transformation_r2k


def r2_operator(k_array, k_weights, r_array, dr, a=500, U_matrix=np.empty(0)):
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
        Parameter in a regulator exponential function exp^(-r^2/a^2).
        If a > 100 fm, we assume the exponential function is 1.
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
    if a < 100:
    
        regulator = np.exp(-r_array**2/a**2)
        
    else:
        
        regulator = 1
        
    # Initialize r^2 in coordinate space first where r^2 is a diagonal matrix
    r2_coordinate_space = np.diag(r_array**2) * regulator
     
    # Transform operator to momentum space
    s_transformation = hankel_transformation_r2k(0, k_array, r_array, dr)
    d_transformation = hankel_transformation_r2k(2, k_array, r_array, dr)

    # Each variable here corresponds to a sub-block of the coupled channel 
    # matrix
    ss_block = s_transformation @ r2_coordinate_space @ s_transformation.T
    dd_block = d_transformation @ r2_coordinate_space @ d_transformation.T
    
    # Get integration measure
    factor_array = np.concatenate(
        (np.sqrt(2*k_weights/np.pi)*k_array,
         np.sqrt(2*k_weights/np.pi)*k_array)
    )
    col, row = np.meshgrid(factor_array, factor_array)
        
    ntot = len(k_array)
        
    # Matrix of zeros (m x m) for coupled-channel operator
    zeros = np.zeros((ntot, ntot))
        
    # Build coupled channel operator with integration measure built-in
    r2_momentum_space = row * col * np.vstack(
        (np.hstack((ss_block, zeros)), np.hstack((zeros, dd_block)))
    )
    
    # Evolve operator?
    if U_matrix.any():
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
    factor_array = np.sqrt(k_weights) * k_array
        
    # Load wave function in momentum space
    u_array = psi[:ntot] / factor_array
    w_array = psi[ntot:] / factor_array
        
    # Interpolate with CubicSpline
    # This is a function of k (float or array)
    u_func = CubicSpline(k_array, u_array) 
    w_func = CubicSpline(k_array, w_array)
        
    # Calculate derivatives of u(k) and w(k) using original k_array
    u_deriv_func = u_func.derivative() # This is a function
    # This is an array with as many points as k_array and k_weights
    u_deriv_array = u_deriv_func(k_array)
    
    # Same but for D-wave
    w_deriv_func = w_func.derivative()
    w_deriv_array = w_deriv_func(k_array)
        
    # RMS radius integrand in momentum-space
    integrand = ((k_array*u_deriv_array)**2 + (k_array*w_deriv_array)**2
                 + 6*w_array**2)
        
    # The sum over the integrand (which is weighted by k_weights) gives the 
    # value of <r^2>
    r2 = np.sum(k_weights*integrand)
        
    # The RMS radius of deuteron is given by 1/2 * sqrt(r^2)
    return 0.5*np.sqrt(r2)


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
    factor_array = np.sqrt(k_weights) * k_array
        
    # Load wave function in momentum space
    u_array = psi[:ntot] / factor_array
    w_array = psi[ntot:] / factor_array
        
    # Interpolate with CubicSpline
    # This is a function of k (float or array)
    u_func = CubicSpline(k_array, u_array) 
    w_func = CubicSpline(k_array, w_array)
        
    # Calculate derivatives of u(k) and w(k) using original k_array
    u_deriv_func = u_func.derivative() # This is a function
    # This is an array with as many points as k_array and k_weights
    u_deriv_array = u_deriv_func(k_array)
    
    # Same but for D-wave
    w_deriv_func = w_func.derivative()
    w_deriv_array = w_deriv_func(k_array)
        
    # Quadrupole moment integrand in momentum-space
    integrand = np.sqrt(8) * ((k_array**2*u_deriv_array*w_deriv_array
                               + 3*k_array*w_array*u_deriv_array)
                              + (k_array*w_deriv_array)**2 + 6*w_array**2)
        
    # The sum over the integrand (which is weighted by k_weights) gives the 
    # value of <Q>
    return -1/20 * np.sum(k_weights*integrand)