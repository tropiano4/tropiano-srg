#!/usr/bin/env python3

"""
File: phase_shifts.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 3, 2019

Compute NN scattering phase shifts.

Last update: March 17, 2022

"""

# To-do: Still having issues with shifts in \pi. Better solution?
# To-do: Try a different implementation of calculating phase shifts?

# Python imports
import numpy as np
import numpy.linalg as la
from scipy.interpolate import RectBivariateSpline

# Imports from A.T. codes
from .tools import (
    build_coupled_channel_matrix, decompose_coupled_channel_matrix
)


def phase_shifts(e_array, V_matrix, k_array, k_weights):
    """
    Calculates NN phase shifts as a function of lab energy for a given 
    potential. Note, this function will not work for coupled-channel.
    
    Parameters
    ----------
    e_array : 1-D ndarray
        Array of lab energies [MeV].
    V_matrix : 2-D ndarray
        Potential matrix [fm].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    
    Returns
    -------
    phase_shifts : 1-D ndarray
        Array of phase shifts [deg].
    
    Notes
    -----
    Difficult to understand what causes shifts in \pi. At the moment, manually
    correcting these shifts is a sloppy fix. Is there a more elegant solution?
    
    """
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_m = 41.47
    
    mtot = len(e_array)
    ntot = len(k_array)
    
    # Maximum momentum value in fm^-1
    k_max = round(max(k_array))

    # Interpolate potential with RectBivariateSpline
    V_func = RectBivariateSpline(k_array, k_array, V_matrix)
    
    # Initialize array for phase shifts
    phase_shifts = np.zeros(mtot)

    # Loop over each lab energy
    for i, ie in enumerate(e_array):

        # Momentum corresponding to center of mass energy E_lab / 2 where the 
        # factor of 41.47 converts from MeV to fm^-1
        k0 = np.sqrt( ie / 2.0 / hbar_sq_over_m )
        
        # Build D_vector
        # First ntot elements of D_vector
        D_vector = 2.0/np.pi * (k_weights*k_array**2) / (k_array**2-k0**2)
        # ntot+1 element of D_vector
        D_last = -2.0/np.pi*k0**2*(np.sum(k_weights / (k_array**2-k0**2)) 
                                   + np.log((k_max+k0)/(k_max-k0))/(2.0*k0))
        # Append ntot+1 element to D_vector
        D_vector = np.append(D_vector, D_last) # Length is now ntot+1
        
        # k0 can be appended to end of k_array regardless of its value
        k_full = np.append(k_array, k0)
        
        # Create meshes for interpolation
        col, row = np.meshgrid(k_full, k_full)
        
        # Append k0 points by using the interpolated potential
        V_matrix = V_func.ev(row, col)
            
        # Build F matrix where F_ij = delta_ij + D_j V_ij
        F_matrix = (np.identity(ntot+1)
                    + np.tile(D_vector, (ntot+1, 1))*V_matrix)

        # Calculate R matrix and define extremes of R_matrix
        R_matrix = la.solve(F_matrix, V_matrix) # Units fm

        phase_shifts[i] = np.arctan(-k0*R_matrix[ntot, ntot])

    # Return phase shifts in degrees
    return np.degrees(phase_shifts)


def coupled_channel_phase_shifts(
        e_array, V_matrix, k_array, k_weights, convention='Stapp'):
    """
    Calculates NN phase shifts as a function of lab energy for a given
    coupled-channel potential. 
    
    Parameters
    ----------
    e_array : 1-D ndarray
        Array of lab energies [MeV].
    V_matrix : 2-D ndarray
        Potential matrix [fm].
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    convention : str, optional
        Phase shift calculation convention 'Stapp' or 'Blatt'.
    
    Returns
    -------
    phase_shifts : 2-D ndarray
        Array of phase shifts delta_a, delta_b, and epsilon (all units [deg]).
        For example, phase_shifts[i, 0] returns delta_a at the ith lab energy.
        For an entire array of one type of phase shift, take 
        phase_shifts[:, j] where j = 0, 1, or 2.
    
    """
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_m = 41.47
    
    mtot = len(e_array)
    ntot = len(k_array)
    
    # Maximum momentum value in fm^-1
    k_max = round(max(k_array))
    
    # Get V sub-blocks
    V11, V12, V21, V22 = decompose_coupled_channel_matrix(V_matrix)
        
    # Interpolate potential with RectBivariateSpline doing each sub-block
    # separately for the coupled-channel potential
    V11_func = RectBivariateSpline(k_array, k_array, V11)
    V12_func = RectBivariateSpline(k_array, k_array, V12)
    V21_func = RectBivariateSpline(k_array, k_array, V21)
    V22_func = RectBivariateSpline(k_array, k_array, V22)
    
    # Initialize array for phase shifts
    phase_shifts = np.zeros((mtot, 3))

    # Loop over each lab energy
    for i, ie in enumerate(e_array):

        # Momentum corresponding to center of mass energy E_lab / 2 where the 
        # factor of 41.47 converts from MeV to fm^-1
        k0 = np.sqrt( ie / 2.0 / hbar_sq_over_m )
        
        # Build D_vector
        # First ntot elements of D_vector
        D_vector = 2.0/np.pi * (k_weights*k_array**2) / (k_array**2-k0**2)
        # ntot+1 element of D_vector
        D_last = -2.0/np.pi*k0**2*(np.sum(k_weights / (k_array**2-k0**2)) 
                                   + np.log((k_max+k0)/(k_max-k0))/(2.0*k0))
        # Append ntot+1 element to D_vector
        D_vector = np.append(D_vector, D_last) # Length is now ntot+1
        
        # k0 can be appended to end of k_array regardless of its value
        k_full = np.append(k_array, k0)
        
        # Create meshes for interpolation
        col, row = np.meshgrid(k_full, k_full)
        
        # Append k0 points by using the interpolated potential
        V11_matrix = V11_func.ev(row, col)
        V12_matrix = V12_func.ev(row, col)
        V21_matrix = V21_func.ev(row, col)
        V22_matrix = V22_func.ev(row, col)
            
        # Build coupled-channel potential with k0 points included
        V_matrix = build_coupled_channel_matrix(V11_matrix, V12_matrix,
                                                V21_matrix, V22_matrix)

        # Build F matrix where F_ij = delta_ij + D_j V_ij
        F_matrix = (np.identity(2*(ntot+1))
                    + np.tile(D_vector, (2*(ntot+1), 2))*V_matrix)

        # Calculate R matrix and define extremes of R_matrix
        R_matrix = la.solve(F_matrix, V_matrix) # Units fm

        # These are scalars!
        R11 = R_matrix[ntot ,ntot]
        R12 = R_matrix[ntot, 2*ntot+1]
        # R21 = R12
        R22 = R_matrix[2*ntot+1, 2*ntot+1]

        # Coupled-channel variables
        eps = 0.5 * np.arctan(2.0*R12/(R11-R22))
        R_eps = (R11-R22) / np.cos(2.0*eps)
        delta_a = -np.arctan(0.5*k0*(R11+R22+R_eps))
        delta_b = -np.arctan(0.5*k0*(R11+R22-R_eps))
            
        # Restrict values on phases
        while delta_a - delta_b <= 0:
            delta_a += np.pi
        while delta_a - delta_b > np.pi/2.0:
            delta_b += np.pi
        
        # Blatt convention
        if convention == 'Blatt':
            
            # Manually fix +/- shifts in \pi
            if delta_b > 0.0: 
                delta_b -= np.pi
                
            phases = np.array((delta_a, delta_b, eps))
        
        # Stapp convention
        else:
        
            eps_bar = 0.5*np.arcsin(np.sin(2.0*eps) * np.sin(delta_a-delta_b))
            delta_bar_a = 0.5 * (
                delta_a + delta_b 
                + np.arcsin(np.tan(2.0*eps_bar)/np.tan(2.0*eps))
            )
        
            delta_bar_b = 0.5 * (
                delta_a + delta_b
                - np.arcsin(np.tan(2.0*eps_bar)/np.tan(2.0*eps))
            )
            
            # Manually fix +/- shifts in \pi
            if delta_b > 0.0: 
                delta_bar_b -= np.pi
                eps_bar *= -1.0
                
            while delta_bar_a < -100.0*np.pi/180.0:
            #while delta_bar_a < 0.0:
                delta_bar_a += np.pi
            if ie > 120.0:
                ang = 80.0*np.pi/180.0
            else:
                ang = np.pi
            while delta_bar_a > ang:
                delta_bar_a -= np.pi
                
            phases = np.array((delta_bar_a, delta_bar_b, eps_bar))
            
        # Append phases to phase_shifts in degrees
        phase_shifts[i, :] = np.degrees(phases)
            
    # End of the loop and return the phases
    # This is an mtot x 3 dimensional array where mtot is the length of e_array
    return phase_shifts