#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: phase_shifts.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     July 11, 2018 
# 
# Calculates NN phase shifts as a function of lab energy for a given potential.
# For details on the calculation, see Phase_Shift_Notes.pdf and 
# PHY989_Project1.pdf in the Notes folder.
#
# Revision history:
#   September 25, 2018 --- Updated to manually correct +/- shifts of pi in
#                          phase shifts.
#   May 28, 2019       --- Replaced SciPy's interp2d with RectBivariateSpline
#                          for interpolation.
#
# Notes:
#   * This code needs to be generalized to non-coupled potentials.
#   * Difficult to understand what causes shifts in pi. At the moment, manually
#     correcting these shifts is a sloppy fix. Is there a more elegant
#     solution?
#
#------------------------------------------------------------------------------


import numpy as np
import numpy.linalg as la
from scipy.interpolate import RectBivariateSpline


class Phase_shifts(object):
    
    
    def __init__(self, V_matrix, k_array, k_weights, convention='Stapp'):
        """
        Define constants and interpolate potential.
        
        Parameters
        ----------
        V : 2-D ndarray
            Potential matrix in units fm.
        k_array : 1-D ndarray
            Momentum array.
        k_weights : 1-D ndarray
            Momentum weights.
        coupled_channel : bool, optional
            True if the channel is coupled channel and false otherwise.
        convention : str, optional
            Phase shift calculation convention - 'Stapp' or 'Blatt'.

        """

        # Save momentum and weights (units are fm^-1)
        self.k_array = k_array
        self.k_weights = k_weights
        
        # Maximum momentum value in fm^-1
        self.kmax = max(k_array)
        
        # Save length of momentum
        N = len(k_array)
        self.N = N
        
        # Interpolate potential with RectBivariateSpline
        # Interpolate each sub-block seperately for coupled-channel potential
        self.V11_func = RectBivariateSpline(k_array, k_array, V_matrix[:N, :N])
        self.V12_func = RectBivariateSpline(k_array, k_array, 
                                            V_matrix[:N, N:2*N])
        self.V21_func = RectBivariateSpline(k_array, k_array, 
                                            V_matrix[N:2*N, :N])
        self.V22_func = RectBivariateSpline(k_array, k_array,
                                            V_matrix[N:2*N, N:2*N])
        
        # Save convention
        self.convention = convention

        
    def delta(self, e):
        """
        Phase shift as a function of lab energy, e [MeV].
        
        Parameters
        ----------
        e : float
            Lab energy in units MeV.
            
        Returns
        -------
        phase_shifts : 1-D ndarray
            delta_a, delta_b, and eps phase shifts in degrees if Blatt 
            convention. delta_bar_a, delta_bar_b, and eps_bar phase shifts in
            degrees if Stapp convention.
            
        """
        
        # Length of momentum array
        N = self.N
        
        # Momentum corresponding to center of mass energy E_lab / 2 where
        # the factor of 41.47 converts from MeV to fm^-1
        k0 = np.sqrt( e / 2.0 / 41.47 )
        
        # Load maximum k value
        kmax = self.kmax
        
        # Build u_j vector
        
        # Load momentum and weights
        k_vec = self.k_array
        w_vec = self.k_weights
        
        # First N elements of u_vec
        u_vec = 2.0/np.pi * ( w_vec * k_vec**2 ) / ( k_vec**2 - k0**2 )
        # N+1 element of u_vec
        u_last = -2.0/np.pi * k0**2 *( np.sum( w_vec / ( k_vec**2-k0**2 ) ) + \
                 np.log( ( kmax + k0 ) / ( kmax - k0 ) ) / ( 2.0*k0 ) )
        # Append N+1 element to u_vec
        u_vec = np.append(u_vec, u_last) # Length is now N+1
        
        # Append k0 to k_array
        k_full = np.append(self.k_array, k0)
        
        # Create meshes for interpolation
        col, row = np.meshgrid(k_full, k_full)
        
        # Append k0 points by using the interpolated potential
        v11 = self.V11_func.ev(row, col)
        v12 = self.V12_func.ev(row, col)
        v21 = self.V21_func.ev(row, col)
        v22 = self.V22_func.ev(row, col)
            
        # Build coupled channel potential with k0 points included
        V_matrix = np.vstack( ( np.hstack( (v11, v12) ), 
                                np.hstack( (v21, v22) ) ) )
        
        # Build A matrix, N+1 x N+1, unitless
        A_matrix = np.identity( 2*(N+1) ) + np.tile( u_vec, (2*(N+1), 2) ) * \
                   V_matrix

        # Calculate R matrices and define extremes of R_matrix
        R_matrix = la.solve(A_matrix, V_matrix) # Units fm

        R11 = R_matrix[N ,N]
        R12 = R_matrix[N, 2*N+1]
        # R21 = R12
        R22 = R_matrix[2*N+1, 2*N+1]

        # Coupled-channel variables
        eps = 0.5 * np.arctan( 2.0 * R12 / ( R11 - R22 ) )
        R_eps = ( R11 - R22 ) / np.cos( 2.0*eps )
        delta_a = -np.arctan( 0.5 * k0 * ( R11 + R22 + R_eps ) )
        delta_b = -np.arctan( 0.5 * k0 * ( R11 + R22 - R_eps ) )
            
        # Restrict values on phases
        while delta_a - delta_b <= 0:
            delta_a += np.pi
        while delta_a - delta_b > np.pi/2.0:
            delta_b += np.pi
        
        # Blatt convention
        if self.convention == 'Blatt':
            
            # Manually fix +/- shifts in pi
            if delta_b > 0.0: 
                delta_b -= np.pi
                
            # Return phase shifts in degrees
            phase_shifts = np.array( (delta_a, delta_b, eps) )
        
        # Stapp convention
        else:
        
            eps_bar = 0.5 * np.arcsin( np.sin( 2.0*eps ) * \
                      np.sin( delta_a - delta_b ) )
            delta_bar_a = 0.5 * ( delta_a + delta_b + np.arcsin(
                          np.tan( 2.0*eps_bar ) / np.tan( 2.0*eps ) ) )
            delta_bar_b = 0.5 * ( delta_a + delta_b - np.arcsin(
                          np.tan( 2.0*eps_bar ) / np.tan( 2.0*eps ) ) )
            
            # Manually fix +/- shifts in pi
            if delta_b > 0.0: 
                delta_bar_b -= np.pi
                eps_bar *= -1.0
                
            while delta_bar_a < -100.0 * np.pi / 180.0:
            #while delta_bar_a < 0.0:
                delta_bar_a += np.pi
            if e > 120.0:
                ang = 80.0 * np.pi / 180.0
            else:
                ang = np.pi
            while delta_bar_a > ang:
                delta_bar_a -= np.pi
                
            phase_shifts = np.array( (delta_bar_a, delta_bar_b, eps_bar) )
            
        # Return phase shifts in degrees
        return phase_shifts * 180.0/np.pi
            
    
    def delta_a_array(self,e_array):
        """
        Calculates an array of delta_a phase shifts for a given array of lab
        energies.
        
        Parameters
        ----------
        e_array : 1-D ndarray
            Array of lab energies in MeV.
            
        Returns
        -------
        deltas : 1-D ndarray
            Array of delta_a(E_lab) in degrees.
            
        """
        
        # Length of lab energies array
        n = len(e_array)
    
        # Initialize deltas array
        deltas = np.zeros(n)
        
        # Loop over each lab energy
        for i in range(n):
            
            e_lab = e_array[i]
        
            # Calculate delta_a
            deltas[i] = self.delta(e_lab)[0]
            
            # Manually correct shifts in +/- pi
            if i > 0:
                
                if deltas[i] < ( deltas[i-1] - 100.0 ):
                    deltas[i] += 180.0
                if deltas[i] > ( deltas[i-1] + 100.0 ):
                    deltas[i] -= 180.0

        return deltas  
    
    
# -----------------------------------------------------------------------------
# Try to do everything above in one function that takes V(k, k'), e_array, and
# convention as arguments
        
    
def phase_shifts(e_array, V_matrix, k_array, k_weights, convention='Stapp'):
    """
    Description.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    
    Notes
    -----
    
    """
    
    return None


def phase_corrector(phase_array):
    """
    Description.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    
    return None