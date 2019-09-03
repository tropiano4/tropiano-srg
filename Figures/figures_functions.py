#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: figures_functions.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 3, 2019
# 
# Functions useful for generating figures.
#
# Revision history:
#   09/03/19 --- Updated kvnn_label_conversion function to include more
#                potentials (e.g. RKE, Gezrelis potentials).
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.interpolate import interp2d


def interpolate_matrix(x_array, M, x_max, ntot=500):
    """
    Interpolate matrix over given array for contour plots. Also adds more
    points to given x_array.
    
    Parameters
    ----------
    x_array : 1-D ndarray
        Array of x points where the matrix M = M(x, x).
    M : 2-D ndarray
        Matrix to be interpolated.
    x_max : float
        Maximum value of x. This functions interpolates M up to this value.
    ntot : int, optional
        Desired length of the interpolated matrix and new x array.
        
    Returns
    -------
    x_array_new : 1-D ndarray
        Array of ntot x points.
    M_new : 2-D ndarray
        Interpolated matrix.
        
    """
    
    # This is a tricky way to select all the x points in the given array that 
    # are within the specified range of 0 to x_max
    
    # List of boolean values (True or False)
    bool_list = x_array <= x_max 
    
    # The number of points in the given momentum array less than k_max
    n = len( list( filter(None, bool_list) ) )
    
    # Resize x_array and M to size that you want to plot
    x_array = x_array[:n]
    M = M[:n, :n]
    
    # Use interp2d to interpolate the truncated matrix
    M_func = interp2d(x_array, x_array, M)

    # New x array for interpolation (dimension ntot x 1)
    x_array_new = np.linspace(0.0, x_max, ntot)
    
    # New matrix (dimension ntot x ntot)
    M_new = M_func(x_array_new, x_array_new) 
    
    # We comment out the RectBivariateSpline method in what follows
    #M_func = RectBivariateSpline(x_array, x_array, M)
    #col_mesh, row_mesh = np.meshgrid(x_array_new, x_array_new)      
    #M_new = M_func.ev(row_mesh,col_mesh) # Spline
    
    return x_array_new, M_new


def kvnn_label_conversion(kvnn):
    """
    Converts a kvnn number to a label for plotting purposes (e.g. kvnn = 6 
    gives 'AV18').
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
        
    Returns
    -------
    label : str
        Label for the potential.
        
    """

    # Argonne v18
    if kvnn == 6:
        label = 'AV18'
    # Entem/Machleidt N3LO (500 MeV cutoff)   
    elif kvnn == 10:
        label = 'EM N' + r'$^3$' + 'LO'
    # RKE N3LO (400, 450, 500 MeV cutoffs)
    elif kvnn in [105, 106, 107]:
        label = 'RKE N' + r'$^3$' + 'LO'
    # RKE N4LO (400, 450, 500 MeV cutoffs)
    elif kvnn in [110, 111, 112]:
        label = 'RKE N' + r'$^4$' + 'LO'
    # Gezrelis N2LO (1 and 1.2 fm cutoff)
    elif kvnn in [222, 224]:
        label = 'Gezerlis N' + r'$^2$' + 'LO'
    # Wendt LO non-local potential
    elif kvnn == 900: # Cutoff 4 fm^-1
        label = r'$\Lambda = 4$' + ' fm' + r'$^{-1}$'
    elif kvnn == 901: # Cutoff 9 fm^-1
        label = r'$\Lambda = 9$' + ' fm' + r'$^{-1}$'
    elif kvnn == 902: # Cutoff 20 fm^-1
        label = r'$\Lambda = 20$' + ' fm' + r'$^{-1}$'

    return label


def replace_periods_with_commas(file_name):
    """
    Replaces all periods in a file name with commas. This is necessary for
    adding figures to LaTeX files which doesn't like periods unless they
    specify a file type. For this reason, do not include the file type 
    extension in the file name (i.e. .jpg, .pdf, etc.)
    
    Parameters
    ----------
    file_name : str
        Original name of the file including periods.
    
    Returns
    -------
    new_file_name : str
        New file name replacing periods with commas.
        
    """
    
    # Initialize new file name
    new_file_name = ''
    
    # Loop over each character in the original file name and append to the new
    # file name replacing periods with commas
    for letter in file_name:
        
        if letter == '.':
            
            letter = ','
            
        new_file_name += letter
        
    return new_file_name