#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: figures_functions.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 3, 2019
# 
# Functions useful for generating figures.
#
#------------------------------------------------------------------------------


import datetime
import numpy as np
from scipy.interpolate import interp2d


def current_date():
    """
    This function is used to save figures in the correct folder: 
    /Figures/Month_Year.
    
    Returns
    -------
    output : str
        Current month and year as 'Month_Year'.
    
    """
    
    now = datetime.datetime.now()

    # Current month
    month = now.month
    # Convert month from integer to string
    if month == 1:
        month = 'January'
    elif month == 2:
        month = 'February'
    elif month == 3:
        month = 'March'
    elif month == 4:
        month = 'April'
    elif month == 5:
        month = 'May'
    elif month == 6:
        month = 'June'
    elif month == 7:
        month = 'July'
    elif month == 8:
        month = 'August'
    elif month == 9:
        month = 'September'
    elif month == 10:
        month = 'October'
    elif month == 11:
        month = 'November'
    elif month == 12:
        month = 'December'
        
    # Current year
    year = now.year
    
    # Return current date as Month_Year (e.g. February_1994)
    return '%s_%d'%(month,year)


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
        label = 'N3LO'
    # Wendt LO non-local potential
    elif kvnn == 900: # Cutoff 4 fm^-1
        label = r'$\Lambda = 4 \/ fm^{-1}$'
    elif kvnn == 901: # Cutoff 9 fm^-1
        label = r'$\Lambda = 9 \/ fm^{-1}$'
    elif kvnn == 902: # Cutoff 20 fm^-1
        label = r'$\Lambda = 20 \/ fm^{-1}$'

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