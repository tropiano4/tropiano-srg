# Created 05/03/19 by A.T. (tropiano.4@osu.edu)

# Functions useful for generating figures.


import datetime
import numpy as np
from scipy.interpolate import interp2d


def current_date():
    '''This function returns the current month and year as string formatted as 
    "Month_Year". This is used to save figures in the correct folder: 
    /Figures/Month_Year.'''
    
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
    '''Interpolate matrix over given array for contour plots. Returns 
    interpolated matrix with new x array with ntot points up to x_max.'''
    
    # Arguments
    
    # x_array (1-D NumPy array): An array of x points (e.g. momentum array)
    # M (2-D NumPy array): A matrix dependent on points x (e.g. Hamiltonian)
    # x_max (float): Maximum value of x for extent of interpolation
    # ntot (integer): Number of points for interpolation
    
    # This is a tricky way to select all the x points in the given array that 
    # are within the specified range of 0 to x_max
    
    # List of boolean values (True or False)
    bool_list = x_array <= x_max 
    
    # The number of points in the given momentum array less than k_max
    n = len( list( filter(None,bool_list) ) )
    
    # Resize x_array and M to size that you want to plot
    x_array = x_array[:n]
    M = M[:n,:n]
    
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
    '''Converts a kvnn number to a label for plotting purposes (e.g. kvnn = 6
    gives 'AV18').'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    
    # Argonne v18
    if kvnn == 6:
        label = 'AV18'
    # Entem/Machleidt N3LO (500 MeV cutoff)   
    elif kvnn == 10:
        label = 'N3LO'
    # Wendt LO non-local potential
    elif kvnn == 119: # Cutoff 4 fm^-1
        label = r'$\Lambda = 4 \/ fm^{-1}$'
    elif kvnn == 120: # Cutoff 9 fm^-1
        label = r'$\Lambda = 9 \/ fm^{-1}$'
    elif kvnn == 121: # Cutoff 20 fm^-1
        label = r'$\Lambda = 20 \/ fm^{-1}$'

    return label