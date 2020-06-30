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
#                potentials (e.g. RKE, Gezerlis potentials). Added a function
#                that returns color arguments for plotting which is useful for
#                putting several curves on the same figure. Also, added a line
#                styles function that works in the same way as xkcd_colors.
#   09/05/19 --- Added generator_label_conversion function. Analogous to
#                kvnn_label_conversion function.
#   09/09/19 --- Added channel_label_conversion function. Analogous to
#                kvnn_label_conversion and generator_label_conversion
#                functions.
#   09/17/19 --- Changed kvnn label for kvnn = 900-902 from Lambda = ... fm^-1
#                to just ... fm^-1.
#   02/26/20 --- Added lambda_label_conversion function to label \lambda or
#                block-diagonal \Lambda_BD.
#   06/03/20 --- Made vector interpolation function analogous to
#                interpolate_matrix function.
#   06/25/20 --- Added a couple more linestyles to line_styles function and
#                updated replace_periods_with_commas to replace_periods.
#   06/30/20 --- Added convert_ticks_to_labels function.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.interpolate import interp1d, interp2d


def channel_label_conversion(channel):
    """
    Converts a channel string argument to a label for plotting purposes.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
        
    Returns
    -------
    label : str
        Label for the channel.
        
    """

    # Make numbers super- and sub-scripts
    label = r'$^%s$' % channel[0] + channel[1] + r'$_%s$' % channel[2]
    
    return label
        

def generator_label_conversion(generator, lambda_bd=0.00):
    """
    Converts a generator string argument to a label for plotting purposes (e.g.
    generator = 'Wegner' gives r'$G = H_D$').
    
    Parameters
    ----------
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling (e.g. 2.00 fm^-1).
        
    Returns
    -------
    label : str
        Label for the generator.
        
    """

    # Block-diagonal
    if generator == 'Block-diag':
        
        # In the case, lambda_bd is not entered as an argument, the function
        # returns only G = H_{BD}
        if lambda_bd == 0.00:
            label = r'$G=H_{BD}$'
            
        # If lambda_bd is already an integer, format as an integer
        elif int(lambda_bd) == lambda_bd:
            label = r'$G=H_{BD}$' + ' (%d fm' % lambda_bd + r'$^{-1}$' + ')'
            
        # Otherwise, present with two decimal places
        else:
            label = r'$G=H_{BD}$' + ' (%.2f fm' % lambda_bd + r'$^{-1}$' + ')'
            
    elif generator == 'Wegner':
        label = r'$G=H_{D}$'
        
    elif generator == 'T':
        label = r'$G=T_{rel}$'

    return label


def interpolate_matrix(x_array, M, x_max, ntot=500):
    """
    Interpolate matrix over given array for contour plots.
    
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
    
    # List of boolean values
    bool_list = x_array <= x_max 
    
    # The number of points in the given momentum array less than x_max
    n = len( list( filter(None, bool_list) ) )
    
    # Resize x_array and M to size that you want to plot
    x_array_resized = x_array[:n]
    M_resized = M[:n, :n]
    
    # Use interp2d to interpolate the truncated matrix
    M_func = interp2d(x_array_resized, x_array_resized, M_resized)

    # New x array for interpolation (dimension ntot x 1)
    x_array_new = np.linspace(x_array_resized[0], x_array_resized[-1], ntot)
    
    # New matrix (dimension ntot x ntot)
    M_new = M_func(x_array_new, x_array_new) 
    
    # We comment out the RectBivariateSpline method in what follows
    #M_func = RectBivariateSpline(x_array_resized, x_array_resized, M_resized)
    #col, row = np.meshgrid(x_array_new, x_array_new)      
    #M_new = M_func.ev(row, col) # Spline
    
    return x_array_new, M_new


def interpolate_vector(x_array, y_array, x_max, ntot=500, order='linear'):
    """
    Interpolate vector given array for plots.
    
    Parameters
    ----------
    x_array : 1-D ndarray
        Array of x points where the vector y_array = y_array(x).
    y_array : 1-D ndarray
        Vector to be interpolated.
    x_max : float
        Maximum value of x. This functions interpolates y_array up to this 
        value.
    ntot : int, optional
        Desired length of the interpolated vector and new x array.
    order : str, optional
        Order for interpolation. Default is 'linear'.
        
    Returns
    -------
    x_array_new : 1-D ndarray
        Array of ntot x points.
    y_array_new : 1-D ndarray
        Interpolated vector.
        
    """
    
    # This is a tricky way to select all the x points in the given array that 
    # are within the specified range of 0 to x_max
    
    # List of boolean values (True or False)
    bool_list = x_array <= x_max 
    
    # The number of points in the given momentum array less than k_max
    n = len( list( filter(None, bool_list) ) )
    
    # Resize x_array and M to size that you want to plot
    x_array_resized = x_array[:n]
    y_array_resized = y_array[:n]
    
    # Use interp1d to interpolate the truncated vector
    y_func = interp1d(x_array_resized, y_array_resized, kind=order)

    # New x array for interpolation (dimension ntot x 1)
    x_array_new = np.linspace(x_array_resized[0], x_array_resized[-1], ntot)
    
    # New y array (dimension ntot x 1)
    y_array_new = y_func(x_array_new) 
    
    return x_array_new, y_array_new


def kvnn_label_conversion(kvnn, full_label=True):
    """
    Converts a kvnn number to a label for plotting purposes (e.g. kvnn = 6 
    gives 'AV18').
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    full_label : bool, optional
        For some labels, there is a shortened version. Set full_label = False
        for the shortened version. For example, kvnn = 902 gives
        '\Lambda = 20 fm^-1' normally, but the shortened version is '20 fm^-1'.
        
    Returns
    -------
    label : str
        Label for the potential.
        
    Notes
    -----
    It might be worth adding the cutoff to the full_label version for all the
    potentials (e.g. EM N3LO 500 MeV).
        
    """

    # Argonne v18
    if kvnn == 6:
        label = 'AV18'
        
    # Entem/Machleidt N3LO (500 MeV cutoff)   
    elif kvnn == 10:
        label = 'EM N' + r'$^3$' + 'LO'
        
    # EMN N4LO (450, 500, 550 MeV cutoffs)
    elif kvnn in [74, 79, 84]:
        label = 'EMN N' + r'$^4$' + 'LO'
        
    # RKE N3LO (400, 450, 500 MeV cutoffs)
    elif kvnn in [105, 106, 107]:
        label = 'RKE N' + r'$^3$' + 'LO'
    # RKE N4LO (400, 450, 500, 550 MeV cutoffs)
    elif kvnn in [110, 111, 112, 113]:
        label = 'RKE N' + r'$^4$' + 'LO'
        
    # Gezerlis N2LO (1 and 1.2 fm cutoff)
    elif kvnn in [222, 224]:
        label = 'Gezerlis N' + r'$^2$' + 'LO'
        
    # Wendt LO non-local potential
    elif kvnn == 900: # Cutoff 4 fm^-1
        if full_label:
            label = r'$\Lambda = 4$' + ' fm' + r'$^{-1}$'
        else:
            label = '4 fm' + r'$^{-1}$'
    elif kvnn == 901: # Cutoff 9 fm^-1
        if full_label:
            label = r'$\Lambda = 9$' + ' fm' + r'$^{-1}$'
        else:
            label = '9 fm' + r'$^{-1}$'
    elif kvnn == 902: # Cutoff 20 fm^-1
        if full_label:
            label = r'$\Lambda = 20$' + ' fm' + r'$^{-1}$'
        else:
            label = '20 fm' + r'$^{-1}$'

    return label


def lambda_label_conversion(lamb, generator='Wegner'):
    """
    Converts a lambda evolution parameter to a label for plotting purposes 
    (e.g. lamb = 2 gives r'$\lambda=2.0 fm^-1$').
    
    Parameters
    ----------
    lamb : float
        Evolution parameter lambda [fm^-1].
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'. Determines whether
        lambda is referring to \lambda or \Lambda_BD.
        
    Returns
    -------
    label : str
        Label for the lambda.
        
    """
    
    # Case for lamb = infinity
    if lamb == np.inf:
        
        # Label \lambda_bd
        if generator == 'Block-diag':
            label = r'$\Lambda_{BD}=\infty$' + ' fm' + r'$^{-1}$'
            
        # Label \lambda
        else:
            label = r'$\lambda=\infty$' + ' fm' + r'$^{-1}$'
            
    # Case for finite lamb
    else:
        
        # Label \lambda_bd
        if generator == 'Block-diag':
            label = r'$\Lambda_{BD}=%.1f$' % lamb + ' fm' + r'$^{-1}$'
            
        # Label \lambda
        else:
            label = r'$\lambda=%.1f$' % lamb + ' fm' + r'$^{-1}$'
            
    return label
        

def line_styles(curve_number):
    """
    Default line styles for plotting curves. Styles are set by the ordering
    of curves on the plot, that is, the first curve corresponds to 'solid', 
    the second to 'dashdot', etc. Do not use this function for plots with more
    than six curves. An error message will display when a number higher than
    six is used and the function will return 'solid'.
    
    Parameters
    ----------
    curve_number : int
        The curve number being assigned a line style.
    
    Returns
    -------
    out : str
        The line style argument to be used in plotting functions. For example,
        
        plt.plot(x, y_1, linestyle=line_style(0))
        plt.plot(x, y_2, linestyle=line_style(3))
        
        plots two curves as solid and dotted lines.
    
    """
    
    # Note, the last two are 'densely dashed' and 'densely dashdotted'
    line_styles = ['solid', 'dashdot', 'dashed', 'dotted', (0, (5, 1)), 
                   (0, (3, 1, 1, 1))]
    
    try:
        
        return line_styles[curve_number]
    
    except IndexError:
        
        error = 'Curve number is too high to match with default line style.'
        suggestion = 'Manually assign styles to each curve.'
        print_statement = error + ' ' + suggestion
        
        print(print_statement)
        
        return 'solid'


def replace_periods(file_name):
    """
    Replaces all periods in a file name with a '_'. This is necessary for 
    adding figures to LaTeX files which don't like periods unless they specify
    a file type. For this reason, do not include the file type extension in 
    the file name (i.e. .jpg, .pdf, etc.)
    
    Parameters
    ----------
    file_name : str
        Original name of the file including periods.
    
    Returns
    -------
    new_file_name : str
        New file name replacing periods with '_'.
        
    """
    
    # Initialize new file name
    new_file_name = ''
    
    # Loop over each character in the original file name and append to the new
    # file name replacing periods with '_'
    for letter in file_name:
        
        if letter == '.':
            
            letter = '_'
            
        new_file_name += letter
        
    return new_file_name


def convert_ticks_to_labels(ticks):
    """
    Converts axes or colorbar ticks to string formatted labels displaying the
    correct number of digits.
    
    Parameters
    ----------
    ticks : 1-D ndarray
        Array of tick values for given axis or colorbar.
        
    Returns
    -------
    tick_labels : list
        List of tick labels which are strings.
        
    Notes
    -----
    This function assumes a linearly-spaced array of ticks. Needs further
    generalization for log-spacing.
        
    """
    
    # Find the spacing between ticks assuming linearly-spaced ticks
    stepsize = ticks[1] - ticks[0]
    
    # Figure out many digits to display for the tick labels
    i = 0
    while round(stepsize, i) != stepsize:
        i += 1
    
    # If i = 0, then display integers
    if i == 0:
        tick_labels = ['%d' % tick for tick in ticks]
    # Otherwise, display floats with the correct number of digits
    elif i == 1:
        tick_labels = ['%.1f' % tick for tick in ticks]
    elif i == 2:
        tick_labels = ['%.2f' % tick for tick in ticks]
    elif i == 3:
        tick_labels = ['%.3f' % tick for tick in ticks]
    elif i == 4:
        tick_labels = ['%.4f' % tick for tick in ticks]
    else:
        tick_labels = ['%.f' % tick for tick in ticks]
        
    return tick_labels
    

def xkcd_colors(curve_number):
    """
    Default curve colors for plotting using xkcd colors. Colors are set by the
    ordering of curves on the plot, that is, the first curve corresponds to
    black, the second to red, etc. Do not use this function for plots with
    more than seven curves. An error message will display when a number higher
    than seven is used and the function will return 'xkcd:black'.
    
    Parameters
    ----------
    curve_number : int
        The curve number being assigned a color.
    
    Returns
    -------
    out : str
        The color argument to be used in plotting functions. For example,
        
        plt.plot(x, y_1, color=xkcd_colors(0))
        plt.plot(x, y_2, color=xkcd_colors(1))
        
        plots two curves as black and red lines.
    
    """
    
    colors = ['xkcd:black', 'xkcd:red', 'xkcd:blue', 'xkcd:green',
              'xkcd:orange', 'xkcd:purple', 'xkcd:grey']
    
    try:
        
        return colors[curve_number]
    
    except IndexError:
        
        error = 'Curve number is too high to match with default color.'
        suggestion = 'Manually assign colors to each curve.'
        print_statement = error + ' ' + suggestion
        
        print(print_statement)
        
        return 'xkcd:black'
