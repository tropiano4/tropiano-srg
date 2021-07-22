#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: figures_functions.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 3, 2019
# 
# Useful functions for generating and labeling figures.
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
#   03/24/21 --- Added nuclei_label_conversion function.
#   04/27/21 --- Added coupled_channel function and convert_number_to_string,
#                and updated channel_label_conversion function and
#                lambda_label_conversion function.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.interpolate import interp1d, interp2d


def channel_label_conversion(channel, label_coupled_channel=True):
    """
    Converts a channel string argument to a label for plotting purposes.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
    label_coupled_channel : bool, optional
        Option to label full coupled-channel (e.g., '3S1-3D1' instead of
        '3S1' only).
        
    Returns
    -------
    label : str
        Label for the channel.
        
    """
    
    # Get S, L, and J as strings where L is given by 'S', 'P', etc.
    S = channel[0]
    L = channel[1]
    J = channel[2]
    
    # Make label with super- and sub-scripts
    # Label coupled channel?
    if label_coupled_channel and coupled_channel(channel):
        
        # Get L + 2 value
        if L == 'S':
            L2 = 'D'
        elif L == 'P':
            L2 = 'F'
        elif L == 'D':
            L2 = 'G'
        elif L == 'F':
            L2 = 'H'
        elif L == 'G':
            L2 = 'I'
        elif L == 'H':
            L2 = 'K'
        elif L == 'I':
            L2 = 'L'
        elif L == 'K':
            L2 = 'M'
        elif L == 'L':
            L2 = 'N'
        elif L == 'M':
            L2 = 'O'
        elif L == 'N':
            L2 = 'Q'
            
        label = r'$^{%s}{\rm %s}_{%s}\endash^{%s}{\rm %s}_{%s}$' % (S, L, J,
                                                                    S, L2, J)
        
    else:
        
        label = r'$^{%s}{\rm %s}_{%s}$' % (S, L, J)
    
    return label


def convert_number_to_string(number):
    """
    Gives the input number with the correct amount of digits meaning plots will
    show '1.35' not '1.35000000023'.

    Parameters
    ----------
    number : float
        Some input float (e.g., \lambda or momentum).

    Returns
    -------
    output : str
        Input number but rounded to the correct number of digits and converted
        to a string.

    """
    
    # Loop over i until the rounded number matches the input number
    # then return the string of the rounded number
    i = 0
    while round(number, i) != number:
        i += 1
        
    return str( round(number, i) )


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
    
    # Initialize tick_labels list
    tick_labels = []
    
    # Keep track of the maximum number of digits to be displayed i
    
    for tick in ticks:
        
        i = 0
        while abs( round(tick, i) - tick ) > 1e-5:
            i += 1
            
        # If digits = 0, then display integers
        if i == 0:
            tick_labels.append('%d' % tick)
        # Otherwise, display floats with the correct number of digits
        elif i == 1:
            tick_labels.append('%.1f' % tick)
        elif i == 2:
            tick_labels.append('%.2f' % tick)
        elif i == 3:
            tick_labels.append('%.3f' % tick)
        elif i == 4:
            tick_labels.append('%.4f' % tick)
        else:
            tick_labels.append('%.f' % tick)
        
    return tick_labels


def coupled_channel(channel):
    """
    Boolean value on whether the given channel is a coupled channel.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
    
    Returns
    -------
    boolean_value : bool
        True if the channel is coupled channel and false otherwise.
        
    """
    
    # List of coupled channels
    coupled_channels = ['3S1', '3P2', '3D3', '3F4', '3G5', '3H6', '3I7',
                        '3K8', '3L9', '3M10', '3N11', '3O12', '3Q13']

    # This is true or false
    boolean_value = channel in coupled_channels
    
    return boolean_value


def generator_label_conversion(generator, lambda_bd=0.00):
    """
    Converts an SRG generator string argument to a label for plotting purposes
    (e.g. generator = 'Wegner' gives r'$G = H_D$').
    
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
            
        # Otherwise, present with lambda_bd value
        else:
            lambda_bd_str = convert_number_to_string(lambda_bd)
            label = r'$G=H_{BD}$' + ' (%s fm'%lambda_bd_str + r'$^{-1}$' + ')'
            
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
    
    # The number of points in x_array less than x_max
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
    (e.g. lamb = 2 gives r'$\lambda=2 fm^-1$').
    
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
        
        # Make \lambda or \lambda_bd a string
        lamb_str = convert_number_to_string(lamb)
        
        # Label \lambda_bd
        if generator == 'Block-diag':
            label = r'$\Lambda_{BD}=%s$' % lamb_str + ' fm' + r'$^{-1}$'
            
        # Label \lambda
        else:
            label = r'$\lambda=%s$' % lamb_str + ' fm' + r'$^{-1}$'
            
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
    line_styles = [ 'solid', 'dashdot', 'dashed', 'dotted', (0, (5, 1) ), 
                    (0, (3, 1, 1, 1) ) ]
    
    try:
        
        return line_styles[curve_number]
    
    except IndexError:
        
        error = 'Curve number is too high to match with default line style.'
        suggestion = 'Manually assign styles to each curve.'
        print_statement = error + ' ' + suggestion
        
        print(print_statement)
        
        return 'solid'
    

def nuclei_label_conversion(nucleus):
    """
    Converts a nucleus string (e.g., 'C12') to a label with the mass number in
    the exponent appearing before the element (e.g., '^{12}C').

    Parameters
    ----------
    nucleus : str
        Specify nucleus (e.g., 'O16', 'Ca40', etc.)

    Returns
    -------
    nucleus_label : str
        Label for nucleus.

    """
    
    # Create mass number and element strings by looping over characters of
    # input nucleus string
    mass_number = ''
    element = ''
    
    for char in nucleus:
        
        try:
            
            number = int(char) # This gives a ValueError if char is a letter
            mass_number += str(number)
        
        except ValueError:
            
            element += char
            
    nucleus_label = r'$^{%s}$' % mass_number + element
            
    return nucleus_label


def replace_periods(file_name):
    """
    Replaces all periods in a file name with a 'p'. This is necessary for 
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
        New file name replacing periods with 'p'.
        
    """
    
    # Initialize new file name
    new_file_name = ''
    
    # Loop over each character in the original file name and append to the new
    # file name replacing periods with 'p'
    for letter in file_name:
        
        if letter == '.':
            
            letter = 'p'
            
        new_file_name += letter
        
    return new_file_name
    

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
              'xkcd:orange', 'xkcd:purple', 'xkcd:grey', 'xkcd:brown',
              'xkcd:pink', 'xkcd:turquoise', 'xkcd:olive', 'xkcd:gold',   
              'xkcd:indigo', 'xkcd:magenta', 'xkcd:tan', 'xkcd:crimson',
              'xkcd:navy', 'xkcd:lime', 'xkcd:plum', 'xkcd:chocolate',
              'xkcd:coral', 'xkcd:darkgreen', 'xkcd:khaki']
    
    try:
        
        return colors[curve_number]
    
    except IndexError:
        
        error = 'Curve number is too high to match with default color.'
        suggestion = 'Manually assign colors to each curve.'
        print_statement = error + ' ' + suggestion
        
        print(print_statement)
        
        return 'xkcd:black'
