#!/usr/bin/env python3

"""
File: figures_graphics.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 3, 2019

Useful functions for plotting figures with matplotlib.

Last update: March 17, 2022

"""

# To-do: Currently keeping setup_rc_params as it was. But there are probably
# many more things to mess with.
# To-do: Could add y_array argument to interpolate matrix. I'm assuming it's
# a matrix M(x, x) currently.
# To-do: Couldn't you just input x_array and M_matrix, interpolate, evaluate
# on new mesh (np.linspace(x_min, x_max, ntot)), and return the new arrays?

# Python imports
import matplotlib as mpl
import numpy as np
from scipy.interpolate import interp1d, interp2d


def setup_rc_params(presentation=False):
    """
    Set matplotlib's rc parameters for figures. Run this function in Jupyter
    notebooks or Python scripts when you want fancy-looking plots.
        
    Parameters
    ----------
    presentation : bool, optional
        Option to enlarge font sizes for presentations.
    
    """

    if presentation:
        fontsize = 14
    else:
        fontsize = 9
    black = 'k'

    mpl.rcdefaults()  # Set to defaults

    mpl.rc('text', usetex=True)
    mpl.rcParams['font.size'] = fontsize
    # This will give an error if you don't have LaTeX
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['font.family'] = 'serif'

    # mpl.rcParams['axes.labelsize'] = fontsize
    mpl.rcParams['axes.edgecolor'] = black
    # mpl.rcParams['axes.xmargin'] = 0
    mpl.rcParams['axes.labelcolor'] = black
    # mpl.rcParams['axes.titlesize'] = fontsize

    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['xtick.labelsize'] = fontsize
    mpl.rcParams['ytick.labelsize'] = fontsize
    mpl.rcParams['xtick.color'] = black
    mpl.rcParams['ytick.color'] = black
    # Make the ticks thin enough to not be visible at the limits of the plot
    mpl.rcParams['xtick.major.width'] = mpl.rcParams['axes.linewidth'] * 0.95
    mpl.rcParams['ytick.major.width'] = mpl.rcParams['axes.linewidth'] * 0.95
    # The minor ticks are little too small, make them both bigger.
    mpl.rcParams['xtick.minor.size'] = 2.4  # Default 2.0
    mpl.rcParams['ytick.minor.size'] = 2.4
    mpl.rcParams['xtick.major.size'] = 3.9  # Default 3.5
    mpl.rcParams['ytick.major.size'] = 3.9
    
    # Added by A.T.
    # Puts tick marks (not labels) on top and right axes
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    
    ppi = 72  # points per inch
    # dpi = 150
    # mpl.rcParams['figure.titlesize'] = fontsize
    mpl.rcParams['figure.dpi'] = 150  # To show up reasonably in notebooks
    mpl.rcParams['figure.constrained_layout.use'] = False
    # 0.02 and 3 points are the defaults:
    # can be changed on a plot-by-plot basis using fig.set_constrained_layout_pads()
    mpl.rcParams['figure.constrained_layout.wspace'] = 0.0
    mpl.rcParams['figure.constrained_layout.hspace'] = 0.0
    mpl.rcParams['figure.constrained_layout.h_pad'] = 3. / ppi  # 3 points
    mpl.rcParams['figure.constrained_layout.w_pad'] = 3. / ppi

    #  mpl.rcParams['legend.title_fontsize'] = fontsize
    #  mpl.rcParams['legend.fontsize'] = fontsize
    # Inherits from axes.edgecolor, to match
    mpl.rcParams['legend.edgecolor'] = 'inherit'  
    # Set facecolor with its own alpha, so edgecolor is unaffected
    mpl.rcParams['legend.fancybox'] = True
    mpl.rcParams['legend.facecolor'] = (1, 1, 1, 0.6)
    # mpl.rcParams['legend.borderaxespad'] = 0.8
    # Do not set overall alpha (affects edgecolor)
    # Handled by facecolor above
    mpl.rcParams['legend.framealpha'] = None
    # This is for legend edgewidth, since it does not have its own option
    mpl.rcParams['patch.linewidth'] = 0.8
    mpl.rcParams['hatch.linewidth'] = 0.5

    # bbox = 'tight' can distort the figure size when saved 
    # (that's its purpose)
    # mpl.rc('savefig', transparent=False, bbox='tight', pad_inches=0.04,
    #        dpi=350, format='png')
    # mpl.rc('savefig', transparent=False, bbox=None, dpi=400, format='png')
    mpl.rc('savefig', bbox='tight', dpi=400)
    

def interpolate_matrix(x_array, M, x_max, ntot=500):
    """
    Interpolate matrix over given array for high resolution contour plots.

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
    ntot = len(list(filter(None, bool_list)))
    
    # Resize x_array and M to size that you want to plot
    x_array_resized = x_array[:ntot]
    M_resized = M[:ntot, :ntot]
    
    # Use interp2d to interpolate the truncated matrix
    M_func = interp2d(x_array_resized, x_array_resized, M_resized)

    # New x array for interpolation (dimension ntot x 1)
    x_array_new = np.linspace(x_array_resized[0], x_array_resized[-1], ntot)
    
    # New matrix (dimension ntot x ntot)
    M_new = M_func(x_array_new, x_array_new) 
    
    return x_array_new, M_new


def interpolate_vector(x_array, y_array, x_max, ntot=500, order='linear'):
    """
    Interpolate vector given array for high resolution plots.
    
    Parameters
    ----------
    x_array : 1-D ndarray
        Array of x points where the vector y_array = y(x).
    y_array : 1-D ndarray
        Vector to be interpolated.
    x_max : float
        Maximum value of x. This functions interpolates y up to this value.
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
    n = len(list(filter(None, bool_list)))
    
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


def line_styles(curve_number):
    """
    Default line styles for plotting figures with many curves. Styles are set
    by the ordering of curves on the plot. Do not use this function for plots
    with more than six curves. An error message will display when a number
    higher than six is used and the function will return 'solid'.
    
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
    line_styles = [
        'solid', 'dashdot', 'dashed', 'dotted', (0, (5, 1) ),
        (0, (3, 1, 1, 1) )
    ]
    
    try:
        
        return line_styles[curve_number]
    
    except IndexError:
        
        error = 'Curve number is too high to match with default line style.'
        suggestion = 'Manually assign styles to each curve.'
        print_statement = error + ' ' + suggestion
        
        print(print_statement)
        
        return 'solid'
    

def xkcd_colors(curve_number):
    """
    Default curve colors for plotting figures with many curves. Colors are set
    by the ordering of curves on the plot. Do not use this function for plots
    with more than 23 curves. An error message will display when a number
    higher than 23 is used and the function will return 'xkcd:black'.
    
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
    
    colors = [
        'xkcd:black', 'xkcd:red', 'xkcd:blue', 'xkcd:green', 'xkcd:orange',
        'xkcd:purple', 'xkcd:grey', 'xkcd:brown', 'xkcd:pink',
        'xkcd:turquoise', 'xkcd:olive', 'xkcd:gold', 'xkcd:indigo',
        'xkcd:magenta', 'xkcd:tan', 'xkcd:crimson', 'xkcd:navy', 'xkcd:lime',
        'xkcd:plum', 'xkcd:chocolate', 'xkcd:coral', 'xkcd:darkgreen',
        'xkcd:khaki'
    ]
    
    try:
        
        return colors[curve_number]
    
    except IndexError:
        
        error = 'Curve number is too high to match with default color.'
        suggestion = 'Manually assign colors to each curve.'
        print_statement = error + ' ' + suggestion
        
        print(print_statement)
        
        return 'xkcd:black'