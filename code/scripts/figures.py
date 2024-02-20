#!/usr/bin/env python3

"""
File: figures.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 3, 2019

Useful functions for plotting figures with matplotlib and adding labels.

Last update: November 29, 2022

"""

# Python imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d, interp2d

# Imports from A.T. codes
from .tools import (
    coupled_channel, convert_l_to_string, convert_number_to_string
)


def set_rc_parameters():
    """
    Set matplotlib's RC parameters for figures.

    Notes
    -----
    * To get usetex=True to work, you need LaTeX, dvipng, and
      Ghostscript>9.0. Will additionally cause issues if you don't have 
      add-on things (e.g., "cm-super: CM-Super family of fonts"). I got
      around this by downloading TeX Live.
    * I've also commented out some axes and title fontsizes to customize on a
      plot-to-plot basis.

    """

    # Defaults for font
    fontsize = 14
    black = 'k'

    mpl.rcdefaults()  # Set to defaults

    # Set LaTeX fonts (this will give an error if you don't have LaTeX)
    mpl.rc('text', usetex=True)
    mpl.rcParams['font.size'] = fontsize
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['font.family'] = 'serif'

    # Settings for axes
    # mpl.rcParams['axes.labelsize'] = fontsize
    mpl.rcParams['axes.edgecolor'] = black
    # mpl.rcParams['axes.xmargin'] = 0
    mpl.rcParams['axes.labelcolor'] = black
    # mpl.rcParams['axes.titlesize'] = fontsize

    # Settings for ticks
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['xtick.labelsize'] = fontsize
    mpl.rcParams['ytick.labelsize'] = fontsize
    mpl.rcParams['xtick.color'] = black
    mpl.rcParams['ytick.color'] = black
    # Make the ticks thin to not be visible at the limits of the plot
    mpl.rcParams['xtick.major.width'] = mpl.rcParams['axes.linewidth'] * 0.95
    mpl.rcParams['ytick.major.width'] = mpl.rcParams['axes.linewidth'] * 0.95
    # The minor ticks are little too small, make them both bigger.
    mpl.rcParams['xtick.minor.size'] = 2.4  # Default 2.0
    mpl.rcParams['ytick.minor.size'] = 2.4
    mpl.rcParams['xtick.major.size'] = 3.9  # Default 3.5
    mpl.rcParams['ytick.major.size'] = 3.9
    # Put ticks on top and right axes as well
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    # Settings for overall figure
    ppi = 72  # points per inch
    # mpl.rcParams['figure.titlesize'] = fontsize
    mpl.rcParams['figure.dpi'] = 150  # To show up reasonably in notebooks
    mpl.rcParams['figure.constrained_layout.use'] = False
    # 0.02 and 3 points are the defaults:
    # can be changed on a plot-by-plot basis using
    # fig.set_constrained_layout_pads()
    mpl.rcParams['figure.constrained_layout.wspace'] = 0.0
    mpl.rcParams['figure.constrained_layout.hspace'] = 0.0
    mpl.rcParams['figure.constrained_layout.h_pad'] = 3. / ppi  # 3 points
    mpl.rcParams['figure.constrained_layout.w_pad'] = 3. / ppi

    # Settings for legend
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

    # Settings for saving figure
    mpl.rc('savefig', bbox='tight', dpi=400)


def interpolate_matrix(
        x_array, y_array, M_array, x_max, y_max, x_tot=500, y_tot=500):
    """
    Interpolate matrix over given array for high resolution contour plots.

    Parameters
    ----------
    x_array : 1-D ndarray
        Array of x points where the matrix M = M(x, y).
    y_array : 1-D ndarray
        Array of y points where the matrix M = M(x, y).
    M_array : 2-D ndarray
        Matrix to be interpolated.
    x_max : float
        Return points up to this maximum value of x.
    y_max : float
        Return points up to this maximum value of y.
    x_tot : int, optional
        Desired number of points in x.
    y_tot : int, optional
        Desired number of points in y.
    
    Returns
    -------
    x_array_new : 1-D ndarray
        Dense array of x points up to x_max with shape (x_tot,).
    y_array_new : 1-D ndarray
        Dense array of y points up to y_max with shape (y_tot,).
    M_array_new : 2-D ndarray
        Dense matrix evaluated up to (x_max, y_max) with
        shape (x_tot, y_tot).
    
    """

    # Use interp2d to interpolate the truncated matrix
    M_func = interp2d(x_array, y_array, M_array)

    # Dense x and y arrays for interpolation
    x_array_dense = np.linspace(x_array[0], x_max, x_tot)
    y_array_dense = np.linspace(y_array[0], y_max, y_tot)

    # Dense matrix
    M_array_dense = M_func(x_array_dense, y_array_dense)

    return x_array_dense, y_array_dense, M_array_dense


def interpolate_vector(x_array, v_array, x_max, x_tot=500, order='linear'):
    """
    Interpolate vector over given array for high resolution plots.

    Parameters
    ----------
    x_array : 1-D ndarray
        Array of x points where the vector v = v(x).
    v_array : 1-D ndarray
        Vector to be interpolated.
    x_max : float
        Return points up to this maximum value of x.
    x_tot : int, optional
        Desired number of points.
    order : str, optional
        Order for interpolation. Default is 'linear'.
    
    Returns
    -------
    x_array_new : 1-D ndarray
        Dense array of x points up to x_max with shape (x_tot,).
    v_array_new : 1-D ndarray
        Dense vector evaluated up to x_max with shape (x_tot,).
    
    """

    # Use interp1d to interpolate the truncated matrix
    v_func = interp1d(x_array, v_array, kind=order)

    # Dense x array for interpolation
    x_array_dense = np.linspace(x_array[0], x_max, x_tot)

    # Dense vector
    v_array_dense = v_func(x_array_dense)

    return x_array_dense, v_array_dense


def line_styles(curve_number):
    """
    Default line styles for plotting figures with many curves. Styles are
    set by the ordering of curves on the plot. Do not use this function
    for plots with more than six curves. An error message will display
    when a number higher than six is used and the function will return
    'solid'.

    Parameters
    ----------
    curve_number : int
        The curve number being assigned a line style.

    Returns
    -------
    out : str
        The line style argument to be used in plotting functions. E.g.,
    
            plt.plot(x, y_1, linestyle=line_style(0))
            plt.plot(x, y_2, linestyle=line_style(3))
    
        plots two curves as solid and dotted lines.

    """

    # Note, the last two are 'densely dashed' and 'densely dashdotted'
    styles = (
        'solid', 'dashdot', 'dashed', 'dotted', (0, (5, 1)),
        (0, (3, 1, 1, 1))
    )

    try:

        return styles[curve_number]

    except IndexError:

        error_message = ("Curve number is too high to match with default "
                         "line style. Manually assign styles to each "
                         "curve.")
        print(error_message)

        return 'solid'


def tab10_colors(curve_number):
    """
    Default curve colors for plotting figures with many curves. Colors are
    set by the ordering of curves on the plot. These colors are the default
    matplotlib colors.

    Parameters
    ----------
    curve_number : int
        The curve number being assigned a color.

    Returns
    -------
    out : str
        The color argument to be used in plotting functions.

    """

    cmap = plt.get_cmap("tab10")

    return cmap(curve_number)


def xkcd_colors(curve_number):
    """
    Default curve colors for plotting figures with many curves. Colors are
    set by the ordering of curves on the plot. Do not use this function
    for plots with more than 23 curves. An error message will display when
    a number higher than 23 is used and the function will return
    'xkcd:black'.

    Parameters
    ----------
    curve_number : int
        The curve number being assigned a color.

    Returns
    -------
    out : str
        The color argument to be used in plotting functions. E.g.,
    
            plt.plot(x, y_1, color=xkcd_colors(0))
            plt.plot(x, y_2, color=xkcd_colors(1))
    
        plots two curves as black and red lines.

    """

    colors = (
        'xkcd:black', 'xkcd:red', 'xkcd:blue', 'xkcd:green', 'xkcd:orange',
        'xkcd:purple', 'xkcd:grey', 'xkcd:brown', 'xkcd:pink',
        'xkcd:turquoise', 'xkcd:olive', 'xkcd:gold', 'xkcd:indigo',
        'xkcd:magenta', 'xkcd:tan', 'xkcd:crimson', 'xkcd:navy',
        'xkcd:lime', 'xkcd:plum', 'xkcd:chocolate', 'xkcd:coral',
        'xkcd:darkgreen', 'xkcd:khaki'
    )

    try:

        return colors[curve_number]

    except IndexError:

        error_message = ("Curve number is too high to match with default "
                         "color. Manually assign colors to each curve.")
        print(error_message)

        return 'xkcd:black'


def label_channel(channel, label_coupled_channel=True):
    """
    Converts a channel string argument to a label.
    
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
        else:
            raise RuntimeError("Orbital angular momentum is too high.")

        label = fr"$^{{{S}}}{{\rm {L}}}_{{{J}}}-^{{{S}}}{{\rm {L2}}}_{{{J}}}$"

    else:

        label = fr"$^{{{S}}}{{\rm {L}}}_{{{J}}}$"

    return label


def label_generator(generator, lambda_bd=None):
    """
    Converts an SRG generator string argument to a label.
    (e.g. generator = 'Wegner' gives r"$G = H_D$").
    
    Parameters
    ----------
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lambda_bd : float, optional
        SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        
    Returns
    -------
    label : str
        Label for the generator.
        
    """

    # Block-diagonal
    if generator == 'Block-diag':

        # In the case, lambda_bd is not entered as an argument, the function
        # returns only G = H_{BD}
        if lambda_bd is None:
            label = r"$G=H_{\rm{BD}}$"

        # Otherwise, present with lambda_bd value
        else:
            lambda_bd_str = convert_number_to_string(lambda_bd)
            label = (r"$G=H_{\rm{BD}}$" + f" ({lambda_bd_str} fm" + r"$^{-1}$"
                     + ")")

    elif generator == 'Wegner':
        label = r"$G=H_{\rm{D}}$"

    elif generator == 'T':
        label = r"$G=T_{\rm{rel}}$"

    else:
        raise RuntimeError("Not a valid generator name.")

    return label


def label_kvnn(kvnn, full_label=True):
    """
    Converts a kvnn number to a label (e.g. kvnn = 6 gives "AV18").
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    full_label : bool, optional
        For some labels, there is a shorter version. Set
        full_label = False for the shorter version. For example,
        kvnn = 902 gives "\Lambda = 20 fm^-1" normally, but the shorter
        version is "20 fm^-1".
        
    Returns
    -------
    output : str
        Label for the potential.
        
    Notes
    -----
    It might be worth adding the cutoff to the full_label version for all the
    potentials (e.g. EM N3LO 500 MeV).
        
    """

    # Paris
    if kvnn == 1:
        return "Paris"

    # Bonn
    elif kvnn == 2:
        return "Bonn"

    # Reid93 potential
    elif kvnn == 3:
        return "Reid93"

    # Nijmegen I potential
    elif kvnn == 4:
        return "Nijmegen I"

    # Nijmegen II potential
    elif kvnn == 5:
        return "Nijmegen II"

    # Argonne v18
    elif kvnn == 6:
        return "AV18"

    # CD-Bonn
    elif kvnn == 7:
        return "CD-Bonn"

    # Entem/Machleidt N3LO (500 MeV cutoff)   
    elif kvnn == 10:
        if full_label:
            return "EM N" + r"$^3$" + "LO 500 MeV"
        else:
            return "EM N" + r"$^3$" + "LO"

    # EMN N4LO (450, 500, 550 MeV cutoffs)
    elif kvnn in [74, 79, 84]:
        if full_label:
            if kvnn == 74:
                return "EMN N" + r"$^4$" + "LO 450 MeV"
            elif kvnn == 79:
                return "EMN N" + r"$^4$" + "LO 500 MeV"
            elif kvnn == 84:
                return "EMN N" + r"$^4$" + "LO 550 MeV"
        else:
            return "EMN N" + r"$^4$" + "LO"

    # SMS N3LO (400, 450, 500 MeV cutoffs)
    elif kvnn in [105, 106, 107]:
        if full_label:
            if kvnn == 105:
                return "SMS N" + r"$^3$" + "LO 400 MeV"
            elif kvnn == 106:
                return "SMS N" + r"$^3$" + "LO 450 MeV"
            elif kvnn == 107:
                return "SMS N" + r"$^3$" + "LO 500 MeV"
        else:
            return "SMS N" + r"$^3$" + "LO"
    # SMS N4LO (400, 450, 500, 550 MeV cutoffs)
    elif kvnn in [110, 111, 112, 113]:
        if full_label:
            if kvnn == 110:
                return "SMS N" + r"$^4$" + "LO 400 MeV"
            elif kvnn == 111:
                return "SMS N" + r"$^4$" + "LO 450 MeV"
            elif kvnn == 112:
                return "SMS N" + r"$^4$" + "LO 500 MeV"
            elif kvnn == 113:
                return "SMS N" + r"$^4$" + "LO 550 MeV"
        else:
            return "SMS N" + r"$^4$" + "LO"

    # Gezerlis N2LO (1 and 1.2 fm cutoff)
    elif kvnn in [222, 224]:
        if full_label:
            if kvnn == 222:
                return "GT+ N" + r"$^2$" + "LO 1 fm"
            elif kvnn == 224:
                return "GT+ N" + r"$^2$" + "LO 1.2 fm"
        else:
            return "GT+ N" + r"$^2$" + "LO"

    # Wendt LO non-local potentials
    elif kvnn == 900:  # Cutoff 4 fm^-1
        return r"$\Lambda = 4$" + " fm" + r"$^{-1}$"
    elif kvnn == 901:  # Cutoff 9 fm^-1
        return r"$\Lambda = 9$" + " fm" + r"$^{-1}$"
    elif kvnn == 902:  # Cutoff 20 fm^-1
        return r"$\Lambda = 20$" + " fm" + r"$^{-1}$"

    else:
        raise RuntimeError("Not a valid kvnn number.")


def label_lambda(lamb, generator='Wegner'):
    """
    Converts a lambda evolution parameter to a label
    (e.g. lamb = 2 gives r"$\lambda=2 fm^-1$").
    
    Parameters
    ----------
    lamb : float
        SRG evolution parameter \lambda or \Lambda_BD [fm^-1].
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'. Determines whether
        lamb is referring to \lambda or \Lambda_BD.
        
    Returns
    -------
    label : str
        Label for the lambda.
        
    """

    if lamb == np.inf:

        # Label \Lambda_BD
        if generator == 'Block-diag':
            label = r"$\Lambda_{\rm{BD}}=\infty$" + " fm" + r"$^{-1}$"

        # Label \lambda
        else:
            label = r"$\lambda=\infty$" + " fm" + r"$^{-1}$"

    # Finite lamb
    else:

        # Make \lambda or \Lambda_BD a string
        lamb_str = convert_number_to_string(lamb)

        # Label \Lambda_BD
        if generator == 'Block-diag':
            label = r"$\Lambda_{\rm{BD}}=%s$" % lamb_str + " fm" + r"$^{-1}$"

        # Label \lambda
        else:
            label = rf"$\lambda={lamb_str}$" + " fm" + r"$^{-1}$"

    return label


def label_nucleus(nucleus_name):
    """
    Converts a nucleus string (e.g., 'C12') to a label with the mass 
    number in the exponent appearing before the element (e.g., "^{12}C").

    Parameters
    ----------
    nucleus_name : str
        Specify the nucleus (e.g., 'O16', 'Ca40', etc.)

    Returns
    -------
    nucleus_label : str
        Label for nucleus.

    """

    # Create mass number and element strings by looping over characters
    # of input nucleus string
    mass_number = ''
    element = ''

    for char in nucleus_name:

        # This gives a ValueError if char is a letter
        try:

            number = int(char)
            mass_number += str(number)

        except ValueError:

            element += char

    nucleus_label = rf"$^{{{mass_number}}}$" + element

    return nucleus_label


def label_sp_state(sp_state):
    """
    Convert single-particle state to a label.
    
    Parameters
    ----------
    sp_state : str
        s.p. state as a string (e.g., '1s0p5').
        
    Returns
    -------
    output : str
        String for figure label (e.g., "1s_{\frac{1}{2}}").
        
    """

    numerator = 2 * int(sp_state[-3]) + 1
    denominator = 2

    return rf"${sp_state[:2]}_{{{numerator}/{denominator}}}$"


def label_nlj_state(n, l, j):
    """
    Convert single-particle state to a label.
    
    Parameters
    ----------
    n : int
        Principal quantum number n = 1, 2, 3, ...
    l : int
        Orbital angular momentum l = 0, 1, 2, ...
    j : float
        Total angular momentum j = 1/2, 3/2, 5/2, ...
        
    Returns
    -------
    output : str
        String for figure label (e.g., "1s_{\frac{1}{2}}").
        
    """
    
    l_str = convert_l_to_string(l)  # E.g., 's', 'p', 'd', ...
    numerator = 2*int(j) + 1
    denominator = 2

    return rf"${n}{l_str}_{{{numerator}/{denominator}}}$"


def label_ticks(ticks):
    """
    Converts axes or colorbar ticks to string formatted labels displaying
    the correct number of digits.
    
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

    tick_labels = []
    # # This won't work for some reason
    # for tick in ticks:
    #     tick_labels.append(convert_number_to_string(tick))

    # Keep track of the maximum number of digits to be displayed i
    for tick in ticks:

        i = 0
        while abs(round(tick, i) - tick) > 1e-5:
            i += 1

        # If digits = 0, then display integers
        if i == 0:
            tick_labels.append("%d" % tick)
        # Otherwise, display floats with the correct number of digits
        elif i == 1:
            tick_labels.append("%.1f" % tick)
        elif i == 2:
            tick_labels.append("%.2f" % tick)
        elif i == 3:
            tick_labels.append("%.3f" % tick)
        elif i == 4:
            tick_labels.append("%.4f" % tick)
        else:
            tick_labels.append("%.f" % tick)

    return tick_labels