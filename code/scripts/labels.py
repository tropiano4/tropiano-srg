#!/usr/bin/env python3

"""
File: labels.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 3, 2019

Useful functions for labeling figures and file names.

Last update: April 27, 2022

"""

# To-do: Change % operator to f-strings?

# Python imports
import numpy as np

# Imports from A.T. codes
from .tools import coupled_channel, convert_number_to_string


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
            
        label = r'$^{%s}{\rm %s}_{%s}-^{%s}{\rm %s}_{%s}$' % (S, L, J,
                                                              S, L2, J)
        
    else:
        
        label = r'$^{%s}{\rm %s}_{%s}$' % (S, L, J)
    
    return label


def label_generator(generator, lambda_bd=0.00):
    """
    Converts an SRG generator string argument to a label.
    (e.g. generator = 'Wegner' gives r'$G = H_D$').
    
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
        if lambda_bd == 0.00:
            label = r'$G=H_{BD}$'
            
        # Otherwise, present with lambda_bd value
        else:
            lambda_bd_str = convert_number_to_string(lambda_bd)
            label = r'$G=H_{BD}$' + f' ({lambda_bd_str} fm' + r'$^{-1}$' + ')'
            
    elif generator == 'Wegner':
        label = r'$G=H_{D}$'
        
    elif generator == 'T':
        label = r'$G=T_{rel}$'

    return label


def label_kvnn(kvnn, full_label=True):
    """
    Converts a kvnn number to a label (e.g. kvnn = 6 gives 'AV18').
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    full_label : bool, optional
        For some labels, there is a shorter version. Set full_label = False
        for the shorter version. For example, kvnn = 902 gives
        '\Lambda = 20 fm^-1' normally, but the shorter version is '20 fm^-1'.
        
    Returns
    -------
    label : str
        Label for the potential.
        
    Notes
    -----
    It might be worth adding the cutoff to the full_label version for all the
    potentials (e.g. EM N3LO 500 MeV).
        
    """

    # Paris
    if kvnn == 1:
        label = 'Paris'

    # Bonn
    elif kvnn == 2:
        label = 'Bonn'

    # Reid93 potential
    elif kvnn == 3:
        label = 'Reid93'

    # Nijmegen I potential
    elif kvnn == 4:
        label = 'Nijmegen I'

    # Nijmegen II potential
    elif kvnn == 5:
        label = 'Nijmegen II'

    # Argonne v18
    elif kvnn == 6:
        label = 'AV18'

    # CD-Bonn
    elif kvnn == 7:
        label = 'CD-Bonn'
        
    # Entem/Machleidt N3LO (500 MeV cutoff)   
    elif kvnn == 10:
        if full_label:
            label = 'EM N' + r'$^3$' + 'LO 500 MeV'
        else:
            label = 'EM N' + r'$^3$' + 'LO'
        
    # EMN N4LO (450, 500, 550 MeV cutoffs)
    elif kvnn in [74, 79, 84]:
        if full_label:
            if kvnn == 74:
                label = 'EMN N' + r'$^4$' + 'LO 450 MeV'
            elif kvnn == 79:
                label = 'EMN N' + r'$^4$' + 'LO 500 MeV'
            elif kvnn == 84:
                label = 'EMN N' + r'$^4$' + 'LO 550 MeV'
        else:
            label = 'EMN N' + r'$^4$' + 'LO'
        
    # SMS N3LO (400, 450, 500 MeV cutoffs)
    elif kvnn in [105, 106, 107]:
        if full_label:
            if kvnn == 105:
                label = 'SMS N' + r'$^3$' + 'LO 400 MeV'
            elif kvnn == 106:
                label = 'SMS N' + r'$^3$' + 'LO 450 MeV'
            elif kvnn == 107:
                label = 'SMS N' + r'$^3$' + 'LO 500 MeV'
        else:
            label = 'SMS N' + r'$^3$' + 'LO'
    # SMS N4LO (400, 450, 500, 550 MeV cutoffs)
    elif kvnn in [110, 111, 112, 113]:
        if full_label:
            if kvnn == 110:
                label = 'SMS N' + r'$^4$' + 'LO 400 MeV'
            elif kvnn == 111:
                label = 'SMS N' + r'$^4$' + 'LO 450 MeV'
            elif kvnn == 112:
                label = 'SMS N' + r'$^4$' + 'LO 500 MeV'
            elif kvnn == 113:
                label = 'SMS N' + r'$^4$' + 'LO 550 MeV'
        else:
            label = 'SMS N' + r'$^4$' + 'LO'
        
    # Gezerlis N2LO (1 and 1.2 fm cutoff)
    elif kvnn in [222, 224]:
        if full_label:
            if kvnn == 222:
                label = 'GT+ N' + r'$^2$' + 'LO 1 fm'
            elif kvnn == 224:
                label = 'GT+ N' + r'$^2$' + 'LO 1.2 fm'
        else:
            label = 'GT+ N' + r'$^2$' + 'LO'
        
    # Wendt LO non-local potentials
    elif kvnn == 900:  # Cutoff 4 fm^-1
        if full_label:
            label = r'$\Lambda = 4$' + ' fm' + r'$^{-1}$'
        else:
            label = '4 fm' + r'$^{-1}$'
    elif kvnn == 901:  # Cutoff 9 fm^-1
        if full_label:
            label = r'$\Lambda = 9$' + ' fm' + r'$^{-1}$'
        else:
            label = '9 fm' + r'$^{-1}$'
    elif kvnn == 902:  # Cutoff 20 fm^-1
        if full_label:
            label = r'$\Lambda = 20$' + ' fm' + r'$^{-1}$'
        else:
            label = '20 fm' + r'$^{-1}$'

    return label


def label_lambda(lamb, generator='Wegner'):
    """
    Converts a lambda evolution parameter to a label
    (e.g. lamb = 2 gives r'$\lambda=2 fm^-1$').
    
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
            label = r'$\Lambda_{BD}=\infty$' + ' fm' + r'$^{-1}$'
            
        # Label \lambda
        else:
            label = r'$\lambda=\infty$' + ' fm' + r'$^{-1}$'
            
    # Finite lamb
    else:
        
        # Make \lambda or \Lambda_BD a string
        lamb_str = convert_number_to_string(lamb)
        
        # Label \Lambda_BD
        if generator == 'Block-diag':
            label = r'$\Lambda_{BD}=%s$' % lamb_str + ' fm' + r'$^{-1}$'
            
        # Label \lambda
        else:
            label = r'$\lambda=%s$' % lamb_str + ' fm' + r'$^{-1}$'
            
    return label
    

def label_nucleus(nucleus_name):
    """
    Converts a nucleus string (e.g., 'C12') to a label with the mass number in
    the exponent appearing before the element (e.g., '^{12}C').

    Parameters
    ----------
    nucleus_name : str
        Specify the nucleus (e.g., 'O16', 'Ca40', etc.)

    Returns
    -------
    nucleus_label : str
        Label for nucleus.

    """
    
    # Create mass number and element strings by looping over characters of
    # input nucleus string
    mass_number = ''
    element = ''
    
    for char in nucleus_name:
        
        # This gives a ValueError if char is a letter
        try:
            
            number = int(char) 
            mass_number += str(number)
        
        except ValueError:
            
            element += char
            
    nucleus_label = r'$^{%s}$' % mass_number + element
            
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
        String for figure label (e.g., '1s_{\frac{1}{2}}').
        
    """
    
    numerator = 2*int(sp_state[-3]) + 1
    denominator = 2
    
    return r'$%s_{%d/%d}$' % (sp_state[:2], numerator, denominator)


def label_ticks(ticks):
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
    
    
    tick_labels = []
    # # This won't work for some reason
    # for tick in ticks:
    #     tick_labels.append(convert_number_to_string(tick))

    # Keep track of the maximum number of digits to be displayed i
    for tick in ticks:
        
        i = 0
        while abs(round(tick,i) - tick) > 1e-5:
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


def replace_periods(file_name):
    """
    Replaces all periods in a file name with a 'p'. This is necessary for 
    adding figures to LaTeX files which don't like periods unless they specify
    a file type. For this reason, do not include the file type extension in 
    the file name (i.e. .jpg, .pdf, etc.)
    
    Parameters
    ----------
    file_name : str
        Original name of the file including periods. Should not include the
        extension (e.g., .pdf).
    
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