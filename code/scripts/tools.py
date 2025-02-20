#!/usr/bin/env python3

"""
File: tools.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 17, 2021

This script contains several useful functions for general purposes.

Last update: February 8, 2023

"""

# Python imports
import numpy as np


def build_coupled_channel_matrix(M11, M12, M21, M22):
    """
    Builds a coupled-channel matrix.

    Parameters
    ----------
    M11 : 2-D ndarray
        Upper left block of full matrix (e.g., 3S1-3S1).
    M12 : 2-D ndarray
        Upper right block of full matrix (e.g., 3S1-3D1).
    M21 : 2-D ndarray
        Lower left block of full matrix (e.g., 3D1-3S1).
    M22 : 2-D ndarray
        Lower right block of full matrix (e.g., 3D1-3D1).

    Returns
    -------
    M_full : 2-D ndarray
        Full coupled-channel matrix. Same units as sub-block inputs. If
        len(M11) = ntot, then len(M_full) = 2*ntot.

    """

    # Build coupled-channel matrix
    M_full = np.vstack((np.hstack((M11, M12)), np.hstack((M21, M22))))

    return M_full


def decompose_coupled_channel_matrix(M_full):
    """
    Extracts the sub-blocks of a coupled-channel matrix.
    
    Parameters
    ----------
    M_full : 2-D ndarray
        Full coupled-channel matrix. Must have even length.
    
    Returns
    -------
    M11 : 2-D ndarray
        Upper left block of full matrix (e.g., 3S1-3S1).
    M12 : 2-D ndarray
        Upper right block of full matrix (e.g., 3S1-3D1).
    M21 : 2-D ndarray
        Lower left block of full matrix (e.g., 3D1-3S1).
    M22 : 2-D ndarray
        Lower right block of full matrix (e.g., 3D1-3D1).
    
    """

    # Full length
    Ntot = len(M_full)

    # Sub-blocks will have exactly half the length of M_full
    ntot = int(Ntot / 2)

    # Get each sub-block
    M11 = M_full[:ntot, :ntot]
    M12 = M_full[:ntot, ntot:]
    M21 = M_full[ntot:, :ntot]
    M22 = M_full[ntot:, ntot:]

    return M11, M12, M21, M22


def channel_L_value(channel):
    """
    Returns the L value associated with the given partial wave channel.

    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').

    Returns
    -------
    L : int
        Total orbital angular momentum associated with partial wave channel.

    """

    # This gives 'S', 'P', etc.
    channel_letter = channel[1]

    if channel_letter == 'S':
        return 0
    elif channel_letter == 'P':
        return 1
    elif channel_letter == 'D':
        return 2
    elif channel_letter == 'F':
        return 3
    elif channel_letter == 'G':
        return 4
    elif channel_letter == 'H':
        return 5
    elif channel_letter == 'I':
        return 6
    elif channel_letter == 'K':
        return 7
    elif channel_letter == 'L':
        return 8
    elif channel_letter == 'M':
        return 9
    elif channel_letter == 'N':
        return 10
    elif channel_letter == 'O':
        return 11
    elif channel_letter == 'Q':
        return 12
    else:
        print("Input channel is outside the range of this function.")
        return None

    
def convert_l_to_string(l):
    """
    Returns the spectroscopic notation of the orbital angular momentum value
    l (e.g., l = 2 returns 'd').

    Parameters
    ----------
    l : int
        Orbital angular momentum of the single-particle (s.p.) state.

    Returns
    -------
    output : str
        Spectroscopic notation of s.p. state orbital angular momentum.

    """

    if l == 0:
        return 's'
    elif l == 1:
        return 'p'
    elif l == 2:
        return 'd'
    elif l == 3:
        return 'f'
    elif l == 4:
        return 'g'
    elif l == 5:
        return 'h'
    elif l == 6:
        return 'i'
    else:
        raise RuntimeError(
            "Input l value is outside the range of this function."
        )

        
def coupled_channel(channel):
    """
    Boolean value on whether the given channel is a coupled-channel.
    
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
    coupled_channels = (
        '3S1', '3D1', '3P2', '3F2', '3D3', '3G3', '3F4', '3H4', '3G5', '3I5',
        '3H6', '3K6', '3I7', '3L7', '3K8', '3M8', '3L9', '3N9', '3M10', '3O10',
        '3N11', '3Q11'
    )

    boolean_value = channel in coupled_channels

    return boolean_value


def convert_number_to_string(number):
    """
    Gives the input number with the correct amount of digits meaning the
    string will show '1.35' not '1.35000000023'.

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

    # First check if the number is an integer
    if number == int(number):
        return str(int(number))

    else:

        # Loop over i until the rounded number matches the input number,
        # then return the string of the rounded number
        i = 0
        while round(number, i) != number:
            i += 1

        return str(round(number, i))


def find_index(x, x_array):
    """
    Finds the index of an element nearest to the input value x in an array.
    For example, say k = 3.0 fm^-1 and k_array does not contain 3.0 exactly.
    Then this function would return an index corresponding to the k value
    nearest to 3.0 in k_array (e.g. k_array[index] = 2.98 fm^-1).
    
    Parameters
    ----------
    x : float
        Any value associated with the input array.
    x_array : 1-D ndarray
        Array of x values. Units of this must match units of x.
        
    Returns
    -------
    index : int
        Index of nearest x value in x_array.
        
    """

    return np.fabs(x_array - x).argmin()


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