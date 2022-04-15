#!/usr/bin/env python3

"""
File: densities.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 27, 2021

Handles nucleonic densities from the directory data/densities.
So far we have densities from the SLy4 Skyrme functional using HFBRAD code,
densities from the Gogny functional, and densities from VMC calculations
(see www.phy.anl.gov/theory/research/density/).

Last update: April 9, 2022

"""

# To-do: Could always rename this script (e.g., nuclei.py).
# To-do: Probably want this to be a class too.
# To-do: Check zero division problem again.
# To-do: How to do AV18 densities (properly). See comment in function.

# Python imports
import numpy as np

# Imports from A.T. codes
from modules.labels import replace_periods


def load_density(nucleon, nucleus_name, Z, N, edf='SLY4'):
    """
    Loads a nucleonic density for the given nucleus. Densities are normalized
    according to
        4*\pi \int_0^\infty dR R^2 \rho_A(R) = Z or N.
    
    Parameters
    ----------
    nucleon : str
        Specify 'proton' or 'neutron'.
    nucleus_name : str
        Specify the nucleus (e.g., 'O16', 'Ca40', etc.)
    Z : int
        Proton number of the nucleus.
    N : int
        Neutron number of the nucleus.
    edf : str, optional
        Name of EDF (e.g., 'SLY4').
        
    Returns
    -------
    R_array : 1-D ndarray
        C.o.M. coordinates [fm].
    rho_array : 1-D ndarray
        Nucleonic density as a function of R [# of nucleons / vol].
                                              
    Notes
    -----
    Momentum distributions code compute intermediate integration arrays in 
    relative k and C.o.M. K which rely on kF(R) values. These values can be
    zero if the density \rho(R) = 0. We must replace zeros in \rho(R) with an
    extremely small number so the codes run correctly to avoid zero division
    errors. (This only happens for edf = 'AV18' densities.)
    
    """


    # Go to directory corresponding to specified nucleus and EDF
    if edf == 'SLY4':
        
        densities_directory = f'../densities/HFBRAD_{edf}/{nucleus_name}/'
        file_name = f'{nucleon}_{N:d}_{Z:d}.dens'
        column_number = 1
        
    elif edf == 'Gogny':
        
        densities_directory = f'../densities/{edf}/{nucleus_name}/'
        file_name = 'DensityQP.dat'
        if nucleon == 'proton':
            column_number = 1
        elif nucleon == 'neutron':
            column_number = 2
    
    # Technically it doesn't make sense to have a case edf == AV18 since
    # AV18 does not use an EDF. It would also make more sense to call it VMC.
    elif edf == 'AV18':
        
        densities_directory = '../densities/{edf}/'
        file_name = '{nucleus_name}_densities_{N:d}_{Z:d}.txt'
        
        # AV18 files either have single \rho column for N=Z nuclei or
        # two columns for proton (1) and neutron (3)
        if N == Z:
            column_number = 1
        else:
            if nucleon == 'proton':
                column_number = 1
            elif nucleon == 'neutron':
                column_number = 3 
        
    # Load file
    table = np.loadtxt(densities_directory + file_name)
    
    R_array = table[:, 0]
    rho_array = table[:, column_number]
    
    # Avoiding zero division errors
    zero_case = rho_array == 0
    rho_array[zero_case] = 1e-30
    
    return R_array, rho_array


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
    l_str : str
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
        print('Input l value is outside the range of this function.')
        return None


def sp_states(nucleus, print_statement=False):
    """
    Return all the occupied single-particle states of the given nucleus.
    
    Parameters
    ----------
    nucleus : tuple
        Details for various nuclei formatted as a tuple:
        (name (str), Z (int), N (int)) (e.g., ('O16', 8, 8)).
    print_statement : bool, optional
        Option to print information for each s.p. state in nucleus.
    
    Returns
    -------
    output : list
        List of lists where the first (second) corresponds to all the occupied
        neutron (proton) s.p. states, which are strings (e.g., '1s0p5' means
        1s with j=1/2).
    
    Notes
    -----
    Currently we're assuming the SLy4 interaction.
    
    """
    
    # Get nucleus name, proton number, and neutron number
    nucleus_name = nucleus[0]
    Z = nucleus[1]
    N = nucleus[2]
    
    # Go to HFBRAD directory
    densities_directory = f'../densities/HFBRAD_SLY4/{nucleus_name}/'
    file_name = f'hfb_{N}_{Z}.spe'
    
    # Open file and add each occupied s.p. state to list
    neutron_states = []
    proton_states = []
    
    f = open(densities_directory + file_name, 'r')
    
    for line in f:
        
        unit = line.strip().split() # Split up row into list
        
        # Make sure you're going through the correct data
        if (len(unit) == 12) and (unit[0] in ('1', '2')):
            
            # Only do occupied states:
            if float(unit[6]) == 1:
            
                # Integer specifying neutron or proton
                nucleon_number = unit[0]
        
                # j value
                j = int(unit[1])/2
        
                # Orbital angular momentum (int)
                l = int(unit[2])
                # Orbital angular momentum (str)
                l_str = convert_l_to_string(l)
        
                # Is this correct? (# of nodes)
                n = unit[11]
            
                # Convert s.p. state to string and append to list
                state_str = f'{n}{l_str}{j:.1f}'
            
                # Add string to neutron or proton list with periods replaced
                # by 'p'
                if nucleon_number == '1': # Neutron
                    neutron_states.append(replace_periods(state_str))
                elif nucleon_number == '2': # Proton
                    proton_states.append(replace_periods(state_str))
            
                # Print information for each state?
                if print_statement:
                
                    info = 'Nuc={:s}, N={:s}, state={:s}'.format(
                        nucleon_number, unit[4], state_str
                    )
                    print(info)
                
    f.close()
    
    return [neutron_states, proton_states]


if __name__ == '__main__':
    
    
    # --- Test densities --- #
    
    # Python imports
    import matplotlib.pyplot as plt
    
    # Imports from A.T. codes
    from modules import figure_graphics as fg

    # Details of example nuclei (format is (nuclei, Z, N) )
    nuclei_details = (
        ('He4', 2, 2), ('C12', 6, 6), ('O16', 8, 8), ('Ca40', 20, 20),
        ('Ca48', 20, 28), ('Fe56', 26, 30), ('Pb208', 82, 126)
    )
    
    # Plot densities as a function of R
    plt.clf()
    for i, nuclei_list in enumerate(nuclei_details):
        
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        
        R_array, rho_array = load_density(nucleus, nucleon, Z, N)
            
        curve_color = fg.xkcd_colors(i)

        plt.plot(R_array, rho_array, label=nucleus, color=curve_color)
        
        normalization = 4*np.pi*np.sum(0.1*R_array**2*rho_array)
        print(normalization)
        
    plt.xlim((0.0, 10.0))
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('R [fm]')
    plt.ylabel(r'$\rho_p(R)$' + ' [fm' + r'$^{-3}$' + ']')
    plt.show()
    
    # Plot proton Fermi momentum for same set of nuclei as function of R
    plt.clf()
    for j, nuclei_list in enumerate(nuclei_details):
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        
        R_array, rho_array = load_density(nucleus, nucleon, Z, N)
        
        kF_array = (3*np.pi**2*rho_array)**(1/3)
        
        curve_color = fg.xkcd_colors(j)

        plt.plot(R_array, kF_array, label=nucleus, color=curve_color)

    plt.xlim((0.0, 15.0))
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('R [fm]')
    plt.ylabel(r'$k_F(R)$' + ' [fm' + r'$^{-1}$' + ']')
    plt.show()