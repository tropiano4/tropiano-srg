#!/usr/bin/env python3

"""
File: densities.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 27, 2021

Handles nucleonic densities from data sub-directories.
So far we have densities from the SLy4 Skyrme functional using HFBRAD code,
densities from the Gogny functional using HFBTHO code, and densities from VMC
calculations (see www.phy.anl.gov/theory/research/density/).

Last update: May 3, 2022

"""

# Python imports
import numpy as np


def load_density(nucleon, nucleus_name, Z=None, N=None, density='Gogny'):
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
    Z : int, optional
        Proton number of the nucleus. This argument is required for SLy4.
    N : int, optional
        Neutron number of the nucleus. This argument is required for SLy4.
    density : str, optional
        Name of density (e.g., 'SLy4').
        
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
    errors. (This only happens for VMC densities.)
    
    """


    # Go to directory corresponding to the specified nucleus and density
    if density == 'SLy4':
        
        densities_directory = f'../data/dft/{density}/{nucleus_name}/'
        file_name = f'{nucleon}_{N:d}_{Z:d}.dens'
        column_number = 1
        
    elif density == 'Gogny':
        
        densities_directory = f'../data/dft/{density}/{nucleus_name}/'
        file_name = 'DensityQP.dat'
        if nucleon == 'proton':
            column_number = 1
        elif nucleon == 'neutron':
            column_number = 2
    
    elif density == 'VMC':
        
        densities_directory = '../data/vmc/densities/'
        file_name = f'{nucleus_name}_single_nucleon.txt'
        
        # VMC files either have single \rho column for N=Z nuclei or
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