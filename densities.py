#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: densities.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 27, 2021
# 
# Loads nucleonic densities from Densities directory. So far relying on only
# Skyrme EDFs from SLy4 using HFBRAD code.
#
# Revision history:
#   03/18/21 --- Added 12C data to Densities. Now shows \rho_proton(r) for
#                12C in plots below.
#   04/16/21 --- Added 56Fe data to Densities. Now shows \rho_proton(r) for
#                56Fe in plots below.
#   05/27/21 --- Split up old lda.py code into this piece and added the
#                averaging over \int dr r^2 to snmd.py and pmd.py separately.
#                Renamed from lda.py to densities.py where lda.py is now in
#                Old_codes.
#   06/17/21 --- Added He4 from www.phy.anl.gov/theory/research/density/.
#   06/30/21 --- Added He8 from www.phy.anl.gov/theory/research/density/.
#
#------------------------------------------------------------------------------


from os import getcwd, chdir
import numpy as np


def load_density(nucleus, nucleon, Z, N, edf='SLY4'):
    """
    Loads a nucleonic density for the given nucleus. Densities are normalized
    according to 4*\pi \int_0^\infty dR R^2 \rho_A(R) = Z or N.
    
    Parameters
    ----------
    nucleus : str
        Specify nucleus (e.g., 'O16', 'Ca40', etc.)
    nucleon : str
        Specify 'proton' or 'neutron'.
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
    
    """

    cwd = getcwd()

    # Go to directory corresponding to specified nucleus
    if edf == 'SLY4':
        
        densities_directory = 'Densities/HFBRAD_%s/%s' % (edf, nucleus)
        file_name = '%s_%d_%d.dens' % (nucleon, N, Z)
        column_number = 1
        
    # Densities from www.phy.anl.gov/theory/research/density/
    elif edf == 'AV18':
        
        densities_directory = 'Densities/%s/%s' % (edf, nucleus)
        file_name = 'densities_%d_%d.txt' % (N, Z)
        
        # AV18 files either have single \rho column for N=Z nuclei or
        # two columns for proton (1) and neutron (3)
        if N == Z:
            column_number = 1
        else:
            if nucleon == 'proton':
                column_number = 1
            elif nucleon == 'neutron':
                column_number = 3 
        
    chdir(densities_directory)
    
    # Load file
    table = np.loadtxt(file_name)

    # Go back to current working directory
    chdir(cwd)
    
    R_array = table[:, 0]
    rho_array = table[:, column_number]
    
    return R_array, rho_array


if __name__ == '__main__':
    
    
    # --- Test densities --- #
    
    import matplotlib.pyplot as plt

    # Details of example nuclei (format is (nuclei, Z, N) )
    nuclei_details = ( ('He4', 2, 2), ('C12', 6, 6), ('O16', 8, 8),
                       ('Ca40', 20, 20), ('Ca48', 20, 28), ('Fe56', 26, 30),
                       ('Pb208', 82, 126) )
    
    # Plot densities as a function of R
    plt.clf()
    for i, nuclei_list in enumerate(nuclei_details):
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        if nucleus == 'He4':
            R_array, rho_array = load_density(nucleus, nucleon, Z, N, 'AV18')
        else:
            R_array, rho_array = load_density(nucleus, nucleon, Z, N)

        plt.plot(R_array, rho_array, label=nucleus)
        
        print( 4*np.pi*np.sum(0.1 * R_array**2 * rho_array) )
        
    plt.xlim( [0.0, 10.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('R [fm]')
    plt.ylabel(r'$\rho_p(R)$' + ' [fm' + r'$^{-3}$' + ']')
    plt.show()
    
    # Plot proton k_F for some nuclei as function of R
    # rho = 2 / (3*\pi^2) k_F^3
    plt.clf()
    for j, nuclei_list in enumerate(nuclei_details):
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        if nucleus == 'He4':
            R_array, rho_array = load_density(nucleus, nucleon, Z, N, 'AV18')
        else:
            R_array, rho_array = load_density(nucleus, nucleon, Z, N)
        kF_array = ( 3*np.pi**2 * rho_array )**(1/3)

        plt.plot(R_array, kF_array, label=nucleus)

    plt.xlim( [0.0, 15.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('R [fm]')
    plt.ylabel(r'$k_F(R)$' + ' [fm' + r'$^{-1}$' + ']')
    plt.show()