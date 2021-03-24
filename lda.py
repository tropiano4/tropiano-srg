#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: lda.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     February 9, 2021
# 
# DESCRIPTION.
#
# Revision history:
#   03/18/21 --- Added 12C data to Densities. Now shows \rho_proton(r) for
#                12C in plots below.
#
#------------------------------------------------------------------------------


from os import getcwd, chdir
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
import observables as ob
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


def hankel_transformation(channel, k_array, k_weights, r_array):
    """
    <r|klm> matrix for given partial wave channel. If len(r_array) = m and len(k_array) = n, then this function 
    returns an m x n matrix.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0').
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
    r_array : 1-D ndarray
        Coordinates array [fm].
        
    Returns
    -------
    M : 2-D ndarray
        Hankel transformation matrix [fm^-3].
        
    Notes
    -----
    The L > 0 transformations may require factors of i or -1. Check conventions.

    """
        
    # L = 0 (0th spherical Bessel function)
    if channel[1] == 'S':
        L = 0
    # L = 1
    elif channel[1] == 'P':
        L = 1
    # L = 2
    elif channel[1] == 'D':
        L = 2
        
    # r_array row vectors and k_array column vectors where both grids are
    # n x m matrices
    k_cols, r_rows = np.meshgrid(k_array, r_array)
    k_weights_cols, _ = np.meshgrid(k_weights, r_array)
        
    M = 2/np.pi * k_cols**2 * k_weights_cols * spherical_jn(L, k_cols * r_rows)

    return M


# Perhaps this function should go in a script within Densities/HFBRAD_SLY4
def load_density(nucleus, nucleon, Z, N):
    """
    Loads a nucleonic density for the given nucleus.
    
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
        
    Returns
    -------
    r_array : 1-D ndarray
        Coordinates array [fm].
    rho_array : 1-D ndarray
        Nucleonic density as a function of r [# of nucleons / vol].
    
    """

    cwd = getcwd()

    # Go to directory corresponding to specified nucleus
    densities_directory = 'Densities/HFBRAD_SLY4/%s' % nucleus
    chdir(densities_directory)
    
    # Load .dens file
    file_name = '%s_%d_%d.dens' % (nucleon, N, Z)
    table = np.loadtxt(file_name)

    chdir(cwd) # Go back to current working directory
    
    r_array = table[:, 0]
    rho_array = table[:, 1]
    
    return r_array, rho_array


class LDA(object):
    
    
    def __init__(self, r_array, rho_p_array, rho_n_array):
        """
        Description.
        
        Parameters
        ----------
        r_array : 1-D ndarray
            Coordinates array [fm].
        rho_p_array : 1-D ndarray
            Proton density as a function of r [# of nucleons / vol].
        rho_n_array : 1-D ndarray
            Neutron density as a function of r [# of nucleons / vol].

        """
        
        # Save the nucleonic densities and coordinates
        self.r_array = r_array
        self.rho_p_array = rho_p_array
        self.rho_n_array = rho_n_array
        

    # BUG TESTING
    def local_density_approximation(self, q_array, func_q, distribution_type):
    #def local_density_approximation(self, q_array):
        """
        Evaluates nuclear-averaged expectation value of the function at a
        range of q values. Function depends on Fermi momentum for protons and
        neutrons.
    
        Parameters
        ----------
        q : 1-D ndarray
            High momentum values [fm^-1].
        func_q : func
            Function that depends on q and, possibly more than 1, k_F.

        Returns
        -------
        expectation_values : 1-D ndarray
            Array of expectation values of the the function evaluated at each
            momentum q.
            
        """
        
        M = len(q_array)
    
        # Load coordinates and nucleonic density
        r_array = self.r_array
        # Number of r_array points
        N = len(r_array)
        
        r2_array = r_array**2
        dr = 0.1 # Spacing between r-points
        # denominator = 4*np.pi * np.sum( dr * r2_array * rho_array )
    
        # \rho(r) = ( 2*k_F(r)^3 ) / ( 3 \pi^2 ) for nucleons (g=4)
        # Evaluate k_F at each point in r_array
        # k_F_array = ( 3*np.pi**2/2 * rho_array )**(1/3)
        
        # pn pair: Two k_F values in this case
        if distribution_type in ['pn', 'p', 'n']:
            
            rho_p_array = self.rho_p_array
            rho_n_array = self.rho_n_array
            
            kFp_array = ( 3*np.pi**2 * rho_p_array )**(1/3)
            kFn_array = ( 3*np.pi**2 * rho_n_array )**(1/3)
        
            expectation_values = np.zeros(M)
            for i, q in enumerate(q_array):
    
                # Now evaluate f(q, kFp, kFn) at each point in q_array and 
                # kF_proton, kF_neutron
                function_array = np.zeros(N)
                for j in range(N):

                    kFp = kFp_array[j]
                    kFn = kFn_array[j]
                    if distribution_type == 'n':
                        function_array[j] = func_q(q, kFn, kFp)
                    else:
                        function_array[j] = func_q(q, kFp, kFn)
                    
                expectation_values[i] = 4*np.pi*dr * \
                                        np.sum( r2_array * function_array )

        else:
            
            if distribution_type == 'pp':
                rho_array = self.rho_p_array
            elif distribution_type == 'nn':
                rho_array = self.rho_n_array
                
            kF_array = ( 3*np.pi**2 * rho_array )**(1/3)

            expectation_values = np.zeros(M)
            for i, q in enumerate(q_array):
    
                # Now evaluate f(q, kF) at each point in q_array and kF
                function_array = np.zeros(N)
                for j, kF in enumerate(kF_array):

                    function_array[j] = func_q(q, kF)
                    
                expectation_values[i] = 4*np.pi*dr * \
                                        np.sum( r2_array * function_array )
  
        return expectation_values
    
    
if __name__ == '__main__':
    
    
    # --- Test densities --- #

    # Details of example nuclei (format is [nuclei, Z, N])
    nuclei_details = [ ['C12', 6, 6], ['O16', 8, 8], ['Ca40', 20, 20],
                       ['Ca48', 20, 28], ['Pb208', 82, 126] ]
    
    # Plot densities as a function of r
    plt.clf()
    for nuclei_list in nuclei_details:
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        r_array, rho_array = load_density(nucleus, nucleon, Z, N)

        plt.plot(r_array, rho_array, label=nucleus)
        
        print(4*np.pi*np.sum(0.1*r_array**2*rho_array))
        
    # --- Show deuteron too --- #
    
    kvnn = 6
    channel = '3S1'
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120
    
    # Load evolved deuteron wave function in momentum-space
    k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
    H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot,
                                      method='srg', generator='Wegner', lamb=lamb)
    U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
    psi_k_unitless = ob.wave_function(H_initial, U=U_matrix)
        
    # Divide out momenta/weights
    factor_array = np.concatenate( (np.sqrt(k_weights) * k_array,
                                    np.sqrt(k_weights) * k_array) ) * np.sqrt(2/np.pi)
    psi_k = psi_k_unitless / factor_array
    
    hank_trans_3S1 = hankel_transformation('3S1', k_array, k_weights, r_array)
    hank_trans_3D1 = hankel_transformation('3D1', k_array, k_weights, r_array)
    
    # sign = -1
    psi_r_3S1 = hank_trans_3S1 @ psi_k[:ntot]
    psi_r_3D1 = hank_trans_3D1 @ psi_k[ntot:]
    
    rho_d = psi_r_3S1**2 + psi_r_3D1**2
    print(4*np.pi*np.sum(0.1*r_array**2*rho_d))
    
    plt.plot(r_array, rho_d, label='Deuteron')
        
    plt.xlim( [0.0, 20.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('r [fm]')
    plt.ylabel(r'$\rho_p(r)$' + ' [fm' + r'$^{-3}$' + ']')
    plt.show()
    
    # Plot k_F for some nuclei as function of r
    # rho = 2 / (3*\pi^2) k_F^3
    plt.clf()
    for nuclei_list in nuclei_details:
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        r_array, rho_array = load_density(nucleus, nucleon, Z, N)
        # kF_array = ( 3*np.pi**2/2 * rho_array )**(1/3)
        kF_array = ( 3*np.pi**2 * rho_array )**(1/3)
    
        plt.plot(r_array, kF_array, label=nucleus)
        
    # Show deuteron too
    plt.plot(r_array, ( 3*np.pi**2 * rho_d )**(1/3), label='Deuteron')
        
    plt.xlim( [0.0, 20.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('r [fm]')
    plt.ylabel(r'$k_F(r)$' + ' [fm' + r'$^{-1}$' + ']')
    plt.show()