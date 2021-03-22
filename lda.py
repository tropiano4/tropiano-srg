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
        if distribution_type in ['pn', 'p']:
            
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
                    function_array[j] = func_q(q, kFp, kFn)
                    
                expectation_values[i] = 4*np.pi*dr * \
                                        np.sum( r2_array * function_array )
                                        
        elif distribution_type == 'n':
            
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
                    function_array[j] = func_q(q, kFn, kFp)
                    
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
        
    plt.xlim( [0.0, 20.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('r [fm]')
    plt.ylabel(r'$k_F(r)$' + ' [fm' + r'$^{-1}$' + ']')
    plt.show()