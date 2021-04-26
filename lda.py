#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: lda.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     February 9, 2021
# 
# Calculates nuclear-averaged momentum distributions under a local density
# approximation. Relies on nucleonic densities from Densities/HFBRAD_SLY4.
#
# Revision history:
#   03/18/21 --- Added 12C data to Densities. Now shows \rho_proton(r) for
#                12C in plots below.
#   04/16/21 --- Added 56Fe data to Densities. Now shows \rho_proton(r) for
#                56Fe in plots below.
#   04/23/21 --- Added method to average over multiple contributions returned
#                from input function. For example, we can now input a function
#                that returns the total, pp, and pn contributions, then average
#                each contribution separately.
#
#------------------------------------------------------------------------------


from os import getcwd, chdir
import numpy as np


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
            Proton density as a function of r [fm^-3].
        rho_n_array : 1-D ndarray
            Neutron density as a function of r [fm^-3].

        """
        
        # Save the r values and nucleonic densities
        self.r_array = r_array
        self.rho_p_array = rho_p_array
        self.rho_n_array = rho_n_array
        

    def local_density_approximation(self, q_array, func_q, distribution_type,
                                    contributions='total'):
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
        distribution_type : str
            Type of momentum distribution (e.g., 'pn'). This determines
            whether the function takes one or two k_F inputs.
        contributions : str, optional
            Option to return different contributions to the momentum
            distribution.
            1. Default is 'total' which only returns the total momentum
               distribution.
            2. Specify 'NN_contributions' for total, pp, and pn (or nn, np)
               where the nucleon-nucleon contributions are isolated in the
               high-q term.
            3. Specify 'q_contributions' for total, 1, \delta U, and 
               \delta U \delta U^\dagger.

        Returns
        -------
        expectation_values : 2-D ndarray
            Array of expectation values of the function for each contribution
            evaluated at each momentum q.
            
        """
        
        # Set shape of return array with 'axes'
        # Only return total contribution at each q
        if contributions == 'total':
            axes = 1
        # Return total, pp (nn), and pn (np) contributions at each q
        elif contributions == 'NN_contributions':
            axes = 3 
        # Return total, 1, \delta U, and \delta U^2 contributions at each q
        elif contributions == 'q_contributions':
            axes = 4
            
        # Load r_array
        r_array = self.r_array
        # Number of r_array points
        mtot = len(r_array)
        r2_grid, _ = np.meshgrid(r_array**2, np.zeros(axes), indexing='ij')
        dr = 0.1 # Spacing between r-points
            
        # Length of q_array
        ntot = len(q_array)
        
        # Evaluate f(q, kF) at each point in q_array and kF
        expectation_values = np.zeros( (ntot, axes) )

        # pn pair or single-nucleon: Two k_F values in this case
        if distribution_type in ['pn', 'p', 'n']:
            
            rho_p_array = self.rho_p_array
            rho_n_array = self.rho_n_array
            
            # Evaluate k_F at each point in r_array
            kFp_array = (3*np.pi**2 * rho_p_array)**(1/3)
            kFn_array = (3*np.pi**2 * rho_n_array)**(1/3)
            
            # Loop over q
            for i, q in enumerate(q_array):
    
                function_array = np.zeros( (mtot, axes) )
                
                # Loop over r for k_F values
                for j, (kFp, kFn) in enumerate( zip(kFp_array, kFn_array) ):

                    if distribution_type == 'n':
                        function_array[j, :] = func_q(q, kFn, kFp,
                                               contributions=contributions)
                    else:
                        function_array[j, :] = func_q(q, kFp, kFn,
                                               contributions=contributions)

                # Integrate over r for each contribution (summing over axis=0)
                expectation_values[i, :] = 4*np.pi*dr * \
                                           np.sum(r2_grid * function_array,
                                                  axis=0)
        
        # pp or nn pair: One k_F value
        else:
            
            if distribution_type == 'pp':
                rho_array = self.rho_p_array
            elif distribution_type == 'nn':
                rho_array = self.rho_n_array
            
            # Evaluate k_F at each point in r_array
            kF_array = (3*np.pi**2 * rho_array)**(1/3)

            # Loop over q
            for i, q in enumerate(q_array):
    
                function_array = np.zeros( (mtot, axes) )
                
                # Loop over r for k_F values
                for j, kF in enumerate(kF_array):

                    function_array[j, :] = func_q(q, kF,
                                                  contributions=contributions)
                    
                expectation_values[i, :] = 4*np.pi*dr * \
                                           np.sum(r2_grid * function_array,
                                                  axis=0)
  
        return expectation_values
    
    
if __name__ == '__main__':
    
    
    # --- Test densities --- #
    
    import matplotlib.pyplot as plt

    # Details of example nuclei (format is (nuclei, Z, N) )
    nuclei_details = ( ('C12', 6, 6), ('O16', 8, 8), ('Ca40', 20, 20),
                       ('Ca48', 20, 28), ('Fe56', 26, 30), ('Pb208', 82, 126) )
    
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
        
        print( 4*np.pi*np.sum(0.1 * r_array**2 * rho_array) )
        
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
        kF_array = ( 3*np.pi**2 * rho_array )**(1/3)
    
        plt.plot(r_array, kF_array, label=nucleus)

    plt.xlim( [0.0, 20.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('r [fm]')
    plt.ylabel(r'$k_F(r)$' + ' [fm' + r'$^{-1}$' + ']')
    plt.show()