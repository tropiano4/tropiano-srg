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
#   ../../.. --- Minor revisions to r^2 operator.
#
#------------------------------------------------------------------------------


from os import getcwd, chdir
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.


# TO DO
# 1. In hindsight, there's no need to interpolate.


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
    
    
    def __init__(self, function, function_name, kvnn, r_array, rho_array, 
                 create_table=False, q_array=None, k_F_array=None):
        """
        Description.
        
        Parameters
        ----------
        function : func
            Function dependent on q [fm^-1] and k_F [fm^-1].
        function_name : str
            Name of the function for file writing/loading purposes.
        kvnn : int
            This number specifies the potential. (We're assuming the function
            will depend on some NN potential.)
        r_array : 1-D ndarray
            Coordinates array [fm].
        rho_array : 1-D ndarray
            Nucleonic density as a function of r [# of nucleons / vol].
        create_table : bool, opt
            Option to create table for the given function with respect to q
            and k_F. Default is False.
        q_array : 1-D ndarray, opt
            Momentum values [fm^-1]. Default is None (assuming table is
            already created).
        k_F_array : 1-D ndarray, opt
            Fermi momentum values [fm^-1]. Default is None (assuming table is
            already created).

        """
        
        # Save the function
        self.function = function
        
        # Save the name of the function and kvnn
        self.function_name = function_name
        self.kvnn = kvnn
        
        # Save the nucleonic densities and coordinates
        self.r_array = r_array
        self.rho_array = rho_array
        
        # # Create the relevant tables for interpolation
        # if create_table:
            
        #     self.create_table(q_array, k_F_array)
            
        # # Load interpolated function for LDA calculation
        # self.load_table()

    
    def create_table(self, q_array, k_F_array):
        """
        Create a table of the function for interpolation over q and k_F.
    
        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        k_F_array : 1-D ndarray
            Fermi momentum values [fm^-1].
        
        """
        
        # Write table in txt file
        f = open( '%s_kvnn_%d.txt' % (self.function_name, self.kvnn), 'w' )
        g = open( 'q_table_kvnn_%d.txt' % self.kvnn, 'w' )
        h = open( 'k_F_table_kvnn_%d.txt' % self.kvnn, 'w' )
        
        for i, q in enumerate(q_array):
            for j, k_F in enumerate(k_F_array):
                f.write( '{:<16.5e}'.format( self.function(q, k_F) ) )
                if i == 0:
                    h.write('%.5f\n'%k_F)
            g.write('%.5f\n'%q)
            f.write('\n')
            
        f.close()
        g.close()
        h.close()
    

    def load_table(self):
        """
        Set the interpolated version of f(q, k_F).

        """

        function_table = np.loadtxt( '%s_kvnn_%d.txt' % (self.function_name, 
                                                     self.kvnn) )
        q_array = np.loadtxt( 'q_table_kvnn_%d.txt' % self.kvnn )
        k_F_array = np.loadtxt( 'k_F_table_kvnn_%d.txt' % self.kvnn )
        
        function_int = RectBivariateSpline(q_array, k_F_array, function_table)
        
        # This is an interpolated version of the function for various q and 
        # k_F values
        self.function_int = function_int
        

    # BUG TESTING
    def local_density_approximation(self, q_array, func_q_k_F):
    #def local_density_approximation(self, q_array):
        """
        Evaluates nuclear-averaged expectation value of the function at a
        range of q values.
    
        Parameters
        ----------
        q : 1-D ndarray
            High momentum values [fm^-1].

        Returns
        -------
        expectation_values : 1-D ndarray
            Array of expectation values of the the function evaluated at each
            relative momentum q.
            
        """
        
        M = len(q_array)
    
        # Load coordinates and nucleonic density
        r_array = self.r_array
        rho_array = self.rho_array
        # Number of r_array points
        N = len(r_array)
        
        r2_array = r_array**2
        dr = 0.1 # Spacing between r-points
        denominator = 4*np.pi * np.sum( dr * r2_array * rho_array )
    
        # \rho(r) = ( 2*k_F(r)^3 ) / ( 3 \pi^2 ) for nucleons (g=4)
        # Evaluate k_F at each point in r_array
        # k_F_array = ( 3*np.pi**2/2 * rho_array )**(1/3)
        k_F_array = ( 3*np.pi**2 * rho_array )**(1/3)
        
        expectation_values = np.zeros(M)
        for i, q in enumerate(q_array):
    
            # Now evaluate f(q, k_F) at each point in q_array and k_F_array
            function_array = np.zeros(N)
            for j, k_F in enumerate(k_F_array):
                # BUG TESTING
                #function_array[j] = self.function_int(q, k_F)
                function_array[j] = func_q_k_F(q, k_F)
        
            # Integrate over r
            integrand = dr * r2_array * function_array * rho_array
            numerator = 4*np.pi * np.sum(integrand)
            # Factor of 4*\pi cancels out
            expectation_values[i] = numerator / denominator
                
        return expectation_values
    
    
if __name__ == '__main__':
    
    
    # --- Test densities --- #

    # Details of example nuclei (format is [nuclei, Z, N])
    nuclei_details = [ ['O16', 8, 8], ['Ca40', 20, 20], ['Ca48', 20, 28],
                       ['Pb208', 82, 126] ]
    
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