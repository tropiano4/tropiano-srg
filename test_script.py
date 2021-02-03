#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a seperate file with an
# extension _testv#.py where v# corresponds to the version number. For example,
# momentum_projection_operator_testv1.py. Use the revision history below to
# document when and why these files are created.
#
# Revision history:
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py in Old_codes.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#   01/26/21 --- Renamed to test_script.py.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Calculate n_{\lambda}(q, Q=0) for different nuclei using LDA. Split into
#   pp, pn, nn, np contributions as in Ryckebusch_2019oya.


from os import getcwd, chdir
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


#### TO DO ####
# 1. Add in low-momentum part
# 2. 1S0 (pn) 3S1-3S1 (pn) and 3S1-3D1 (pn) 1S0 (nn) contributions next
# 3. Make separate folder for tables


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
    
    def __init__(self, pair, kvnn, lamb, r_array, rho_array, 
                 create_table=False):
        """
        Loads the initial and evolved Hamiltonians with the momentum mesh.

        Parameters
        ----------
        pair : str
            Type of pair: 'pp', 'pn', 'nn', or 'np'.
        kvnn : int
            This number specifies the potential.
        lamb : float
            Evolution parameter lambda [fm^-1].
        r_array : 1-D ndarray
            Coordinates array [fm].
        rho_array : 1-D ndarray
            Nucleonic density as a function of r [# of nucleons / vol].
        create_table : bool, opt
            Option to create table for n_lambda_pair as a function of q and
            k_F.
            
        """
        
        # Copying this from the function before
        
        self.pair = pair
        self.kvnn = kvnn
        self.r_array = r_array
        self.rho_array = rho_array
        
        if pair == 'pp' or pair == 'nn':
            
            channel = '1S0'
            
            # Load initial and evolved Hamiltonians with momentum mesh
            # This function will be called in an integration routine so you really
            # don't want to load the momentum mesh here. Ideally it'd reference
            # a variable within a class self.k_array and self.k_weights
            k_array, k_weights = vnn.load_momentum(kvnn, channel)
            self.k_array = k_array
            self.k_weights = k_weights
            # Length of k_array
            ntot = len(k_array)
            self.ntot = ntot
            factor_array = k_array * np.sqrt( 2*k_weights / np.pi )
            self.factor_array = factor_array # [fm^-3/2]
            row, col = np.meshgrid(factor_array, factor_array)
            self.row = row
            self.col = col
        
            # Units are MeV
            H_initial = vnn.load_hamiltonian(kvnn, channel) # MeV
            # Use band-diagonal evolution as default
            H_evolved = vnn.load_hamiltonian(kvnn, channel, method='srg', 
                                             generator='Wegner', lamb=lamb) 
        
            # Calculate SRG transformation
            U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
        
            # Subtract out identity matrix and divide out momentum and weights
            I = np.eye(ntot, ntot)
            delta_U_matrix = (U_matrix - I) / row / col # [fm^3]
            
            self.delta_U_matrix = delta_U_matrix
            
        if create_table:
            
            self.create_n_lambda_pair_table()
            
        # Load n_lambda_pair_function
        self.load_n_lambda_pair_table()
        
    
    def n_lambda_1S0(self, q, k_F):
        """
        Gives the N,Np contribution to the SRG-evolved pair momentum
        distribution under the following approximations:
            1. Q_vector = 0.
            2. q >> k_F.
            3. Taking only S-waves: 1S0, 3S1-3S1, 3S1-3D1.
    
        Parameters
        ----------
        q : float
            High momentum value [fm^-1].
        k_F : float
            Fermi momentum [fm^-1].

        Returns
        -------
        pair_contribution : float
            Pair momentum distribution for given pair at q and Q=0 [fm^3].
    
        Notes
        -----
        1. Make functions for each factor with arguments T, M_T, etc.
        2. Should lamb be replaced by k_F as an argument instead? Or should 
           lamb vary as k_F does (that would require many more evolution 
           snapshots)?
        3. Should it be (U_matrix - I) / row / col or 
           U_matrix / row / col - I?
        4. Check the normalization of this function. Integrating over d3q and 
           d3Q should give 1, correct?
        5. Should factor of 1/8 be 1/4 in term four?
    
        """
        
        # T = 1
        # M_T = 1 or M_T = -1 (doesn't matter!)
        # Make function that takes T, M_T as argument?
        # Clebsch-Gordan coefficients for isospin are 1 here
        # < 1/2 1/2 | 1 1 > = 1 or < -1/2 -1/2 | 1 -1 > = 1
        isospin_factor = 1
        
        # Select the allowed value of S and M_S
        # S = 0
        # M_S = 0
        # Make function that takes S, M_S as argument?
        
        # We sum over m_s, m_s', m_s_1, m_s_2 restricted to the following
        # values:
        #   m_s = up, m_s' = down, m_s_1 = up, m_s_2 = down
        #   m_s = up, m_s' = down, m_s_1 = down, m_s_2 = up
        #   m_s = down, m_s' = up, m_s_1 = up, m_s_2 = down
        #   m_s = down, m_s' = up, m_s_1 = down, m_s_2 = up
        # Each time we get factors of 1/\sqrt(2) or -1/\sqrt(2) but
        # appearing in even powers - thus, it's always a factor of 1/4.
        # But there are 4 possible combinations and we are summing over them!
        spin_factor = 1
        
        # Select the associated partial wave channel
        # channel = '1S0'
        # This means L=0, M_L=0, J=0, M_J=0
        # < M_L M_S | J M_J> factors are always 1 in this case
        total_ang_momentum_factor = 1
        # Is this ever not 1?
        
        # Expansion of k_vector and q_vector factors
        k_expansion_factor = 1 / (2*np.pi**2) # ( \sqrt(2/\pi)*Y00 )^2
        q_expansion_factor = 1 / (2*np.pi**2) # ( \sqrt(2/\pi)*Y00 )^2
        
        # Load delta_U and momentum mesh
        k_array = self.k_array
        delta_U_matrix = self.delta_U_matrix
        
        # Index of q in k_array
        q_index = find_q_index(q, k_array)
        
        # Do first three terms if q < k_F
        if q < k_F:
            
            # First term in Overleaf notes
            # Sum is over m_t, m_t', then m_s, m_s'
            # Factor of 2 for m_s and m_s' combinations
            factor_1 = 2
            I = np.eye(self.ntot, self.ntot) / self.row / self.col # [fm^3]
            # pair_contribution_1 = 1/2 * I[q_index, q_index] * factor_1
            pair_contribution_1 = 0
            
            # Second and third term in Overleaf notes
            factor_2 = isospin_factor * spin_factor * \
                       total_ang_momentum_factor * q_expansion_factor
                       
            delta_U_qq = delta_U_matrix[q_index, q_index] # [fm^3]
            delta_U_dag_qq = delta_U_matrix.T[q_index, q_index] # [fm^3]
            
            # Minus sign???
            # pair_contribution_2 = -1/4 * (delta_U_qq + delta_U_dag_qq) * \
            #                       factor_2 # [fm^3]
            pair_contribution_2 = 0
                                  
        else:
            # \theta functions give 0 here
            pair_contribution_1 = 0
            pair_contribution_2 = 0
    
        # Fourth term in Overleaf notes
        # Let's combine all the relevant factors here to be concise
        factor_3 = isospin_factor * spin_factor * total_ang_momentum_factor * \
                   k_expansion_factor * q_expansion_factor
                 
        # Find the matrix element of delta_U(k, q)
        delta_U_vector = delta_U_matrix[:self.ntot, q_index]
        # Do the same for delta_U^{\dagger}
        delta_U_dag_vector = delta_U_matrix.T[q_index, :self.ntot]
        
        # The overall contribution to the evolved pair momentum distribution
        # as a function of k (which is summed to k_F later) [fm^6]
        pair_contribution_k = delta_U_vector * delta_U_dag_vector * factor_3
        # 1/(2*\pi^2)^2 * \delta U_1S0(k,q) * \delta U_1S0^{\dagger}(q,k)
                              
        # Sum over k up to k_F (this should probably be an intergral?)
        # \sum_k k_i^2 w_i -> \int dk k^2
        # I have \sum_k... But the weights and k factors should have come from
        # U O U^t (two integrals over k and k')
        # One of the k' is done by \delta function
        # There is one more. So yes, do 2 / \pi \sum_k k_i^2 w_i
        integrand = self.factor_array * pair_contribution_k # [fm^3]
        
        # Set k_F cutoff in integration limit where we prevent undershooting
        # k_F (meaning \Lambda will never be less than k_F)
        cutoff_index = find_q_index(k_F, k_array)
        while k_array[cutoff_index] < k_F:
            cutoff_index += 1
        integrand_resized = integrand[:cutoff_index+1]

        pair_contribution_3 = 1/8 * np.sum(integrand_resized) # [fm^3]
        # pair_contribution_3 = 0
        
        # Sum all the terms for the total contribution
        pair_contribution = pair_contribution_1 + pair_contribution_2 + \
                            pair_contribution_3

        return pair_contribution
    
    
    def create_n_lambda_pair_table(self):
        """
        Create a table for interpolation over q and k_F.
        
        """
        
        # Set q_array and k_F_array values for table
        # Do momentum values from 1-6 fm^-1 and k_F values from 0-1.5 fm^-1
        k_array = self.k_array
        # k_min_index = find_q_index(0.9, k_array)
        # k_max_index = find_q_index(6.1, k_array)
        # q_array = k_array[k_min_index:k_max_index+1]
        q_array = k_array
        k_F_array = np.linspace(0.0, 1.5, 30)
        
        # Write table in txt file
        f = open( 'n_lambda_%s_table_kvnn_%d.txt' % (self.pair, self.kvnn), 
                                                     'w' )
        g = open('q_table_kvnn_%d.txt' % self.kvnn, 'w')
        h = open('k_F_table_kvnn_%d.txt' % self.kvnn, 'w')
        
        for i, q in enumerate(q_array):
            for j, k_F in enumerate(k_F_array):
                f.write( '{:<16.5e}'.format( self.n_lambda_1S0(q, k_F) ) )
                if i == 0:
                    h.write('%.5f\n'%k_F)
            g.write('%.5f\n'%q)
            f.write('\n')
            
        f.close()
        g.close()
        h.close()
                
        
    def load_n_lambda_pair_table(self):
        """
        Set the interpolated version of n_lambda_pair as a function of q and
        k_F.

        """

        n_lambda_pair_table = np.loadtxt( 'n_lambda_%s_table_kvnn_%d.txt' % \
                                         (self.pair, self.kvnn) )
        q_array = np.loadtxt( 'q_table_kvnn_%d.txt' % self.kvnn )
        k_F_array = np.loadtxt( 'k_F_table_kvnn_%d.txt' % self.kvnn )
        
        n_lambda_pair_func = RectBivariateSpline(q_array, k_F_array, 
                                                 n_lambda_pair_table)
        
        # This is a function of q and k_F: n(q, k_F)
        self.n_lambda_pair_func = n_lambda_pair_func
        

    def local_density_approximation(self, q_array):
        """
        Evaluates nuclear-averaged expectation value of an SRG-evolved pair 
        momentum distribution at a range of q values.
    
        Parameters
        ----------
        q : 1-D ndarray
            High momentum values [fm^-1].

        Returns
        -------
        n_lambda_pair_exp_array : 1-D ndarray
            Array of expectation values of the SRG-evolved pair momentum 
            distribution evaluated at each relative momentum q with total 
            momentum Q=0.
            
        """
        
        M = len(q_array)
    
        r_array = self.r_array
        rho_array = self.rho_array
        # Number of r_array points
        N = len(r_array)
        
        r2_array = r_array**2
        dr = 0.1 # Spacing between r-points
        denominator = 4*np.pi * np.sum( dr * r2_array * rho_array )
    
        # \rho(r) = ( 2*k_F(r)^3 ) / ( 3 \pi^2 ) for nucleons (g=4)
        # Evaluate k_F at each point in r_array
        k_F_array = ( 3*np.pi**2/2 * rho_array )**(1/3)
        
        n_lambda_pair_exp_array = np.zeros(M)
        for i, q in enumerate(q_array):
    
            # Now evaluate n_lambda_pair at each point in q_array and k_F_array
            n_lambda_pair = np.zeros(N)
            for j, k_F in enumerate(k_F_array):
                n_lambda_pair[j] = self.n_lambda_pair_func(q, k_F)
        
            # Integrate over r
            integrand = dr * r2_array * n_lambda_pair * rho_array
            numerator = 4*np.pi * np.sum(integrand)
            # Factor of 4*\pi cancels out
            n_lambda_pair_exp_array[i] = numerator / denominator
            
        return n_lambda_pair_exp_array


if __name__ == '__main__':
    
    pair = 'pp'
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    
    # --- Test densities --- #
    
    # Confused by the N and Z for Pb208?
    
    # Details of example nuclei (format is [nuclei, Z, N])
    nuclei_details = [ ['O16', 8, 8], ['Ca40', 20, 20], ['Ca48', 20, 28] ]
    
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
        kF_array = ( 3*np.pi**2/2 * rho_array )**(1/3)
    
        plt.plot(r_array, kF_array, label=nucleus)
        
    plt.xlim( [0.0, 20.0] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('r [fm]')
    plt.ylabel(r'$k_F(r)$' + ' [fm' + r'$^{-1}$' + ']')
    plt.show()
    
    
    # --- Write the relevant tables --- #
    # Comment this out when you're done
    
    for nuclei_list in nuclei_details:
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        r_array, rho_array = load_density(nucleus, nucleon, Z, N)
        lda = LDA(pair, kvnn, lamb, r_array, rho_array, create_table=True)


    # --- Plot n_lambda_pair_exp for pp pairs --- #
    
    q_array = np.linspace(0.3, 6.0, 100)
    
    plt.clf()
    for nuclei_list in nuclei_details:
        
        # Plot density for some nuclei here
        nucleus = nuclei_list[0]
        nucleon = 'proton'
        Z = nuclei_list[1]
        N = nuclei_list[2]
        r_array, rho_array = load_density(nucleus, nucleon, Z, N)
        lda = LDA(pair, kvnn, lamb, r_array, rho_array)

        n_lambda_pair_exp_array = lda.local_density_approximation(q_array)
    
        # plt.plot(q_array, n_lambda_pair_exp_array, label=nucleus)
        plt.semilogy(q_array, n_lambda_pair_exp_array, label=nucleus)
        
    plt.xlim( [min(q_array), max(q_array)] )
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('q [fm' + r'$^{-1}$' + ']')
    plt.ylabel(r'$<n_{\lambda}^{pp}(q,Q=0)>$' + ' [fm' + r'$^3$' + ']')
    plt.title('kvnn = %d, ' % kvnn + r'$\lambda=%.2f$' % lamb + 
              ' [fm' + r'$^{-1}$' + ']')
    plt.show()