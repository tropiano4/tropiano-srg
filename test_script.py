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
#   03/16/21 --- Testing pair momentum distributions for different nuclei
#                using LDA with simple filled Fermi sea single-particle
#                momentum distributions. Created pmd_simple_test.py in
#                Old_codes.
#   03/16/21 --- Testing pair momentum distribution for deuteron starting from
#                second quantization derivation and using expansion of
#                U(k, k'). Created pmd_deuteron_test.py in Old_codes.
#   03/16/21 --- Testing single-nucleon momentum distributions for different
#                nuclei using LDA. Created single_particle_momentum_dist.py.
#                in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Calculate n_{\lambda}(q) for different nuclei using LDA.


import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
# Scripts made by A.T.
from lda import load_density, LDA
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class pair_momentum_distributions(object):
    
    
    def __init__(self, kvnn, lamb, kmax=0.0, kmid=0.0, ntot=0):
        """
        Description.

        Parameters
        ----------
        kvnn : TYPE
            DESCRIPTION.
        lamb : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        channels = ['1S0', '3S1']
        
        d = {}
        
        for channel in channels:
            
            # Store relevant pieces in dictionary for each channel
            d[channel] = {}
            
            # Load momentum mesh and store in dictionary
            k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid,
                                                   ntot)
            ntot = len(k_array)
            d[channel]['k_array'] = k_array
            d[channel]['k_weights'] = k_weights
            d[channel]['ntot'] = ntot
            
            # Load SRG transformation in 1S0 and 3S1-3D1 channels
            H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
            H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
                                             method='srg', generator='Wegner', 
                                             lamb=lamb)
            # Load U(k, k') [unitless]
            U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
            
            # Divide out momentum/weights factor
            factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
            if channel == '3S1':
                factor_array = np.concatenate( (factor_array, factor_array) )
            row, col = np.meshgrid(factor_array, factor_array)
            
            # Converting to units [fm^3]
            I_matrix_unitless = np.eye( len(factor_array), len(factor_array) )
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            
            # # Store both matrices in the dictionary
            # d[channel]['I'] = I_matrix_unitless
            # d[channel]['delta_U'] = delta_U_matrix_unitless
            
            # U_matrix = U_matrix_unitless / row / col
            I_matrix = I_matrix_unitless / row / col
            delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # Store both matrices in the dictionary
            d[channel]['I'] = I_matrix # [fm^3]
            d[channel]['delta_U'] = delta_U_matrix # [fm^3]
            
        # Save dictionary
        self.d = d
        
    
    def n_lambda_pp(self, q, k_F):
        # don't forget docstring
        # should be option somewhere to angle-average or set Q=0
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['1S0']['k_array']
        k_weights = self.d['1S0']['k_weights']
        ntot = self.d['1S0']['ntot']
        #I_matrix = self.d['1S0']['I']
        delta_U_matrix = self.d['1S0']['delta_U'] # fm^3
        
        # Find index of q in k_array
        q_index = find_q_index(q, k_array)
        # weight = 2/np.pi * k_array[q_index]**2 * k_weights[q_index] # fm^-3
        
        # T = 1, M_T = 1/2 + 1/2 = 1
        
        # Only one k_F in this case
        
        # First three terms (combining second and third into one term linear
        # in \delta U)
        if q < k_F:
            
            # Factor of 2 for sum over \sigma and \sigma', units fm^3
            #first_term = 2 / weight
            first_term = 2
            # Units fm^3
            middle_terms = 2 / np.pi * delta_U_matrix[q_index, q_index]
            
            low_q_contribution = first_term + middle_terms
            
        else:
            
            low_q_contribution = 0
            
        # Fourth term
        
        # Index of where to stop the integral in the momentum mesh k_array
        k_F_cutoff = find_q_index(k_F, k_array)
        if k_array[k_F_cutoff] > k_F:
            k_F_cutoff -= 1 # Prevents overshooting
            
        # fourth_term_integrand = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot]
        # # Integration measure 2/\pi k^2 dk is implicitly in the delta_U's here

        # high_q_contribution = 1/4 * 2/np.pi * \
        #                       np.sum( fourth_term_integrand[:k_F_cutoff+1] )
        
        fourth_term_integrand = 2/np.pi * k_array**2 * k_weights * \
                                delta_U_matrix[:ntot, q_index] * \
                                delta_U_matrix.T[q_index, :ntot]
        
        high_q_contribution = 1/4 * 2/np.pi * \
                              np.sum( fourth_term_integrand[:k_F_cutoff+1] )
                              
        # return low_q_contribution
        # return high_q_contribution
        return low_q_contribution + high_q_contribution
    
    
    def n_lambda_pn(self, q, kF1, kF2):
        # don't forget docstring
        # should be option somewhere to angle-average or set Q=0
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['1S0']['k_array']
        k_weights = self.d['1S0']['k_weights']
        ntot = self.d['1S0']['ntot']
        #I_matrix = self.d['1S0']['I']
        delta_U_1S0 = self.d['1S0']['delta_U'] # fm^3
        delta_U_3S1 = self.d['3S1']['delta_U'] # fm^3
        
        # Find index of q in k_array
        q_index = find_q_index(q, k_array)
        # weight = 2/np.pi * k_array[q_index]**2 * k_weights[q_index] # fm^-3
        
        k_F = min( (kF1, kF2) )
        
        # T = 1 or 0, M_T = 1/2 + 1/2 = 0
        # We get 1S0 and 3S1-3D1 contributions here

        # First three terms (combining second and third into one term linear
        # in \delta U)
        # This if statement is the two \theta functions
        if q < k_F:
            
            # Factor of 2 for sum over \sigma and \sigma', units fm^3
            #first_term = 2 / weight
            first_term = 2
            # Units fm^3
            middle_terms = 2 / np.pi * ( delta_U_1S0[q_index, q_index] / 2 + \
                           3/2 * delta_U_3S1[q_index, q_index] )            
            
            low_q_contribution = first_term + middle_terms
            
        else:
            
            low_q_contribution = 0
            
        # Fourth term
        
        # Index of where to stop the integral in the momentum mesh k_array
        k_F_cutoff = find_q_index(k_F, k_array)
        if k_array[k_F_cutoff] > k_F:
            k_F_cutoff -= 1 # Prevents overshooting
            
        # fourth_term_integrand = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot]
        # # Integration measure 2/\pi k^2 dk is implicitly in the delta_U's here

        # high_q_contribution = 1/4 * 2/np.pi * \
        #                       np.sum( fourth_term_integrand[:k_F_cutoff+1] )
        
        fourth_term_integrand = 2/np.pi * k_array**2 * k_weights * (\
                                delta_U_1S0[:ntot, q_index] * \
                                delta_U_1S0.T[q_index, :ntot]/4 + \
                                3/4 * delta_U_3S1[:ntot, q_index] * \
                                delta_U_3S1.T[q_index, :ntot] + \
                                3/4 * delta_U_3S1[:ntot, ntot+q_index] * \
                                delta_U_3S1.T[ntot+q_index, :ntot] )
        
        high_q_contribution = 1/2 * 2/np.pi * \
                              np.sum( fourth_term_integrand[:k_F_cutoff+1] )
                              
        # return low_q_contribution
        # return high_q_contribution
        return low_q_contribution + high_q_contribution


if __name__ == '__main__':

    
    # --- Set up --- #
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    # # Gezerlis N2LO 1 fm with \lambda=1.35 fm^-1
    # kvnn = 222
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120
    # kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    pmd = pair_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
    
    # Details of example nuclei (format is [nuclei, Z, N])
    # nuclei_list = ['O16', 8, 8]
    # nuclei_list = ['Ca40', 20, 20]
    # nuclei_list = ['Ca48', 20, 28]
    # nuclei_list = ['Pb208', 82, 126]
    
    d = {}
    
    nuclei_details = [ ['O16', 8, 8], ['Ca40', 20, 20], ['Ca48', 20, 28],
                       ['Pb208', 82, 126] ]
    for nuclei_list in nuclei_details:

        nucleus = nuclei_list[0]
        Z = nuclei_list[1]
        N = nuclei_list[2]
        A = N + Z
        r_array, rho_p_array = load_density(nucleus, 'proton', Z, N)
        r_array, rho_n_array = load_density(nucleus, 'neutron', Z, N)
        
        d[nucleus] = {}
    
        # Call LDA class
        lda = LDA(r_array, rho_p_array, rho_n_array)
    
    
    # --- Calculate momentum distributions --- #
    
    # n(q, Q=0) pp-pair
    # Should this be divided by # of pairs?
        n_pp_array = lda.local_density_approximation( q_array, pmd.n_lambda_pp, 
                                                     'pp' )
        n_pn_array = lda.local_density_approximation( q_array, pmd.n_lambda_pn, 
                                                     'pn' )
        
        d[nucleus] = n_pp_array / n_pn_array
        
        
    
    ### BUG TESTING ###
    # n_pp_array *= (Z-1)/(A-1)
    # n_pp_array /= factor_array
    # n_pp_array *= factor_array
    
    
    # --- Plot --- #
    
    # # n(q, Q=0) pairs
    #     plt.clf()
    #     f, ax = plt.subplots( 1, 1, figsize=(4, 4) )
    
    #     plt.semilogy(q_array, n_pp_array, color='blue', linestyle='solid',
    #                  label=nucleus+' (pp)')
    #     plt.semilogy(q_array, n_pn_array, color='red', linestyle='solid',
    #                  label=nucleus+' (pn)')
    #     plt.ylabel(r'$<n_{\lambda}(q,Q=0)>$' + ' [fm' + r'$^{6}$' + ']')
    #     ax.set_xlim( [min(q_array), 5] )
    #     ax.set_ylim( [1e-3, 1e4])
    
    #     ax.legend(loc=0, frameon=False)
    #     ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    #     lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    #     lambda_label_location = 'center right'
    #     anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
    #                               frameon=False)
    #     ax.add_artist(anchored_text)
    
    #     plt.show()
    

    # ratio of n(q, Q=0) pairs
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )
    
    plt.semilogy(q_array, d['O16'], color='black',
                 linestyle='solid', label='O16')
    plt.semilogy(q_array, d['Ca40'], color='red',
                 linestyle='dashdot', label='Ca40')
    plt.semilogy(q_array, d['Ca48'], color='blue',
                 linestyle='dashed', label='Ca48')
    plt.semilogy(q_array, d['Pb208'], color='green',
                 linestyle='dotted', label='Pb208')
    plt.ylabel(r'$<n_{\lambda}(q,Q=0)>$' + ' (pp/pn ratios)')
    ax.set_xlim( [min(q_array), 5] )
    ax.set_ylim( [1e-2, 2e0])
    
    ax.legend(loc='lower right', frameon=False)
    ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    lambda_label_location = 'lower left'
    anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
                                  frameon=True)
    ax.add_artist(anchored_text)
    
    plt.show()