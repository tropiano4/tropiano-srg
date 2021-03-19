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
#   Calculate n_{\lambda}(q) for deuteron using LDA.


import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
import lda
import observables as ob
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation
from pmd_deuteron_test import deuteron_pair_momentum_distribution_v1


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


class deuteron(object):
    
    
    def __init__(self, r_array, kvnn, lamb, kmax=0.0, kmid=0.0, ntot=0):
        
        self.r_array = r_array

        k_array, k_weights = vnn.load_momentum(kvnn, '3S1', kmax, kmid, ntot)
        self.k_array = k_array
        self.k_weights = k_weights
        self.ntot = len(k_array)
        
        channel = '3S1'
        
        # Load SRG transformation in 1S0 and 3S1-3D1 channels
        H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
        H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
                                         method='srg', generator='Wegner', 
                                         lamb=lamb)
        # Load U(k, k') [unitless]
        U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
            
        # Divide out momentum/weights factor
        factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
        factor_array_long = np.concatenate( (factor_array, factor_array) )
        row, col = np.meshgrid(factor_array_long, factor_array_long)
            
        # Converting to units [fm^3]
        I_matrix_unitless = np.eye( 2*self.ntot, 2*self.ntot )
        delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            
        # # Store both matrices in the dictionary
        # d[channel]['delta_U'] = delta_U_matrix_unitless
            
        # U_matrix = U_matrix_unitless / row / col
        delta_U_matrix = delta_U_matrix_unitless / row / col
            
        # Save matrix
        self.delta_U = delta_U_matrix # [fm^3]
        
        # get deuteron density
        
        psi_k_unitless = ob.wave_function(H_initial, U=U_matrix_unitless)
    
        # Divide out momenta/weights
        psi_k = psi_k_unitless / factor_array_long
    
        hank_trans_3S1 = hankel_transformation('3S1', k_array, k_weights, r_array)
        hank_trans_3D1 = hankel_transformation('3D1', k_array, k_weights, r_array)
    
        # sign = -1
        psi_r_3S1 = hank_trans_3S1 @ psi_k[:ntot]
        psi_r_3D1 = hank_trans_3D1 @ psi_k[ntot:]
    
        self.rho_array = psi_r_3S1**2 + psi_r_3D1**2


    def n_lambda_pn(self, q, k_F):
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.k_array
        k_weights = self.k_weights
        ntot = self.ntot
        delta_U_3S1 = self.delta_U # fm^3
        
        # Find index of q in k_array
        q_index = find_q_index(q, k_array)
        # weight = 2/np.pi * k_array[q_index]**2 * k_weights[q_index] # fm^-3

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
            
            # # 2J+1 factor
            # middle_terms = 2 / np.pi * 3/2 * delta_U_3S1[q_index, q_index] 
            
            # no 2J+1 factor
            middle_terms = 2 / np.pi * 1/2 * delta_U_3S1[q_index, q_index]        
            
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
        
        # # 2J+1 factors
        # fourth_term_integrand = 2/np.pi * k_array**2 * k_weights * (\
        #                         3/4 * delta_U_3S1[:ntot, q_index] * \
        #                         delta_U_3S1.T[q_index, :ntot] + \
        #                         3/4 * delta_U_3S1[:ntot, ntot+q_index] * \
        #                         delta_U_3S1.T[ntot+q_index, :ntot] )
        
        # no 2J+1 factors
        fourth_term_integrand = 2/np.pi * k_array**2 * k_weights * (\
                                1/4 * delta_U_3S1[:ntot, q_index] * \
                                delta_U_3S1.T[q_index, :ntot] + \
                                1/4 * delta_U_3S1[:ntot, ntot+q_index] * \
                                delta_U_3S1.T[ntot+q_index, :ntot] )
        
        high_q_contribution = 1/2 * 2/np.pi * \
                              np.sum( fourth_term_integrand[:k_F_cutoff+1] )
                              
        # return low_q_contribution
        # return high_q_contribution
        return low_q_contribution + high_q_contribution


    def local_density_approximation(self, q_array):
        
        M = len(q_array)
    
        # Number of r_array points
        r_array = self.r_array
        rho_array = self.rho_array
        N = len(r_array)
        
        r2_array = r_array**2
        dr = 0.1 # Spacing between r-points
        # denominator = 4*np.pi * np.sum( dr * r2_array * rho_array )

        # pn pair
        kF_array = ( 3*np.pi**2 * rho_array )**(1/3)
        
        expectation_values = np.zeros(M)
        for i, q in enumerate(q_array):
    
            # Now evaluate f(q, kF) at each point in q_array and kF_array
            function_array = np.zeros(N)
            for j, k_F in enumerate(kF_array):

                function_array[j] = self.n_lambda_pn(q, k_F)
                    
            expectation_values[i] = 4*np.pi*dr * np.sum( r2_array * function_array )
  
        return expectation_values
    
    
if __name__ == '__main__':


    # --- How would deuteron work in LDA? --- #
    
    # # just need the r_array
    # r_array, _ = lda.load_density('O16', 'proton', 8, 8)
    # dr = 0.1
    
    # kvnn = 6
    # channel = '3S1'
    # lamb = 1.35
    # kmax, kmid, ntot = 10.0, 2.0, 120
    
    # # Load evolved deuteron wave function in momentum-space
    # k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
    # H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    # H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot,
    #                                  method='srg', generator='Wegner', lamb=lamb)
    # U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
    # psi_k_unitless = ob.wave_function(H_initial, U=U_matrix)
        
    # # Divide out momenta/weights
    # factor_array = np.concatenate( (np.sqrt(k_weights) * k_array,
    #                                 np.sqrt(k_weights) * k_array) ) * np.sqrt(2/np.pi)
    # psi_k = psi_k_unitless / factor_array
    
    # hank_trans_3S1 = hankel_transformation('3S1', k_array, k_weights, r_array)
    # hank_trans_3D1 = hankel_transformation('3D1', k_array, k_weights, r_array)
    
    # # sign = -1
    # psi_r_3S1 = hank_trans_3S1 @ psi_k[:ntot]
    # psi_r_3D1 = hank_trans_3D1 @ psi_k[ntot:]
    
    # rho_d = psi_r_3S1**2 + psi_r_3D1**2
    
    # # This is normalized to 1
    # normalization = np.sum(r_array**2 * rho_d) * dr
    # print(normalization)
    
    # kF_pn_array = ( 3*np.pi**2 * rho_d )**(1/3)
    # print(kF_pn_array)
    
    # --- Set up --- #
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120
    # kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    pmd = deuteron_pair_momentum_distribution_v1(kvnn, lamb, kmax, kmid, ntot)
    
    
    # --- Use deuteron to check pair momentum distribution --- #
    
    # Load evolved wave function here (unitless)
    H_initial = vnn.load_hamiltonian(kvnn, '3S1', kmax, kmid, ntot)
    H_evolved = vnn.load_hamiltonian(kvnn, '3S1', kmax, kmid, ntot, 
                                     method='srg', generator='Wegner', 
                                     lamb=lamb)
    U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
    
    psi_deuteron = ob.wave_function(H_initial, U=U_matrix_unitless)
    
    n_d_1 = np.zeros(ntot)
    n_d_23 = np.zeros(ntot)
    n_d_4 = np.zeros(ntot)
    for iq, q in enumerate(q_array):
        first_term, second_third_term, fourth_term = pmd.n_deuteron(q, psi_deuteron)
        n_d_1[iq] = first_term
        n_d_23[iq] = second_third_term
        n_d_4[iq] = fourth_term
    n_d_total_array = n_d_1 + n_d_23 + n_d_4
    
    # Correct momentum distribution for comparison
    n_d_exact = ( ob.wave_function(H_initial)[:ntot]**2 + \
                  ob.wave_function(H_initial)[ntot:]**2 ) / factor_array
        
    # Now calculate LDA deuteron
    
    # just need the r_array
    r_array, _ = lda.load_density('O16', 'proton', 8, 8)
    
    lda_deuteron = deuteron(r_array, kvnn, lamb, kmax, kmid, ntot)
    n_d_lda = lda_deuteron.local_density_approximation(q_array)
    # factors you might be missing
    overall_factor = 1/(2*np.pi)**3 * 4*np.pi / 2
    n_d_lda *= overall_factor


    # --- Plot n_lambda^d(q) --- #
    
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )

    ax.semilogy(q_array, n_d_total_array, color='black', linestyle='solid',
                label='High Res.')
    ax.semilogy(q_array, n_d_1, color='blue', linestyle='dotted', 
                label='1') 
    ax.semilogy(q_array, abs(n_d_23), color='purple', linestyle='dotted', 
                label=r'$|\delta U|$') 
    ax.semilogy(q_array, n_d_4, color='red', linestyle='dotted', 
                label=r'$|\delta U \delta U^{\dagger}|$') 
    ax.semilogy(q_array, n_d_exact, color='green', linestyle='dashed', 
                label='exact') 
    ax.semilogy(q_array, n_d_lda, color='gray', linestyle='dashdot',
                label='LDA')
    ax.set_ylabel(r'$<n^{\lambda}_d(q)>$' + ' [fm' + r'$^3$' + ']')
    ax.set_xlim( [min(q_array), 4] )
    ax.set_ylim( [1e-5, 1e3] )
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    lambda_label_location = 'lower right'
    anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
                                  frameon=False)
    ax.add_artist(anchored_text)

    plt.show()