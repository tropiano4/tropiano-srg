#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a separate file with an
# extension _test.py. For example, momentum_projection_operator_testv1.py (v1
# means 'version 1'). Use the revision history below to document when and why
# these files are created.
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
#                nuclei using LDA. Created single_particle_momentum_dist.py
#                in Old_codes.
#   04/14/21 --- Creating AV18 SRG evolution figure for APS April Meeting
#                presentation.
#                potential_contours_kvnn_6_channel_1P1_Wegner_lamb_1p5.png in
#                Figures/Operator_evolution/Old_figures/Potential_contours.
#   04/28/21 --- Testing normalization and contributions of \delta U, etc. or
#                pp/pn to single-nucleon momentun distributions. Created
#                lda_normalizations_test.py in Old_codes.
#   05/04/21 --- Testing higher partial waves of SRG transformations: 3P2-3F2
#                and 3D3-3G3 have numerical artifacts. Created
#                high_partial_waves_srg_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Test cancellation of
#     \delta U + \delta U^\dagger + 1/2 \sum \delta U \delta U^\dagger
#   for 1S0 and 3S1-3D1 channels.


import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Potentials.vsrg_macos import vnn
# from SRG.srg_unitary_transformation import SRG_unitary_transformation
from dmd import deuteron_momentum_distributions


# # Set up
# kvnn = 6
# channel = '1S0'
# # channel = '3S1'
# kmax, kmid, ntot = 15.0, 3.0, 120
# lamb = 1.35

# # Load momentum and weights
# k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
# # For dividing out momenta/weights
# factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
# # For coupled-channel matrices
# if vnn.coupled_channel(channel):
#     factor_array = np.concatenate( (factor_array, factor_array) )
        
# # Load SRG transformation
# H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
# H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, method='srg',
#                                  generator='Wegner', lamb=lamb)
# # Load U(k, k') [unitless]
# U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)

# # Isolate 2-body term and convert to fm^3
# if vnn.coupled_channel(channel):
#     I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
# else:
#     I_matrix_unitless = np.eye(ntot, ntot)
# row, col = np.meshgrid(factor_array, factor_array)

# delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
# # delta_U_matrix = delta_U_matrix_unitless / row / col

# # Compute \delta U \delta U^\dagger term
# delU2_unitless = delta_U_matrix_unitless @ delta_U_matrix_unitless.T
# # delU2 = delU2_unitless / row / col

# # Calculate part that should cancel
# # Factor of 1/2?
# cancel_matrix_unitless = delta_U_matrix_unitless + \
#                          delta_U_matrix_unitless.T + delU2_unitless
# print(cancel_matrix_unitless)

# Notes for next test
# 1. Compare deuteron LDA vs deuteron exact (plot and calculation
# normalizations)

kvnn = 6
kmax, kmid, ntot = 15.0, 3.0, 120
lamb = 1.35

# Load momentum and weights
q_array, q_weights = vnn.load_momentum(kvnn, '3S1', kmax, kmid, ntot)
factor_array_exact = 2/np.pi * q_weights * q_array**2
factor_array_lda = 4*np.pi/(2*np.pi)**3 * q_weights * q_array**2
factor_diff = 4*np.pi/(2*np.pi)**3 * np.pi/2

dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)

# Initialize arrays
n_tot_exact_array = np.zeros(ntot)
n_1_exact_array = np.zeros(ntot)
n_delU_exact_array = np.zeros(ntot)
n_delU2_exact_array = np.zeros(ntot)
n_tot_lda_array = np.zeros(ntot)
n_1_lda_array = np.zeros(ntot)
n_delU_lda_array = np.zeros(ntot)
n_delU2_lda_array = np.zeros(ntot)

for iq, q in enumerate(q_array):
    
    n_tot_exact, n_1_exact, n_delU_exact, n_delU2_exact = dmd.n_lambda_pair_exact(
                                            q, contributions='q_contributions')
    n_tot_exact_array[iq] = n_tot_exact
    n_1_exact_array[iq] = n_1_exact
    n_delU_exact_array[iq] = n_delU_exact
    n_delU2_exact_array[iq] = n_delU2_exact
    
    expectation_values = dmd.local_density_approximation(q_array, 'pair',
                                               contributions='q_contributions')
    n_tot_lda_array[iq] = expectation_values[iq, 0]
    n_1_lda_array[iq] = expectation_values[iq, 1]
    n_delU_lda_array[iq] = expectation_values[iq, 2]
    n_delU2_lda_array[iq] = expectation_values[iq, 3]
    
# Print normalizations
print('_'*50)
print('Exact total = %.5f' % np.sum(n_tot_exact_array*factor_array_exact))
print('Exact 1 term = %.5f' % np.sum(n_1_exact_array*factor_array_exact))
print('Exact \delta U term = %.5f' % np.sum(n_delU_exact_array*factor_array_exact))
print('Exact \delta U^2 term = %.5f' % np.sum(n_delU2_exact_array*factor_array_exact))
print('LDA total = %.5f' % np.sum(n_tot_lda_array*factor_array_lda))
print('LDA 1 term = %.5f' % np.sum(n_1_lda_array*factor_array_lda))
print('LDA \delta U term = %.5f' % np.sum(n_delU_lda_array*factor_array_lda))
print('LDA \delta U^2 term = %.5f' % np.sum(n_delU2_lda_array*factor_array_lda))

plt.clf()
plt.semilogy(q_array, n_tot_exact_array, label='Exact total',
             linestyle='solid', color='k')
plt.semilogy(q_array, n_1_exact_array, label='Exact 1', linestyle='dotted',
             color='b')
plt.semilogy(q_array, abs(n_delU_exact_array), linestyle='dotted', color='g',
             label='Exact ' + r'$|\delta U|$')
plt.semilogy(q_array, n_delU2_exact_array, linestyle='dotted', color='r',
             label='Exact ' + r'$\delta U^2$')

plt.semilogy(q_array, n_tot_lda_array*factor_diff, label='LDA total',
             linestyle='dashdot', color='tab:gray')
plt.semilogy(q_array, n_1_lda_array*factor_diff, label='LDA 1', linestyle='dashed',
             color='tab:blue')
plt.semilogy(q_array, abs(n_delU_lda_array)*factor_diff, linestyle='dashed',
             color='tab:green', label='LDA ' + r'$|\delta U|$')
plt.semilogy(q_array, n_delU2_lda_array*factor_diff, linestyle='dashed',
             color='tab:red', label='LDA ' + r'$\delta U^2$')

# legend_size = 14
# legend_location = 'upper left'
# plt.legend(bbox_to_anchor=(1.05, 1), loc=legend_location, borderaxespad=0.,
#            fontsize=legend_size)
legend_size = 12
legend_location = 'upper right'
plt.legend(loc=legend_location, fontsize=legend_size, ncol=2)

x_label = 'q [fm' + r'$^{-1}$' + ']'
plt.xlabel(x_label)
plt.ylabel(r'$n_d(q)$'+' [fm'+r'$^3$' + ']')

plt.xlim((0.0, 6.0))
plt.ylim((1e-5, 1e3))

plt.show()
