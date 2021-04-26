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
#                nuclei using LDA. Created single_particle_momentum_dist.py
#                in Old_codes.
#   04/14/21 --- Creating AV18 SRG evolution figure for APS April Meeting
#                presentation.
#                potential_contours_kvnn_6_channel_1P1_Wegner_lamb_1p5.png in
#                Figures/Operator_evolution/Old_figures/Potential_contours.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Trying to get correct deuteron SNMD normalization in LDA for reliable A/d
#   ratios.


import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from dmd import deuteron_momentum_distributions
from Potentials.vsrg_macos import vnn


# Load AV18 data
av18_data = np.loadtxt('av18_deuteron_data_temp.txt')
q_array_av18 = av18_data[:, 0]
n_p_array_av18 = av18_data[:, 3]
normalization_av18 = 4*np.pi/(2*np.pi)**3 * \
                     np.sum(0.1*q_array_av18**2*n_p_array_av18)
print('AV18 normalization = %.5f' % normalization_av18)

# Compute proton momentum distribution under LDA
kvnn, lamb, kmax, kmid, ntot = 6, 1.35, 10.0, 2.0, 120
q_array, q_weights = vnn.load_momentum(kvnn, '3S1', kmax, kmid, ntot)
dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
n_p_array = dmd.local_density_approximation(q_array, 'single-nucleon')
normalization = 4*np.pi/(2*np.pi)**3 * np.sum(q_weights*q_array**2*n_p_array)
print('LDA normalization = %.5f' % normalization)

# Plot
plt.semilogy(q_array, n_p_array, 'r-', label='LDA')
plt.semilogy(q_array_av18, n_p_array_av18, 'k:', label='AV18')

plt.xlabel('q fm^-1')
plt.ylabel('n(q) fm^3')
plt.xlim( (0.0, 5.0) )
plt.ylim( (1e-3, 1e4) )
plt.legend(loc='upper right')