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


import numpy as np
# Scripts made by A.T.
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


# Set up
kvnn = 6
channel = '1S0'
# channel = '3S1'
kmax, kmid, ntot = 15.0, 3.0, 120
lamb = 1.35

# Load momentum and weights
k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
# For dividing out momenta/weights
factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
# For coupled-channel matrices
if vnn.coupled_channel(channel):
    factor_array = np.concatenate( (factor_array, factor_array) )
        
# Load SRG transformation
H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, method='srg',
                                 generator='Wegner', lamb=lamb)
# Load U(k, k') [unitless]
U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)

# Isolate 2-body term and convert to fm^3
if vnn.coupled_channel(channel):
    I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
else:
    I_matrix_unitless = np.eye(ntot, ntot)
row, col = np.meshgrid(factor_array, factor_array)

delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
# delta_U_matrix = delta_U_matrix_unitless / row / col

# Compute \delta U \delta U^\dagger term
delU2_unitless = delta_U_matrix_unitless @ delta_U_matrix_unitless.T
# delU2 = delU2_unitless / row / col

# Calculate part that should cancel
# Factor of 1/2?
cancel_matrix_unitless = delta_U_matrix_unitless + \
                         delta_U_matrix_unitless.T + delU2_unitless
print(cancel_matrix_unitless)

# Notes for next test
# 1. Compare deuteron LDA vs deuteron exact (plot and calculation
# normalizations)