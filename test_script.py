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
#                momentum_projection_operator_testv1.py based off last tests in
#                this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py.
#   01/26/21 --- Renamed to test_script.py.
#   03/16/21 --- Testing pair momentum distributions for different nuclei
#                using LDA with simple filled Fermi sea single-particle
#                momentum distributions. Created pmd_simple_test.py.
#   03/16/21 --- Testing pair momentum distribution for deuteron starting from
#                second quantization derivation and using expansion of
#                U(k, k'). Created pmd_deuteron_test.py.
#   03/16/21 --- Testing single-nucleon momentum distributions for different
#                nuclei using LDA. Created single_particle_momentum_dist.py.
#   04/14/21 --- Creating AV18 SRG evolution figure for APS April Meeting
#                presentation.
#                potential_contours_kvnn_6_channel_1P1_Wegner_lamb_1p5.png in
#                figures/operator_evolution/old/potential_contours.
#   04/28/21 --- Testing normalization and contributions of \delta U, etc. or
#                pp/pn to single-nucleon momentun distributions. Created
#                lda_normalizations_test.py.
#   05/04/21 --- Testing higher partial waves of SRG transformations: 3P2-3F2
#                and 3D3-3G3 have numerical artifacts. Created
#                high_partial_waves_srg_test.py.
#   06/10/21 --- Verifying \theta functions averaging in snmd.py and dmd.py by
#                comparing numerical functions to analytic evaluation of 
#                \int d3K \int d3k \theta(kF-|K/2+k|) \theta(kF-|K/2-k|).
#                Created theta_functions_test.py.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Testing Fourier transformation of AV18 wave functions.


import numpy as np

rho_data = np.loadtxt('densities/AV18/he4_densities_2_2.txt')
n_data = np.loadtxt('data/qmc/AV18_He4_snmd.txt')

r_array = rho_data[:, 0]
dr = r_array[1]-r_array[0]
rho_array = rho_data[:, 1]
print( np.sum(4*np.pi*dr*r_array**2*rho_array) )

k_array = n_data[:, 0]
dk = k_array[1]-k_array[0]
n_array = n_data[:, 1]
print( np.sum(4*np.pi/(2*np.pi)**3*dk*k_array**2*n_array) )

# Test 1: stupid way
k_mesh, r_mesh = np.meshgrid(k_array, r_array, indexing='ij')
_, rho_mesh = np.meshgrid(k_array, rho_array, indexing='ij')
exponential_mesh = np.exp( 1j * k_mesh * r_mesh )
n_array_new = np.sum( 4*np.pi*dr*r_mesh**2*rho_mesh*exponential_mesh, axis=-1)
print(n_array)
print( abs(n_array_new) )
print(n_array/abs(n_array_new))

# Are the densities and single-nucleon momentum distributions really just
# Fourier transforms of each other? I doubt it.
