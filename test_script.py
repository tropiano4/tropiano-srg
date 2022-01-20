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
#   Test the ratio of the cross section for the high energy nuclear
#   photoeffect to that for the deuteron photoeffect using just the squares of
#   the quasi-deuteron and deuteron wave functions as in Levinger:1951vp.

#   Here we will assume the quasi-deuteron wave function squared is given by
#   the nuclear momentum distribution squared using only 3S1-3D1.


import numpy as np
# Scripts made by A.T.
from densities import load_density
from dmd import deuteron_momentum_distributions
from misc.integration import gaussian_quadrature_mesh
# from pmd import pair_momentum_distributions
from potentials.vsrg_macos import vnn
from snmd import single_nucleon_momentum_distributions
import time


# Summary of what's done in Levinger:1951vp
#   - They take \sigma_qd / \sigma_d = |\psi_k|^2 / |\psi_d|^2
#     where \psi_k is the quasideuteron wave function.
#   - There is dependence on several parameters in the ratio: \alpha (inverse
#     scattering length), r_0 (effective range), k (relative momentum), and
#     v (volume of the nucleus).
#   - For k-dependence they average over all possible values of k assuming
#     Fermi distributions for proton and neutron k values up to some maximum
#     wave number 1 fm^-1.
#   - Note, this is less than \lambda though. Try a few different things.

# --- set-up --- #

# Potential number
kvnn = 6 # AV18

# Channels to include in calculation
channels = ['3S1'] # Triplet S only
# channels = ['1S0', '3S1']

# SRG \lambda value and generator
lamb = 1.35
generator = 'Wegner'

# Details of the momentum mesh
kmax, kmid, ntot = 15.0, 3.0, 120
# Momenta and weights
k_array, k_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)

# Nucleus formatted as (Name, Z, N)
# nucleus = ('He4', 2, 2)
# nucleus = ('C12', 6, 6)
# nucleus = ('O16', 8, 8)
# nucleus = ('Ca40', 20, 20)
# nucleus = ('Ca48', 20, 28)
nucleus = ('Pb208', 82, 126)

nucleus_name = nucleus[0]
Z = nucleus[1]
N = nucleus[2]
A = N+Z

# Energy density functional
edf = 'SLY4'

# Get densities
R_array, rho_p_array = load_density(nucleus_name, 'proton', Z, N, edf)
R_array, rho_n_array = load_density(nucleus_name, 'neutron', Z, N, edf)
dR = R_array[2] - R_array[1]


# --- Calculate quasi-deuteron momentum distribution --- #

# # Initialize pair-nucleon momentum distribution class
# pmd = pair_momentum_distributions(kvnn, channels, lamb, kmax, kmid, ntot)
        
# # Quasideuteron momentum distribution as a function of k
# qd_array = 2*pmd.n_Q0(k_array, R_array, dR, rho_p_array, rho_n_array)

# Initialize single-nucleon momentum distribution class
snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax, kmid,
                                             ntot)
        
# Try using interpolated version first
try:
    n_p_func, _, _, _ = snmd.n_lambda_interp(nucleus_name, 'proton', Z, N, edf)
    n_n_func, _, _, _ = snmd.n_lambda_interp(nucleus_name, 'neutron', Z, N,
                                             edf)
# Need to generate the files first
except OSError:
    t0 = time.time()
    snmd.write_file(nucleus_name, 'proton', Z, N, edf)
    snmd.write_file(nucleus_name, 'neutron', Z, N, edf)
    t1 = time.time()
    mins = (t1-t0)/60
    print(f'Done after {mins:.5f} minutes.')
    n_p_func, _, _, _ = snmd.n_lambda_interp(nucleus_name, 'proton', Z, N, edf)
    n_n_func, _, _, _ = snmd.n_lambda_interp(nucleus_name, 'neutron', Z, N, edf)

# Calculate each distribution for momenta k
n_p_array = n_p_func(k_array)
n_n_array = n_n_func(k_array)

# # Combine proton+neutron and divide by A for quasideuteron
# n_qd_array = (n_p_array+n_n_array) / A # This is what is done for a_2

# # Normalization
# norm_qd = 4*np.pi/(2*np.pi)**3 * np.sum(k_weights*k_array**2*n_qd_array)
# print(f'Normalization of n_qd(q)={norm_qd:.5f}')


# --- Calculate deuteron momentum distribution --- #

# Initialize deuteron momentum distributions class
dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot,
                                      interp=True)
    
# Get interpolated functions of deuteron momentum distribution
# Ignore the 1, \delta U, and \delta U^2 isolated contributions (take total only)
n_d_func, _, _, _ = dmd.n_lambda_interp()
    
# # Calculate deuteron momentum distribution
# n_d_array = n_d_func(k_array)

# # Normalization 
# norm_d = 4*np.pi/(2*np.pi)**3 * np.sum(k_weights*k_array**2*n_d_array)
# print(f'Normalization of n_d(q)={norm_d:.5f}')


# --- Evaluate ratio --- #

# Levinger gets \sigma_qd / \sigma_d = 6.4*(N*Z)/A = 1.6*A (if N=Z for latter)
factor = N*Z/A
# The way I'm evaluating this, I should expect something near 1.6

# Take one value in k (loop over several options)
k_values = np.linspace(1.0, 3.0, 3)
print('--- Evaluate ratio at one value of k ---\n')
for k in k_values:
    ratio = ( n_p_func(k)+n_n_func(k) ) / ( factor * n_d_func(k) )
    print(f'k={k:.1f}, ratio={ratio:.5f}')

# Integrate over an interval in k [0, k_max] where k_max is the above k_values
print('\n--- Evaluate ratio integrating over [0, k_max] ---\n')
for k_max in k_values:
    
    # Get Gaussian-quadrature mesh
    k_array, k_weights = gaussian_quadrature_mesh(k_max, 100)
    
    # # Check it gives back a_2 (you have to change channels to include 1S0)
    # # This is confirmed
    # k_array, k_weights = gaussian_quadrature_mesh(4.5, 40, xmin=3.8)
    
    # Compute numerator
    n_p_array = n_p_func(k_array)
    n_n_array = n_n_func(k_array)

    # Combine proton+neutron and divide by A for quasideuteron
    n_qd_array = (n_p_array+n_n_array) / factor
    numerator = np.sum(k_weights * k_array**2 * n_qd_array)
    
    # Compute denominator
    n_d_array = n_d_func(k_array)
    denominator = np.sum(k_weights * k_array**2 * n_d_array)
    
    ratio = numerator / denominator
    print(f'k_max={k_max:.1f}, ratio={ratio:.5f}')