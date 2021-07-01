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
#   06/10/21 --- Verifying \theta functions averaging in snmd.py and dmd.py by
#                comparing numerical functions to analytic evaluation of 
#                \int d3K \int d3k \theta(kF-|K/2+k|) \theta(kF-|K/2-k|).
#                Created theta_functions_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Looking for bug in He4, He8 snmd.py code.


import numpy as np
import time
# Scripts made by A.T.
from densities import load_density
from Misc.integration import gaussian_quadrature_mesh
from pmd import pair_momentum_distributions
from Potentials.vsrg_macos import vnn
from snmd import single_nucleon_momentum_distributions


# Set-up
kvnn = 6
channels = ('1S0', '3S1')
lamb = 1.35
kmax, kmid, ntot = 15.0, 3.0, 120

# Initialize classes
pmd = pair_momentum_distributions(kvnn, channels, lamb, kmax, kmid, ntot)
snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax, kmid,
                                             ntot)

# Get nucleonic densities
nucleus = 'He4'
Z = 2
N = 2
# nucleus = 'He8'
# Z = 2
# N = 6
R_array, rho_p_array = load_density(nucleus, 'proton', Z, N, 'AV18')
rho_n_array = rho_p_array
dR = R_array[1] - R_array[0]

# Get momentum (channel argument doesn't matter here)
q_array, _ = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)

# Evaluate kF values at each point in R_array to set max value of Q
kFp_array = (3*np.pi**2 * rho_p_array)**(1/3)
kFn_array = (3*np.pi**2 * rho_n_array)**(1/3)

# Get C.o.M. momentum for pair distribution
Q_max = max(kFp_array) + max(kFn_array)
ntot_Q = 50
Q_array, _ = gaussian_quadrature_mesh(Q_max, ntot_Q)

# Evaluate p distribution for each q
t0 = time.time()
n_p_array = snmd.n_total(q_array, R_array, dR, rho_p_array, rho_n_array)
t1 = time.time()
mins = (t1-t0)/60
print('Proton momentum distribution done after %.5f minutes.' % mins)

# Evaluate pn distribution for each q and Q
t0 = time.time()
n_pn_array = pmd.n_total(q_array, Q_array, R_array, dR, rho_p_array,
                         rho_n_array)
t1 = time.time()
mins = (t1-t0)/60
print('pn pair momentum distribution done after %.5f minutes.' % mins)