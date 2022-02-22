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
#   Generate momentum distribution files for other potentials.


import numpy as np
# Scripts made by A.T.
from dmd import deuteron_momentum_distributions
from pmd import pair_momentum_distributions
from snmd import single_nucleon_momentum_distributions
import time


# --- Set-up --- #

# kvnns = [111, 113, 222, 224]
kvnns = [110, 111, 112, 113]
channels = ['1S0', '3S1']
# channels = ['1S0', '3S1', '3P0', '3P1', '1P1']
kmax, kmid, ntot = 15.0, 3.0, 120
lamb = 1.35
# nuclei = ( ('He4', 2, 2), ('C12', 6, 6), ('O16', 8, 8), ('Ca40', 20, 20),
#            ('Ca48', 20, 28), ('Fe56', 26, 30), ('Pb208', 82, 126) )
# nuclei = [
#     ('Li7', 3, 4), ('Ti48', 22, 26), ('Ag107', 47, 60), ('Sn118', 50, 68),
#     ('Ce140', 58, 82), ('Ta181', 73, 108), ('U238', 92, 146)
# ]
# Nuclei from Gogny densities
nuclei = [
    ('He4', 2, 2), ('Li7', 3, 4), ('Be9', 4, 5), ('C12', 6, 6), ('O16', 8, 8),
    ('Al27', 13, 14), ('Ca40', 20, 20), ('Ca48', 20, 28), ('Ti48', 22, 26),
    ('Fe56', 26, 30), ('Cu63', 29, 34), ('Ag107', 47, 60), ('Sn118', 50, 68),
    ('Ce140', 58, 82), ('Ta181', 73, 108), ('Au197', 79, 118),
    ('Pb208', 82, 126), ('U238', 92, 146)
]
# edf = 'SLY4'
edf = 'Gogny'


# # --- SRG --- #
# from run_srg import run_srg

# lambda_array = np.array( [6.0, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.35] )

# for kvnn in kvnns:
#     for channel in channels:
#         d = run_srg(kvnn, channel, kmax, kmid, ntot, 'Wegner', lambda_array)


# --- Momentum distributions --- #

for kvnn in kvnns:
    
    print(f'Starting kvnn = {kvnn}.')
    
    # Start timer
    t0 = time.time()
    
    # --- Deuteron --- #
    
    # Initialize deuteron momentum distribution class
    dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
    
    # Write deuteron momentum distribution file
    dmd.write_file()
    
    print('Done with deuteron.')
    
    # Initialize single-nucleon and pair momentum distributions classes
    snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax,
                                                  kmid, ntot)
    pmd = pair_momentum_distributions(kvnn, channels, lamb, kmax, kmid, ntot)
    
    for nucleus in nuclei:
        
        # Details of nucleus
        nucleus_name = nucleus[0]
        Z = nucleus[1]
        N = nucleus[2]
        
        # --- Single-nucleon --- #
        snmd.write_file(nucleus_name, 'proton', Z, N, edf)
        snmd.write_file(nucleus_name, 'neutron', Z, N, edf)
    
        # --- Pair --- #
        pmd.write_file(nucleus_name, 'pn', Z, N, edf)
        pmd.write_file(nucleus_name, 'pp', Z, N, edf)
        pmd.write_file(nucleus_name, 'nn', Z, N, edf)
        
        print(f'Done with {nucleus_name}.')
    
    # End timer
    t1 = time.time()
    mins = (t1-t0)/60
    
    # Print time for one kvnn
    print(f'kvnn = {kvnn} done after {mins:.2f} minutes.\n')