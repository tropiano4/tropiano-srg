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
#   Generating momentum distributions for Gogny densities.


import time
# Scripts made by A.T.
# from dmd import deuteron_momentum_distributions
from figures import figures_functions as ff
from pmd import pair_momentum_distributions
from snmd import single_nucleon_momentum_distributions


# - Set-up - #
# Potentials
# kvnns_list = [6, 222, 224]
kvnns_list = [6]

# Generate single-nucleon and pair momentum distributions

# Channels to include in calculation (S-waves only or higher partial waves)
# channels_list_av18 = [ ('1S0', '3S1'), ('1S0', '3S1', '3P0', '1P1', '3P1') ]
# channels_list_gez = [ ('1S0', '3S1') ]
channels_list_gogny = [ ('1S0', '3S1') ]

# SRG \lambda values
# lambdas_list_6_222 = [1.35, 1.5, 2.0, 3.0, 6.0]
# lambdas_list_224 = [1.35]
lambdas_list_gogny = [1.35]
    
# Momentum mesh details
kmax, kmid, ntot = 15.0, 3.0, 120 # Default

# Nuclei to calculate
# nuclei_list = [ ('He4', 2, 2), ('He8', 2, 6), ('Be9', 4, 5), ('C12', 6, 6),
#                 ('O16', 8, 8), ('Ca40', 20, 20), ('Ca48', 20, 28),
#                 ('Fe56', 26, 30), ('Pb208', 82, 126) ]
nuclei_list_gogny = [ ('He4', 2, 2), ('Be9', 4, 5), ('C12', 6, 6),
                      ('O16', 8, 8), ('Al27', 13, 14), ('Ca40', 20, 20),
                      ('Ca48', 20, 28), ('Fe56', 26, 30), ('Cu63', 29, 34),
                      ('Au197', 79, 118), ('Pb208', 82, 126)]
    

# - Generate all data for single-nucleon and pair momentum distributions - #
for kvnn in kvnns_list:
        
    # if kvnn == 6:
    #     channels_list = channels_list_av18
    # else:
    #     channels_list = channels_list_gez
    channels_list = channels_list_gogny
    t0_k = time.time()
    
    for ic, channels in enumerate(channels_list):
        
        # if kvnn == 224:
        #     lambdas_list = lambdas_list_224
        # else:
        #     lambdas_list = lambdas_list_6_222
        lambdas_list = lambdas_list_gogny
        t0_c = time.time()
            
        for lamb in lambdas_list:
                
            t0_l = time.time()
            
            # Initialize classes
            snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb,
                                                         kmax, kmid, ntot)
            pmd = pair_momentum_distributions(kvnn, channels, lamb, kmax, kmid,
                                              ntot)
                
            # for nuclei in nuclei_list:
            for nuclei in nuclei_list_gogny:
                    
                t0_N = time.time()

                nucleus = nuclei[0]
                Z = nuclei[1]
                N = nuclei[2]
                # if Z < 6:
                #     edf = 'AV18'
                # else:
                #     edf = 'SLY4'
                edf = 'Gogny'
                    
                # Write single-nucleon files for given nucleus
                snmd.write_file(nucleus, 'proton', Z, N, edf)
                snmd.write_file(nucleus, 'neutron', Z, N, edf)
                    
                # # Write pair files for given nucleus
                # pmd.write_file(nucleus, 'pp', Z, N, edf)
                # pmd.write_file(nucleus, 'nn', Z, N, edf)
                # pmd.write_file(nucleus, 'pn', Z, N, edf)
                    
                # Time for each nucleus to run
                t1_N = time.time()
                mins_N = (t1_N-t0_N)/60
                print('\t\t\tDone with %s after %.5f minutes.' % (nucleus,
                                                                  mins_N) )
            
            # Time for each \lambda to run   
            t1_l = time.time()
            mins_l = (t1_l-t0_l)/60
            print( '\n\t\tDone with \lambda=%s after %.5f minutes.\n' % (
                   ff.convert_number_to_string(lamb), mins_l) )
               
        # Time for each channels to run
        t1_c = time.time()
        hours_c = (t1_c-t0_c)/3600
        if ic == 0:
            print('\tDone with S-waves after %.3f hours.\n' % hours_c)
        else:
            print('\tDone with all channels after %.3f hours.\n' % hours_c)
        
    # Time for each potential to run
    t1_k = time.time()
    hours_k = (t1_k-t0_k)/3600
    print( 'Done with kvnn=%d after %.5f hours.\n' % (kvnn, hours_k) )
    
    
# # - Generate all data for deuteron momentum distributions - #
# for kvnn in kvnns_list:
        
#     t0_k = time.time()
        
#     for lamb in lambdas_list_6_222:
                
#         t0_l = time.time()
        
#         # Initialize class
#         dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
        
#         # Write deuteron files
#         dmd.write_file()
        
#         # Time for each \lambda to run   
#         t1_l = time.time()
#         mins_l = (t1_l-t0_l)/60
#         print( '\n\t\tDone with \lambda=%s after %.5f minutes.\n' % (
#                ff.convert_number_to_string(lamb), mins_l) )
        
#     # Time for each potential to run
#     t1_k = time.time()
#     mins_k = (t1_k-t0_k)/60
#     print( 'Done with kvnn=%d after %.5f minutes.\n' % (kvnn, mins_k) )
