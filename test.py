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
#   08/28/19 --- Testing the deuteron momentum distribution using projection
#                operators from operators.py
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   09/10/19 --- Using this script to run SRG evolution on several potentials.
#
#------------------------------------------------------------------------------


import numpy as np
from run_srg import run_srg


kvnn_list = [10, 106, 222]
channel_list = ['1S0', '1P1']
ntot = 120

generator = 'Block-diag'
lambda_array = np.array( [6.0, 3.0, 2.0, 1.5, 1.0] )
# Save the evolved Hamiltonian?
save = True


for kvnn in kvnn_list:
    
    if kvnn == 10:
        kmax = 30.0
        kmid = 4.0
    elif kvnn == 106:
        kmax = 8.0
        kmid = 2.0
    else:
        kmax = 10.0
        kmid = 2.0
    
    for channel in channel_list:
        
        for lamb in lambda_array:
            
            # Evolve Hamiltonian
            d = run_srg(kvnn, channel, kmax, kmid, ntot, generator,
                        lambda_array, lamb, save)