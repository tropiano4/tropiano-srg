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
# Scripts made by A.T.
from run_srg import run_srg


#kvnn = 902
kvnn = 10
# Could add the RKE 1-20 fm^-1 LO potentials as 903, 904, ...
channel = '3S1'
#kmax = 30.0
kmax = 8.0
#kmid = 4.0
kmid = 2.0
ntot = 120

#lambda_array = np.array( [6.0, 3.0, 2.0, 1.5, 1.0] )
lambda_array = np.array( [6.0, 3.0, 2.0, 1.5] )
#generator_list = ['Wegner', 'T']

#for generator in generator_list:
    
    #d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
                #save=True)
                
generator = 'Block-diag'

for lamb in lambda_array:
    
    d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
                lambda_bd=lamb, save=True)