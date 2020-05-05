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
#
#------------------------------------------------------------------------------


# Description of this test:
#   Testing mesh converting function.


from run_srg import run_srg
import numpy as np


kvnns = [900, 901, 902]
channel = '3S1'
kmax = 10.0
kmid = 2.0
ntot = 120

generators = ['Wegner', 'Block-diag']
lambda_array = np.array( [10.0, 2.8, 2.0, 1.2] )

for kvnn in kvnns:
    for generator in generators:
        if generator == 'Block-diag':
            for lamb in lambda_array:
                d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, 
                            np.array( [lambda_array[-1]] ), lambda_bd=lamb,
                            save=True)
        else:
            d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, 
                        lambda_array, save=True)