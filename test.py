#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes.
#
# Last thing tested:
#  Running run_magnus.py for the following: kvnn = 900-901, G= 'Wegner' and 
#  'T', and k_magnus = 2, 6, and 10.
#
#------------------------------------------------------------------------------


import numpy as np
from run_magnus import run_magnus


# Potential details
kvnn_list = [900, 901]
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120
    
# Evolution details
generator_list = ['Wegner', 'T']
lambda_array = np.array( [10.0] )
ds = 1e-6
#lambda_array = np.array( [2.8, 2.0, 1.2] )
#ds = 1e-5
k_magnus_list = [2, 6, 10]    

# Save the evolved Hamiltonian and omega?
save = True
    
# Loop over each option
for kvnn in kvnn_list:
    for generator in generator_list:
        for k_magnus in k_magnus_list:
            
            # Evolve Hamiltonian
            d = run_magnus(kvnn, channel, kmax, kmid, ntot, generator, 
                           lambda_array, k_magnus, ds, save)