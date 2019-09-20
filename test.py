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
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
from SRG_codes import srg_wegner


kvnn = 902
# Could add the RKE 1-20 fm^-1 LO potentials as 903, 904, ...
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120

lambda_array = np.array( [6.0, 3.0, 2.0, 1.5, 1.0] )
generator = 'Wegner'
#generator = 'T'

# Load initial Hamiltonian, kinetic energy, and weights
H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
T_rel = lp.load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)     
k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
# h-bar^2 / M [MeV fm^2] for conversion from MeV to scattering units
hbar_sq_over_M = 41.47
    
# Initial value of lambda in units fm^-1 (most potentials in 
# Potentials/vsrg_macos are generated at lambda = 12 fm^-1 but check 
# run_generate_vsrg_vlowk.pl to be sure)
lambda_initial = 12.0

evolve = srg_wegner.SRG(H_initial)
# NEED TO ADD ODEINT TO T IF IT WORKS
#evolve = srg_kinetic_energy.SRG(H_initial, T_rel)
    
# Time the evolution and return dictionary d of evolved Hamiltonians where
# the keys are lambda values
t0 = time.time() # Start time
d = evolve.evolve_hamiltonian(lambda_initial, lambda_array, method='odeint')
t1 = time.time() # End time
    
# Print details
mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
print('_'*85)
print( 'H(s) done evolving to final lambda = %.1f fm^-1 after %f minutes'
      % (lambda_array[-1], mins) )
print('_'*85)
print('\nSpecifications:\n')
print( 'kvnn = %d, channel = %s' % (kvnn, channel) )
print( 'kmax = %.1f, kmid = %.1f, ntot = %d' % (kmax, kmid, ntot) )
print( 'method = srg, generator = %s' % generator )

# Save evolved potential for each lambda value
for lamb in lambda_array:

    H_evolved = d[lamb] # Scattering units here [fm^-2]
    # Subtract off kinetic energy (need to convert T_rel from MeV to fm^-2)
    V_evolved = H_evolved - T_rel / hbar_sq_over_M
    # Save evolved potential
    lp.save_potential(k_array, k_weights, V_evolved, kvnn, channel, kmax, kmid,
                      ntot, 'srg', generator, lamb)