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
#   09/24/19 --- Comparing wave functions from different SRG-evolved potentials
#                by looking at momentum distribution functions.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   10/22/19 --- Comparing the initial and SRG block-diagonal evolved deuteron
#                wave functions squared. Possible connection between V_low-k
#                and block-diagonal SRG.
#   10/29/19 --- Testing r^2 operator and RMS half-radius of deuteron.
#
#------------------------------------------------------------------------------


import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
import observables as ob
import operators as op
#from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# Set-up
kvnn = 111 # RKE N4LO (450 MeV)
channel = '3S1'
kmax = 8.0
kmid = 2.0
ntot = 120

k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)

# Initial Hamiltonian and wave function
H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
psi_initial = ob.wave_function(H_initial)

# SRG-evolved
#generator = 'Wegner'
#lamb = 1.5
#H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg', 
                                #generator, lamb)
#U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
#psi_evolved = ob.wave_function(H_initial, U=U_matrix)
# Should work the same way
#psi_evolved = ob.wave_function(H_evolved)

# Loop over r-dependent stuff for various r_max values

r_min = 0.005
r_max_array = np.arange(10.0, 32.0, 2.0)
N = len(r_max_array)
dr = 0.005

rms_radius_array = np.zeros(N)

for i in range(N):
    
    r_max = r_max_array[i]
    r_array = np.arange(r_min, r_max + dr, dr)
    
    # r^2 operator
    r2_operator_init = op.r2_operator(k_array, k_weights, r_array, dr)
    #r2_operator_evolved = op.r2_operator(k_array, k_weights, r_array, dr,
                                        #U=U_matrix)
    
    rms_radius_array[i] = ob.rms_radius_from_rspace(psi_initial,
                                                    r2_operator_init, k_array,
                                                    k_weights)    
    
# Compare to exact value
rms_radius_exact_array = np.repeat(1.966, N)
    
# Plot figure
plt.plot(r_max_array, rms_radius_exact_array, 'r-')
plt.plot(r_max_array, rms_radius_array, 'k:o')
plt.ylabel(r'$r_d$'+' [fm]')
plt.xlabel(r'$r_{max}$'+' [fm]')
plt.show()
