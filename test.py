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
#------------------------------------------------------------------------------


import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
import observables as ob
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Set up --- #

# Potential details
kvnn_list = [10, 106]
channel = '3S1'
kmax = 10.0
kmid = 2.0
ntot = 120

# SRG details
generator = 'Wegner'
lamb = 1.5

# Deuteron energy in MeV
eps = -2.22


# --- Load wave functions, calculate momentum distribution, and store --- #

d = {}
d['init'] = {}
d['srg'] = {}
for kvnn in kvnn_list:
    
    # Load momentum and weights
    k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # Load initial Hamiltonian
    H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    
    # Compute initial wave function
    psi_initial = ob.wave_function(H_initial, eps)
    u_initial = psi_initial[:ntot] # 3S1 component
    w_initial = psi_initial[ntot:] # 3D1 component
    
    # Initial momentum distribution (divide by momenta and weights for 
    # mesh-independent result)
    d['init'][kvnn] = ( u_initial**2 + w_initial**2 ) / \
                         ( k_array**2 * k_weights )
    
    # Build unitary transformation
    H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg',
                                    generator, lamb)
    U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
        
    # Compute evolved wave function
    psi_evolved = ob.wave_function(H_initial, eps, U=U_matrix)
    u_evolved = psi_evolved[:ntot] # 3S1 component
    w_evolved = psi_evolved[ntot:] # 3D1 component
        
    # Evolved momentum distribution (divide by momenta and weights for 
    # mesh-independent result)
    d['srg'][kvnn] = ( u_evolved**2 + w_evolved**2 ) / \
                     ( k_array**2 * k_weights )
    
# Calculate ratio and absolute difference
diff_init = abs( d['init'][kvnn_list[0]] - d['init'][kvnn_list[1]] )
ratio_init = d['init'][kvnn_list[0]] / d['init'][kvnn_list[1]]             
diff_srg = abs( d['srg'][kvnn_list[0]] - d['srg'][kvnn_list[1]] )
ratio_srg = d['srg'][kvnn_list[0]] / d['srg'][kvnn_list[1]]


# --- Plot --- #

# Absolute difference
plt.clf()
plt.semilogy(k_array, diff_init, color='xkcd:blue', label='Initial')
plt.semilogy(k_array, diff_srg, color='xkcd:red', label='SRG')
plt.xlim( (0.0, 3.0) )
plt.ylim( (1e-9, 1e0) )
plt.title('Absolute difference in '+r'$|\phi_d(k)|^2$')
plt.xlabel('k [fm' + r'$^{-1}$' + ']')
plt.ylabel('EM - RKE [fm' + r'$^3$' + ']')
plt.legend()
plt.show()

# Ratio
plt.clf()
plt.plot(k_array, ratio_init, color='xkcd:blue', label='Initial')
plt.plot(k_array, ratio_srg, color='xkcd:red', label='SRG')
plt.plot(k_array, np.ones(ntot), 'k:')
plt.xlim( (0.0, 3.0) )
plt.ylim( (0.0, 3.0) )
#plt.xlim( (0.0, 10.0) )
#plt.ylim( (0.0, 10.0) )
plt.title('Ratio of '+r'$|\phi_d(k)|^2$')
plt.xlabel('k [fm' + r'$^{-1}$' + ']')
plt.ylabel('EM / RKE')
plt.legend()
plt.show()