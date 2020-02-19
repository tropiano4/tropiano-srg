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
#
#------------------------------------------------------------------------------


import matplotlib.pyplot as plt
import numpy as np
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# Potential specifications
#kvnn = 79
#kvnn = 111
kvnn = 222
channel = '3S1'
kmax = 10.0
kmid = 2.0
ntot = 120

# SRG specifications
generator = 'Wegner'
#generator = 'Block-diag'
lambda_bd = 2.00
lamb = 2.0


# Load momentum
k_array, _ = lsp.load_momentum(kvnn, channel, kmax, kmid, ntot)
# Load initial and evolved Hamiltonians
H_initial = lsp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
H_evolved = lsp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg', 
                                 generator, lamb, lambda_bd)
# Unitary transformation
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)

# k_0, k_i values
k_0 = 0.1
k_0_index = op.find_q_index(k_0, k_array)
k_values = np.array([0.5, 1.0, 1.5, 3.0])

# Initialize dictonary to store ratios
d = {}

# Loop over k_i values and add to ratio arrays
for k_i in k_values:
    
    k_i_index = op.find_q_index(k_i, k_array)
    numerator_array = U_matrix[k_i_index, :ntot]
    denominator_array = U_matrix[k_0_index, :ntot]
    d[k_i] = abs( numerator_array / denominator_array )


# Plot figures
plt.plot(k_array, d[k_values[0]], label=r'$k_i=%.1f$'%k_values[0])
plt.plot(k_array, d[k_values[1]], label=r'$k_i=%.1f$'%k_values[1])
plt.plot(k_array, d[k_values[2]], label=r'$k_i=%.1f$'%k_values[2])
plt.plot(k_array, d[k_values[3]], label=r'$k_i=%.1f$'%k_values[3])
plt.xlim([0.0, 5.0])
plt.ylim([0.0, 60.0])
plt.xlabel('q')
plt.ylabel(r'$|U(k_i,q)/U(k_0,q)|$')
plt.legend(loc='upper right')
plt.show()