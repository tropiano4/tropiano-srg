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
#   Compute <p^2> for the spurious state for kvnn=901


import numpy as np
import observables as ob
from Potentials.vsrg_macos import load_save_potentials as lsp


# Potential specifications
kvnn = 901
channel = '3S1'

# Load momentum
k_array, _ = lsp.load_momentum(kvnn, channel)
k_array_long = np.concatenate( (k_array, k_array) )

# Build p^2 operator in momentum-space
p2_operator = np.diag( k_array_long**2 )

# Spurious bound state energy
#eps = -2000
eps = -2.22

# Load initial Hamiltonian
H_matrix = lsp.load_hamiltonian(kvnn, channel)

# Spurious bound state wave function (unitless)
psi_s = ob.wave_function(H_matrix, eps)

p2_value = psi_s.T @ p2_operator @ psi_s # Units should be fm^-2
print(p2_value)