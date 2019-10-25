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
#
#------------------------------------------------------------------------------


# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
import observables as ob
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


kvnn = 111 # RKE N4LO (450 MeV)
channel = '3S1'
kmax = 8.0
kmid = 2.0
ntot = 120

generator = 'Block-diag'
lamb = 1.5
#lamb = 1.0
#lamb = 3.0
#lamb = 6.0
#Lamb = 2.0
Lamb = 3.0
#Lamb = 1.5

#eps = -2.22
eps = 200.0

# Load initial Hamiltonian (and momentum)
k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
H_init = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)

# Calculate initial wave function squared
psi_init = ob.wave_function(H_init, eps)
u_init = psi_init[:ntot]
w_init = psi_init[ntot:]
phi_squared_init = ( u_init**2 + w_init**2 ) / ( k_array**2 * k_weights )

# Evolved
H_evol = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg', generator,
                             lamb, Lamb)
U_matrix = SRG_unitary_transformation(H_init, H_evol)

psi_evol = ob.wave_function(H_init, eps, U_matrix)
u_evol = psi_evol[:ntot]
w_evol = psi_evol[ntot:]
phi_squared_evol = ( u_evol**2 + w_evol**2 ) / ( k_array**2 * k_weights )

# Print
for i, j, k in zip(k_array, phi_squared_init, phi_squared_evol):
    
    print('k = %.5f fm^-1,   ratio = %.8f' % (i, k/j))
