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
#   Testing the deuteron momentum distribution using projection operators from
#   operator.py.
#
#------------------------------------------------------------------------------


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# Load momentum distributions

# Details of potential
kvnn = 10
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120

# SRG generator
generator = 'Wegner'
# SRG evolution value
lamb = 1.5

# Deuteron energy in MeV
eps = -2.22

    
# Load initial Hamiltonian, momentum, and weights
H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot) 
    
# Load unitary transformation
# SRG calls function which builds U(s) out of un-evolved and evolved
# eigenvectors
H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg',
                                generator, lamb)
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
        
# Compute initial wave functions
psi_initial = ob.wave_function(H_initial, eps)
u_initial = psi_initial[:ntot] # 3S1 component
w_initial = psi_initial[ntot:] # 3D1 component
    
# Initial and evolved momentum distribution (divide by momenta and weights for
# mesh-independent result)
phi_squared_initial = ( u_initial**2 + w_initial**2 ) / \
                      ( k_array**2 * k_weights )

# Load evolved wave function
psi_evolved = ob.wave_function(H_initial, eps, U_matrix) # Unitless

# Initialize evolved momentum distribution
phi_squared_evolved = np.zeros(ntot)

i = 0
for q in k_array:
    
    # Load evolved momentum projection operator
    momentum_proj_op = op.momentum_projection_operator(q, k_array, k_weights,
                                                       U_matrix) # Units fm^3
    
    phi_squared_evolved[i] = psi_evolved.T @ momentum_proj_op @ psi_evolved
    
    i += 1
    

# Plot specifications

# Limits of x and y axes
xlim = [0.0, 4.0]
ylim = [1e-5, 1e3]
    
# Labels
k_label = r'$k \/ [fm^{-1}]$'
phi_squared_label = r'$\phi_d^2 \/ \/ [fm^3]$'
lamb_label = r'$\lambda=%.1f \/ \/ fm^{-1}$'
generator_label = r'$G=H_{D}$'
    
# Fontsize for labels
legend_label_size = 16
k_label_size = 18
phi_squared_label_size = 20
generator_label_size = 20
    
# Location of labels
legend_label_location = 1
generator_label_location = 3
    
    
# Plot the semi-log figure
plt.close('all')

f, ax = plt.subplots()
    
ax.semilogy(k_array, phi_squared_evolved, 'r-', label=lamb_label%lamb)
ax.semilogy(k_array, phi_squared_initial, 'k:', 
            label=r'$\lambda=\infty \/ \/ fm^{-1}$')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.legend(loc=legend_label_location, frameon=False, fontsize=legend_label_size)
ax.set_xlabel(k_label, fontsize=k_label_size)
ax.set_ylabel(phi_squared_label, fontsize=phi_squared_label_size)
anchored_text = AnchoredText(generator_label, 
                             prop=dict(size=generator_label_size),
                             loc=generator_label_location, frameon=False)
ax.add_artist(anchored_text)

plt.show()