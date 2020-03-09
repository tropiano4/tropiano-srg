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


# Description of this test:
#   Testing out how to consistently divide out momenta and weights in operators
#   and wave functions.


#from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from Figures import figures_functions as ff
import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# TILDE MEANS THAT THE WAVE FUNCTION OR OPERATOR INCLUDES FACTORS OF K AND DK


# --- Set up --- #

# Specify potential and SRG-evolution here
kvnn = 111
channel = '3S1'
ntot = 120
method = 'srg'
generator = 'Wegner'
lamb = 1.5

# Load momentum and weights, and define factor_array
k_array, k_weights = lsp.load_momentum(kvnn, channel)
factor_array = np.sqrt(2/np.pi) * k_array * np.sqrt(k_weights)
factor_array_long = np.concatenate( (factor_array, factor_array) )
row, col = np.meshgrid(factor_array_long, factor_array_long)

# Load un-evolved and evolved Hamiltonian in units MeV
H_initial = lsp.load_hamiltonian(kvnn, channel)
H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg',
                                 generator=generator, lamb=lamb)
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)


# --- Wave function --- #

# Load un-evolved and evolved wave functions
psi_initial_tilde = ob.wave_function(H_initial)
psi_evolved_tilde = ob.wave_function(H_evolved)

# Convert to fm^3/2 and split in 3S1 and 3D1 components
psi_initial = psi_initial_tilde / factor_array_long
u_initial = psi_initial[:ntot]
w_initial = psi_initial[ntot:]
psi_evolved = psi_evolved_tilde / factor_array_long
u_evolved = psi_evolved[:ntot]
w_evolved = psi_evolved[ntot:]

# Calculate momentum distributions
phi_squared_initial = u_initial**2 + w_initial**2
phi_squared_evolved = u_evolved**2 + w_evolved**2
phi_squared_old = phi_squared_initial * 2/np.pi

# Plot wave functions
plt.close('all')
f, ax = plt.subplots(figsize=(4, 4)) 
    
ax.semilogy(k_array, phi_squared_evolved, color=ff.xkcd_colors(1),
            label='Wegner', linestyle='dashdot', linewidth=2.0)
ax.semilogy(k_array, phi_squared_initial, color=ff.xkcd_colors(0),
            label='Initial', linestyle='dotted', linewidth=2.0)
ax.semilogy(k_array, phi_squared_old, color=ff.xkcd_colors(2),
            label='Initial * 2/pi', linestyle='dotted', linewidth=2.0)
# Specify axes limits
ax.set_xlim([0.0, 4.0])
ax.set_ylim([1e-5, 1e3])
ax.legend()
    
# Set axes labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']')
ax.set_ylabel(r'$\phi_d^2$' + ' [fm' + r'$^3$' + ']')
plt.show()


# --- Momentum projection operator --- #

# Set momentum value
q = 3.00

# Load operator
momentum_proj_tilde = op.momentum_projection_operator(q, k_array, k_weights, 
                                                      channel, U=U_matrix)
momentum_proj = momentum_proj_tilde / row / col

# Interpolate
k_array_int, momentum_proj_int = ff.interpolate_matrix(k_array, momentum_proj,
                                                       4.0)

# Plot operators as contours
plt.close('all')
f, ax = plt.subplots(figsize=(4, 4))

# Settings for contour
mx = 0.01
mn = -0.01
levels_number = 61
levels = np.linspace(mn, mx, levels_number)
levels_ticks = np.linspace(mn, mx, 9)
levels_ticks_strings = ['%.3f' % tick for tick in levels_ticks]
    
c = ax.contourf(k_array_int, k_array_int, momentum_proj_int, levels, 
                cmap='jet', extend='both')

# Specify axes limits
ax.set_xlim( [0.0, 4.0] )
ax.set_ylim( [0.0, 4.0] )
# Switch from bottom to top
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.tick_params(labeltop=True)
                                         
# Set axes labels
ax.set_xlabel("k' [fm" + r'$^{-1}$' + ']')
ax.set_ylabel('k [fm' + r'$^{-1}$' + ']')
                                         
# Invert y-axis
plt.gca().invert_yaxis()
# Set colorbar axe
f.subplots_adjust(right=0.8) # Adjust for colorbar space
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
# Set colorbar
cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
cbar.ax.set_yticklabels(levels_ticks_strings)

plt.show()