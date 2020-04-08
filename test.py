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
#   Generating a few figures relevant to r^2 operator evolution to understand
#   how r^2 evolves for different potentials and SRG generators.


import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from Figures import figures_functions as ff
import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Set up --- #

# Specify potential and SRG-evolution here
# RKE N4LO 450 MeV potential
kvnn = 111
# AV18
#kvnn = 6
channel = '3S1'
# Block-diagonal evolution
generator = 'Block-diag'
if kvnn == 111:
    lamb = 1.0
elif kvnn == 6:
    lamb = 1.5
lambda_bd = 2.0

# Plotting specifications that are dependent on the settings above
# Momentum settings
k_label = 'k [fm' + r'$^{-1}$' + ']'
kp_label = "k' [fm" + r'$^{-1}$' + ']'
k_max_1 = 0.4
k_max_2 = 4.0


# --- Main calculations --- #

# Load momentum and weights arrays
k_array, k_weights = lsp.load_momentum(kvnn, channel)
factor_array = np.concatenate( ( k_array * np.sqrt(k_weights), 
                                  k_array * np.sqrt(k_weights) ) ) * \
               np.sqrt(2/np.pi)
row, col = np.meshgrid(factor_array, factor_array)
ntot = len(k_array)

# Specify coordinates array
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

# Get initial and evolved Hamiltonians
H_initial = lsp.load_hamiltonian(kvnn, channel)
H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg',
                                 generator=generator, lamb=lamb, 
                                 lambda_bd=lambda_bd)

# Calculate unitary transformation
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)

# Load initial and evolved wave functions and divide out momenta/weights
psi_initial = ob.wave_function(H_initial) / factor_array
psi_evolved = ob.wave_function(H_initial, U=U_matrix) / factor_array

# Calculate momentum distributions
phi2_initial = psi_initial[:ntot]**2 + psi_initial[ntot:]**2
phi2_evolved = psi_evolved[:ntot]**2 + psi_evolved[ntot:]**2

# Print relative difference in (evolved - initial) / initial (percent) at
# k = 0 fm^-1
rel_diff = abs(phi2_evolved[0] - phi2_initial[0]) / phi2_initial[0] * 100.0
print(rel_diff)

# Calculate r^2 operators
r2_initial = op.r2_operator(k_array, k_weights, r_array, dr) / row / col
r2_evolved = op.r2_operator(k_array, k_weights, r_array, dr, U=U_matrix) \
             / row / col
        
# Take r^2 difference and re-size to 3S1 - 3S1 sub-block
r2_diff = (r2_evolved - r2_initial)[:ntot, :ntot]

# Matrix elements of <\psi|r^2|\psi>
psi_row, psi_col = np.meshgrid(psi_evolved, psi_evolved)
integrand = abs( psi_row * r2_evolved * psi_col )


# --- Plot momentum distributions --- #

# First plot momentum distributions on semi-log scale
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))

initial_label = 'Initial'
evolved_label = r'$\Lambda_{BD}=%.1f$' % lambda_bd + ' fm' + r'$^{-1}$'

ax.semilogy(k_array, phi2_evolved, color='xkcd:red', label=evolved_label, 
            linestyle='solid', linewidth=2.0)
ax.semilogy(k_array, phi2_initial, color='xkcd:black', label=initial_label, 
            linestyle='dotted', linewidth=2.0)

# Label axes and set legend
ax.set_xlabel(k_label)
ax.set_ylabel(r'$\phi_d^2(k)$' + ' [fm' + r'$^3$' + ']')
ax.legend(loc='upper right', frameon=False)

# Set axes limits
ax.set_xlim([0.0, 4.0])
ax.set_ylim([1e-5, 1e3])

# Save figure
file_name = 'momentum_distributions_test.pdf'
f.savefig(file_name, bbox_inches='tight')
plt.show()

# Now plot difference of evolved - initial
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))

ax.plot(k_array, phi2_evolved - phi2_initial, color='xkcd:blue', 
        linestyle='solid', linewidth=2.0)

# Label axes and set legend
ax.set_xlabel(k_label)
ax.set_ylabel(r'$\Delta \phi_d^2(k)$' + ' [fm' + r'$^3$' + ']')

# Set axes limits
ax.set_xlim([0.0, 0.4])
if kvnn == 111:
    ax.set_ylim([0.0, 0.8])
elif kvnn == 6:
    ax.set_ylim([0.0, 1.5])

# Save figure
file_name = 'momentum_distributions_difference_test.pdf'
f.savefig(file_name, bbox_inches='tight')
plt.show()


# --- Plot r^2 operator --- #

for k_max in [k_max_1, k_max_2]:

    levels_number = 61
    if k_max < 3.0:
        mx = 8e2
        mn = -8e2
    else:
        mx = 5e1
        mn = -5e1
        
    levels = np.linspace(mn, mx, levels_number)
    levels_ticks = np.linspace(mn, mx, 9)
    
    if k_max < 3.0:
        levels_ticks_strings = ['%.0f' % tick for tick in levels_ticks]
    else:
        levels_ticks_strings = ['%.1f' % tick for tick in levels_ticks]

    # Interpolate for better looking figure
    k_array_int, r2_diff_int = ff.interpolate_matrix(k_array, r2_diff, k_max)

    # Plot the difference of the r^2 operator
    plt.close('all')
    f, ax = plt.subplots(figsize=(4, 3.5))
    
    c = ax.contourf(k_array_int, k_array_int, r2_diff_int, levels, 
                    cmap='turbo', extend='both')

    # Set axes label
    ax.set_xlabel(kp_label)
    ax.set_ylabel(k_label)
                        
    # Specify axes limits
    ax.set_xlim( [0.0, k_max] )
    ax.set_ylim( [0.0, k_max] )
                     
    # Switch from bottom to top
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=True)

    # Invert y-axis
    plt.gca().invert_yaxis()
                                                                        
    # Set colorbar axe
    f.subplots_adjust(right=0.8) # Adjust for colorbar space
    cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
                                         
    # Set colorbar
    cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
    cbar.ax.set_yticklabels(levels_ticks_strings)

    # Save figure
    file_name = 'r2_operator_difference_test_kmax_%.1f.pdf' % k_max
    f.savefig(file_name, bbox_inches='tight')
    plt.show()


# --- Plot matrix elements of expectation value --- #

mx = 8
mn = -3
levels_number = 61
levels = np.logspace(mn, mx, levels_number)
levels_ticks = np.logspace(mn, mx, 12)
levels_ticks_strings = [r'$10^{%d}$' % step for step in range(mn, mx+1)]
colorbar_norm = colors.LogNorm(vmin=mn, vmax=mx)

# Interpolate for better looking figure
k_array_int, integrand_int = ff.interpolate_matrix(k_array, integrand, k_max)

# Plot the difference of the r^2 operator
plt.close('all')
f, ax = plt.subplots(figsize=(4, 3.5))
    
c = ax.contourf(k_array_int, k_array_int, integrand_int, levels, cmap='Blues', 
                norm=colors.LogNorm(), extend='both')

# Set axes label
ax.set_xlabel(kp_label)
ax.set_ylabel(k_label)
                        
# Specify axes limits
ax.set_xlim( [0.0, k_max] )
ax.set_ylim( [0.0, k_max] )
                     
# Switch from bottom to top
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.tick_params(labeltop=True)

# Invert y-axis
plt.gca().invert_yaxis()
                                                                        
# Set colorbar axe
f.subplots_adjust(right=0.8) # Adjust for colorbar space
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
                                         
# Set colorbar
cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
cbar.ax.set_yticklabels(levels_ticks_strings)

# Save figure
file_name = 'r2_integrand_test_kmax_%.1f.pdf' % k_max
f.savefig(file_name, bbox_inches='tight')
plt.show()