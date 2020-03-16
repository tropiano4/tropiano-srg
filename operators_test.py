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


from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from Figures import figures_functions as ff
import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# TILDE MEANS THAT THE WAVE FUNCTION OR OPERATOR INCLUDES FACTORS OF K AND DK


# --- Set up --- #

cwd = getcwd()

# Specify potential and SRG-evolution here
kvnn = 10
channel = '3S1'
# Test two different meshes
kmax_1 = 10.0
kmid_1 = 2.0
kmax_2 = 30.0
kmid_2 = 4.0
ntot = 120
method = 'srg'
generator = 'Wegner'
lamb = 1.5

# Legend label
legend_label = r'$k_{\rm max}=%.1f$' + ' fm' + r'$^{-1}$'

# Load momentum and weights, and define factor_array
k_array_1, k_weights_1 = lsp.load_momentum(kvnn, channel, kmax=kmax_1, 
                                           kmid=kmid_1, ntot=ntot)
factor_array_1 = np.sqrt(2/np.pi) * k_array_1 * np.sqrt(k_weights_1)
factor_array_long_1 = np.concatenate( (factor_array_1, factor_array_1) )
row_1, col_1 = np.meshgrid(factor_array_long_1, factor_array_long_1)
k_array_2, k_weights_2 = lsp.load_momentum(kvnn, channel, kmax=kmax_2, 
                                           kmid=kmid_2, ntot=ntot)
factor_array_2 = np.sqrt(2/np.pi) * k_array_2 * np.sqrt(k_weights_2)
factor_array_long_2 = np.concatenate( (factor_array_2, factor_array_2) )
row_2, col_2 = np.meshgrid(factor_array_long_2, factor_array_long_2)

# Specify coordinates array
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

# Load Hamiltonians and unitary transformations
H_initial_1 = lsp.load_hamiltonian(kvnn, channel, kmax=kmax_1, kmid=kmid_1, 
                                   ntot=ntot)
H_evolved_1 = lsp.load_hamiltonian(kvnn, channel, kmax=kmax_1, kmid=kmid_1,
                                   ntot=ntot, method='srg', 
                                   generator=generator, lamb=lamb)
U_matrix_1 = SRG_unitary_transformation(H_initial_1, H_evolved_1)
H_initial_2 = lsp.load_hamiltonian(kvnn, channel, kmax=kmax_2, kmid=kmid_2, 
                                   ntot=ntot)
H_evolved_2 = lsp.load_hamiltonian(kvnn, channel, kmax=kmax_2, kmid=kmid_2,
                                   ntot=ntot, method='srg',
                                   generator=generator, lamb=lamb)
U_matrix_2 = SRG_unitary_transformation(H_initial_2, H_evolved_2)

# Load wave functions
psi_tilde_1 = ob.wave_function(H_initial_1)
psi_tilde_2 = ob.wave_function(H_initial_2)

# Load evolved momentum projection operators
q = 3.00
momentum_proj_op_tilde_1 = op.momentum_projection_operator(q, k_array_1, 
                                                           k_weights_1, 
                                                           channel, U=U_matrix_1)
momentum_proj_op_tilde_2 = op.momentum_projection_operator(q, k_array_2, 
                                                           k_weights_2,
                                                           channel, U=U_matrix_2)

# Load evolved r^2 operators
r2_op_tilde_1 = op.r2_operator(k_array_1, k_weights_1, r_array, dr, U=U_matrix_1)
r2_op_tilde_2 = op.r2_operator(k_array_2, k_weights_2, r_array, dr, U=U_matrix_2)


# --- Check normalization and plot momentum distributions --- #

normalization_1 = psi_tilde_1.T @ psi_tilde_1
print('kmax=%.1f normalization = %.3f' % (kmax_1, normalization_1))
normalization_2 = psi_tilde_2.T @ psi_tilde_2
print('kmax=%.1f normalization = %.3f' % (kmax_2, normalization_2))

# Convert to fm^3/2 and split in 3S1 and 3D1 components
psi_1 = psi_tilde_1 / factor_array_long_1
phi_squared_1 = psi_1[:ntot]**2 + psi_1[ntot:]**2
psi_2 = psi_tilde_2 / factor_array_long_2
phi_squared_2 = psi_2[:ntot]**2 + psi_2[ntot:]**2

# Plot momentum distributions
plt.close('all')
f, ax = plt.subplots(figsize=(4, 4))
    
ax.semilogy(k_array_1, phi_squared_1, color=ff.xkcd_colors(0),
            label=legend_label % kmax_1, linestyle='solid', linewidth=2.0)
ax.semilogy(k_array_2, phi_squared_2, color=ff.xkcd_colors(1),
            label=legend_label % kmax_2, linestyle='dashdot', linewidth=2.0)

# Specify axes limits
ax.set_xlim([0.0, 4.0])
ax.set_ylim([1e-5, 1e3])
ax.legend(loc='upper right', fontsize=13, frameon=False)

# Set axes labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=18)
ax.set_ylabel(r'$\phi_d^2$' + ' [fm' + r'$^3$' + ']', fontsize=20)

# Enlarge axes tick marks
ax.tick_params(labelsize=14)

chdir('Figures/Operator_evolution')
f.savefig('momentum_distribution_test.png', bbox_inches='tight')
chdir(cwd)
plt.show()


# --- Momentum projection operator --- #

# Convert momentum projection operator to units fm^6
# Load operator
momentum_proj_op_1 = ( momentum_proj_op_tilde_1 / row_1 / col_1 )[:ntot, :ntot]
momentum_proj_op_2 = ( momentum_proj_op_tilde_2 / row_2 / col_2 )[:ntot, :ntot]

# Interpolate
k_array_int_1, momentum_proj_int_1 = ff.interpolate_matrix(k_array_1, 
                                                           momentum_proj_op_1,
                                                           4.0)
k_array_int_2, momentum_proj_int_2 = ff.interpolate_matrix(k_array_2, 
                                                           momentum_proj_op_2,
                                                           4.0)

# Axes ticks
axes_max = 4.0
axes_stepsize = 1.0 # Step-size in labeling tick marks
axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)

# Plot operators as contours
plt.close('all')
f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True,
                             figsize=(8, 3.5))

# Settings for contour
mx = 0.01
mn = -0.01
levels_number = 61
levels = np.linspace(mn, mx, levels_number)
levels_ticks = np.linspace(mn, mx, 9)
levels_ticks_strings = ['%.3f' % tick for tick in levels_ticks]
    
ax1.contourf(k_array_int_1, k_array_int_1, momentum_proj_int_1, levels, 
             cmap='jet', extend='both')
c = ax2.contourf(k_array_int_2, k_array_int_2, momentum_proj_int_2, levels, 
                 cmap='jet', extend='both')

# Specify axes limits
ax1.set_xlim( [0.0, axes_max] )
ax1.set_ylim( [0.0, axes_max] )
ax2.set_xlim( [0.0, axes_max] )
ax2.set_ylim( [0.0, axes_max] )
# Switch from bottom to top
ax1.xaxis.set_label_position('top')
ax1.xaxis.tick_top()
ax1.xaxis.set_ticks(axes_ticks)
ax1.xaxis.set_ticklabels(axes_ticks)
ax1.tick_params(labeltop=True, labelsize=14)
ax1.yaxis.set_ticks(axes_ticks)
ax1.yaxis.set_ticklabels(axes_ticks)
ax2.xaxis.set_label_position('top')
ax2.xaxis.tick_top()
ax2.xaxis.set_ticks(axes_ticks)
ax2.xaxis.set_ticklabels(axes_ticks)
ax2.tick_params(labeltop=True, labelsize=14)
                                         
# Set axes labels
ax1.set_xlabel("k' [fm" + r'$^{-1}$' + ']', fontsize=18)
ax1.set_ylabel('k [fm' + r'$^{-1}$' + ']', fontsize=18)
ax2.set_xlabel("k' [fm" + r'$^{-1}$' + ']', fontsize=18)

# Add kmax labels
anchored_text_1 = AnchoredText(legend_label % kmax_1, loc='lower left', 
                               prop=dict(size=13))
ax1.add_artist(anchored_text_1)
anchored_text_2 = AnchoredText(legend_label % kmax_2, loc='lower left', 
                               prop=dict(size=13))
ax2.add_artist(anchored_text_2)
                                         
# Invert y-axis
plt.gca().invert_yaxis()
# Set colorbar axe
f.subplots_adjust(right=0.8) # Adjust for colorbar space
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
# Set colorbar
cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
cbar.ax.set_yticklabels(levels_ticks_strings)
# Set colorbar label
cbar.ax.set_title('[fm' + r'$^6$' + ']', fontsize=14)

chdir('Figures/Operator_evolution')
f.savefig('momentum_projection_operator_test.png', bbox_inches='tight')
chdir(cwd)
plt.show()


# --- r^2 operator --- #

# Convert r^2 operator to units fm^5
# Load operator
r2_op_1 = ( r2_op_tilde_1 / row_1 / col_1 )[:ntot, :ntot]
r2_op_2 = ( r2_op_tilde_2 / row_2 / col_2 )[:ntot, :ntot]

# Interpolate
k_array_int_1, r2_op_int_1 = ff.interpolate_matrix(k_array_1, r2_op_1, 0.4)
k_array_int_2, r2_op_int_2 = ff.interpolate_matrix(k_array_2, r2_op_2, 0.4)

# Axes ticks
axes_max = 0.4
axes_stepsize = 0.1 # Step-size in labeling tick marks
axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)

# Plot operators as contours
plt.close('all')
f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True,
                             figsize=(8, 3.5))

# Settings for contour
mx = 6e6
mn = -6e6
levels_number = 61
levels = np.linspace(mn, mx, levels_number)
levels_ticks = np.linspace(mn, mx, 9)
levels_ticks_strings = ['%.1e' % tick for tick in levels_ticks]
    
ax1.contourf(k_array_int_1, k_array_int_1, r2_op_int_1, levels, cmap='jet', 
             extend='both')
c = ax2.contourf(k_array_int_2, k_array_int_2, r2_op_int_2, levels, cmap='jet', 
                 extend='both')

# Specify axes limits
ax1.set_xlim( [0.0, axes_max] )
ax1.set_ylim( [0.0, axes_max] )
ax2.set_xlim( [0.0, axes_max] )
ax2.set_ylim( [0.0, axes_max] )
# Switch from bottom to top
ax1.xaxis.set_label_position('top')
ax1.xaxis.tick_top()
# ax1.xaxis.set_ticks(axes_ticks)
# ax1.xaxis.set_ticklabels(axes_ticks)
ax1.tick_params(labeltop=True, labelsize=14)
# ax1.yaxis.set_ticks(axes_ticks)
# ax1.yaxis.set_ticklabels(axes_ticks)
ax2.xaxis.set_label_position('top')
ax2.xaxis.tick_top()
# ax2.xaxis.set_ticks(axes_ticks)
# ax2.xaxis.set_ticklabels(axes_ticks)
ax2.tick_params(labeltop=True, labelsize=14)
                                         
# Set axes labels
ax1.set_xlabel("k' [fm" + r'$^{-1}$' + ']', fontsize=18)
ax1.set_ylabel('k [fm' + r'$^{-1}$' + ']', fontsize=18)
ax2.set_xlabel("k' [fm" + r'$^{-1}$' + ']', fontsize=18)

# Add kmax labels
anchored_text_1 = AnchoredText(legend_label % kmax_1, loc='lower right', 
                               prop=dict(size=13))
ax1.add_artist(anchored_text_1)
anchored_text_2 = AnchoredText(legend_label % kmax_2, loc='lower right', 
                               prop=dict(size=13))
ax2.add_artist(anchored_text_2)
                                         
# Invert y-axis
plt.gca().invert_yaxis()
# Set colorbar axe
f.subplots_adjust(right=0.8) # Adjust for colorbar space
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
# Set colorbar
cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
cbar.ax.set_yticklabels(levels_ticks_strings)
# Set colorbar label
cbar.ax.set_title('[fm' + r'$^5$' + ']', fontsize=14)

chdir('Figures/Operator_evolution')
f.savefig('r2_operator_test.png', bbox_inches='tight')
chdir(cwd)
plt.show()