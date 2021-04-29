#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a separate file with an
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
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#   01/26/21 --- Renamed to test_script.py.
#   03/16/21 --- Testing pair momentum distributions for different nuclei
#                using LDA with simple filled Fermi sea single-particle
#                momentum distributions. Created pmd_simple_test.py in
#                Old_codes.
#   03/16/21 --- Testing pair momentum distribution for deuteron starting from
#                second quantization derivation and using expansion of
#                U(k, k'). Created pmd_deuteron_test.py in Old_codes.
#   03/16/21 --- Testing single-nucleon momentum distributions for different
#                nuclei using LDA. Created single_particle_momentum_dist.py
#                in Old_codes.
#   04/14/21 --- Creating AV18 SRG evolution figure for APS April Meeting
#                presentation.
#                potential_contours_kvnn_6_channel_1P1_Wegner_lamb_1p5.png in
#                Figures/Operator_evolution/Old_figures/Potential_contours.
#   04/28/21 --- Testing normalization and contributions of \delta U, etc. or
#                pp/pn to single-nucleon momentun distributions. Created
#                lda_normalizations_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Looking at higher partial waves for potential bug in 3P2 and 3D3 channels.
#   Seems 3S1-3D1, 3P2-3F2, 3D3-3G3 channels accumulate large numbers at high
#   momentum in
#       < k | \delta U | k' > and < k | \delta U \delta U^\dagger | k' >.
#   Results:
#     1. Numerical artifact at end of AV18 \delta U(k,k) in 3S1 channel,
#        kmax=10 fm^-1 but doesn't seem to disrupt calculations. (Minor for
#        kmax=15 fm^-1.)
#     2. Numerical artifacts at end of AV18 \delta U(k,k) in 3P2, 3D3
#        channels, kmax=10, 15. (Minor for kmax=30.)
#     3. Numerical artifacts at front of AV18, N2LO \delta U(k,k) in 3P2, 3D3
#        channels, kmax=10, 15, 30 fm^-1 and screws up calculations.


import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


def delta_U_contours(q_array, matrix, axes_lim=(0.0, 5.0), colorbar_limits=(-0.1, 0.1)):
    
    # --- Set-up --- #
    
    # Size of figure
    row_number = 1
    col_number = 1
    figure_size = (4*col_number, 3.5*row_number) # Extra width for colorbar
    
    # Axes limits
    axes_max = axes_lim[1]
    
    # Axes ticks, labels, and fontsizes
    x_label = "k' [fm" + r'$^{-1}$' + ']'
    y_label = 'k [fm' + r'$^{-1}$' + ']'
    axes_label_size = 18
    if 15 < axes_max <= 30:
        axes_stepsize = 5.0 # Step-size in labeling tick marks
    elif 10 < axes_max <= 15:
        axes_stepsize = 3.0
    elif 5 < axes_max <= 10:
        axes_stepsize = 2.0
    else:
        axes_stepsize = 1.0
    axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
    axes_ticks_strings = ff.convert_ticks_to_labels(axes_ticks)
    axes_tick_size = 18
    
    # Colorbar ticks, label, and fontsize
    mn = colorbar_limits[0]
    mx = colorbar_limits[1]
    levels_number = 61
    levels = np.linspace(mn, mx, levels_number)
    levels_ticks = np.linspace(mn, mx, 9)
    levels_ticks_strings = ff.convert_ticks_to_labels(levels_ticks)
    colorbar_tick_size = 18

    # Color scheme for contour plots
    color_style = 'turbo'


    # --- Plot data --- #
    
    # Initialize figure
    plt.close('all')
    f, ax = plt.subplots(row_number, col_number, figsize=figure_size)
    
    c = ax.contourf(q_array, q_array, matrix, levels, cmap=color_style,
                    extend='both')
        
    # Specify axes limits
    ax.set_xlim( axes_lim )
    ax.set_ylim( axes_lim )
                     
    # Specify axes tick marks
    ax.xaxis.set_ticks(axes_ticks)
    ax.xaxis.set_ticklabels(axes_ticks_strings)
    # Switch from bottom to top
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=True, labelsize=axes_tick_size)

    # Set x-axis label
    ax.set_xlabel(x_label, fontsize=axes_label_size)
                                         
    # Specify axes tick marks
    ax.yaxis.set_ticks(axes_ticks)
    ax.yaxis.set_ticklabels(axes_ticks_strings)
    ax.tick_params(labelsize=axes_tick_size)
                                      
    # Set y-axis label
    ax.set_ylabel(y_label, fontsize=axes_label_size)

    # Invert y-axis
    plt.gca().invert_yaxis()
                                         
    # Set colorbar axes
    f.subplots_adjust(right=0.8) # Adjust for colorbar space
    cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
                                         
    # Set tick size and labels
    cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
    cbar.ax.tick_params(labelsize=colorbar_tick_size)
    cbar.ax.set_yticklabels(levels_ticks_strings)

    return f, ax


# Load \delta U term here for 3P2
kvnn = 6
# kvnn = 222
# channel = '3S1'
# channel = '3P2'
channel = '3D3'
lamb = 1.35
# kmax, kmid, ntot = 10.0, 2.0, 120
kmax, kmid, ntot = 15.0, 3.0, 120
# kmax, kmid, ntot = 30.0, 4.0, 120

# Load and save momentum arrays for integration
k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
# For dividing out momenta/weights
factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
# For coupled-channel matrices
factor_array_cc = np.concatenate( (factor_array, factor_array) )
        
        
# --- Evaluate matrix elements --- #
        
# Load SRG transformation
H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, method='srg',
                                 generator='Wegner', lamb=lamb)
# Load U(k, k') [unitless]
U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)

# Isolate 2-body term and convert to fm^3
I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
row, col = np.meshgrid(factor_array_cc, factor_array_cc)

delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
delta_U_matrix = delta_U_matrix_unitless / row / col
            
# 2J+1 factor
J = int( channel[-1] )
            
# Add to the pp and pn terms
# Coupled-channel
# First L of coupled-channel
# Isospin CG's=1/\sqrt(2) for pn
deltaU = (2*J+1)/2 * ( delta_U_matrix[:ntot, :ntot] + delta_U_matrix[ntot:, ntot:] )
deltaU2 = (2*J+1)/4 * ( delta_U_matrix[:ntot, :ntot]**2 + \
                        delta_U_matrix[:ntot, ntot:]**2 + \
                        delta_U_matrix[ntot:, :ntot]**2 + \
                        delta_U_matrix[ntot:, ntot:]**2 )

# --- Plot \delta U and \delta U^2 --- #

axes_lim = (0.0, kmax)

if channel == '3S1':
    c_lim = (-0.4, 0.4)
elif channel == '3P2':
    c_lim = (-0.2, 0.2)
elif channel == '3D3':
    c_lim = (-0.04, 0.04)

f, ax = delta_U_contours(k_array, deltaU, axes_lim, c_lim)

# Add matrix labels
matrix_label_size = 20
matrix_label_location = 'lower left'
matrix_label = r'$\delta$'+'U ' + ff.channel_label_conversion(channel)
anchored_text = AnchoredText(matrix_label, loc=matrix_label_location,
                             prop=dict(size=matrix_label_size))
ax.add_artist(anchored_text)

plt.show()


if channel == '3S1':
    c_lim = (-0.1, 0.1)
elif channel == '3P2':
    c_lim = (-0.02, 0.02)
elif channel == '3D3':
    c_lim = (-0.004, 0.004)

f, ax = delta_U_contours(k_array, deltaU2, axes_lim, c_lim)

# Add matrix labels
matrix_label_size = 20
matrix_label_location = 'lower left'
matrix_label = r'$\delta$'+'U'+r'$^2$' ' ' + ff.channel_label_conversion(channel)
anchored_text = AnchoredText(matrix_label, loc=matrix_label_location,
                             prop=dict(size=matrix_label_size))
ax.add_artist(anchored_text)

plt.show()


# --- Plot diagonal terms of \delta U --- #

plt.close('all')
f, ax = plt.subplots(1, 1, figsize=(4,4))

ax.plot(k_array, np.diag(deltaU))

ax.set_xlim( axes_lim )
if channel == '3S1':
    ax.set_ylim( (-0.5, 0.01) )
elif channel == '3P2':
    ax.set_ylim( (-0.5, 0.01) )
elif channel == '3D3':
    ax.set_ylim( (-0.2, 0.01) )

ax.set_xlabel('k ' + r'$fm^{-1}$')
ax.set_ylabel(r'$\delta$'+'U(k,k)')

# Add matrix labels
matrix_label_size = 16
matrix_label_location = 'lower center'
matrix_label = ff.channel_label_conversion(channel)
anchored_text = AnchoredText(matrix_label, loc=matrix_label_location,
                             prop=dict(size=matrix_label_size), frameon=False)
ax.add_artist(anchored_text)