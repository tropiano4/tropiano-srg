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
#                coordinate-space).
#
#------------------------------------------------------------------------------


from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Set up --- #

cwd = getcwd()

# Use RKE N3LO as test
#kvnn = 106
#kmax = 10.0
#kmid = 2.0

# Use EM N3LO as test
#kvnn = 10
#kmax = 30.0
#kmax = 10.0
#kmid = 4.0
#kmid = 2.0

# Use Gezerlis et al. as test
kvnn = 222
kmax = 10.0
kmid = 2.0

channel = '3S1'
ntot = 120

# SRG details
generator = 'Wegner'
lambda_array = np.array( (6.0, 3.0, 2.0, 1.5) )

# Load Hamiltonians and momentum
H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
# The arrays below are used later to present a mesh-independent result
factor_array = np.concatenate( (k_array * np.sqrt(k_weights), k_array * np.sqrt(k_weights) ) )
row, col = np.meshgrid(factor_array, factor_array)

# Initial operator
constant = 1.0
operator_initial = np.ones( (2*ntot, 2*ntot) ) * constant * row * col


# --- Calculate evolved operator --- #

d = {}

for lamb in lambda_array:
    
    # Load evolved Hamiltonian
    H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg',
                                    generator, lamb)

    # Load SRG unitary transformation
    U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
    
    # Calculate evolved operator
    evolved_operator = U_matrix @ operator_initial @ U_matrix.T
    
    # Re-size to same dimension of k_array and store in dictionary
    d[lamb] = ( evolved_operator / row / col )[:ntot, :ntot]


# --- Plot figure --- #
    
# Limits of axes on contours (units are fm^-1)
axes_max = 10.0
        
# Specifications of x and y axes
# Step-size in labeling tick marks
axes_stepsize = 2.0
# x and y axes ticks
axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
        
# Labels
axes_label = 'k [fm' + r'$^{-1}$' + ']'
lambda_label = r'$\lambda=%.1f$' + ' fm' + r'$^{-1}$'
generator_label = ff.generator_label_conversion(generator)
    
# Fontsize for labels and tick marks
axes_label_size = 18
lambda_label_size = 17
generator_label_size = 17
axes_tick_size = 18
colorbar_tick_size = 18
        
# Limits of colorbar
mx = 4.0
mn = 0.0
        
# Location of labels
generator_label_location = 'upper right'
lambda_label_location = 'lower right'
        
# Color scheme for contour plots
color_style = 'Blues'
    
# Size of figure
row_number = 1
col_number = len(lambda_array)
figure_size = (4*col_number, 3.5*row_number)

# Initialize figure
plt.close('all')
f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True,
                      figsize=figure_size)

i = 0 # Sub-plot number
for lamb in lambda_array:
    
    # Add sub-plot to figure (for last sub-plot, must specify colorbar c)
    if i != ( len(lambda_array) - 1):
        axs[i].pcolormesh(k_array, k_array, d[lamb], cmap=color_style,
                          vmin=mn, vmax=mx, rasterized=True)
    else:
        c = axs[i].pcolormesh(k_array, k_array, d[lamb], cmap=color_style,
                              vmin=mn, vmax=mx, rasterized=True)
    
    # Specify axes tick marks
    axs[i].xaxis.set_ticks(axes_ticks)
    axs[i].xaxis.set_ticklabels(axes_ticks)
        
    # Specify axes limits
    axs[i].set_xlim( (0, axes_max) )
    axs[i].set_ylim( (0, axes_max) )
        
    # Position of x-axis label and tick marks
    axs[i].xaxis.set_label_position('top')
    axs[i].xaxis.tick_top()
    axs[i].tick_params(labeltop=True, labelsize=axes_tick_size)
        
    # Prevent overlapping x-axis tick marks unless it's the last sub-plot
    if i != ( len(lambda_array) - 1 ):
        xticks = axs[i].xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
            
    # Set axes labels
    axs[i].set_xlabel(axes_label, fontsize=axes_label_size)
    # Only specify y axis tick marks, set label, and add generator label as anchored text for 1st sub-plot
    if i == 0:
        axs[i].yaxis.set_ticks(axes_ticks)
        axs[i].yaxis.set_ticklabels(axes_ticks)
        axs[i].set_ylabel(axes_label, fontsize=axes_label_size)
        generator_anchored_text = AnchoredText(generator_label, prop=dict(size=generator_label_size),
                                               loc=generator_label_location)
        axs[i].add_artist(generator_anchored_text)
            
    # Add lambda label as anchored text
    lambda_anchored_text = AnchoredText(lambda_label % lambda_array[i], 
                                        prop=dict(size=lambda_label_size),
                                        loc=lambda_label_location)
    axs[i].add_artist(lambda_anchored_text)
            
    i += 1
        
    
# --- Set figure specifications and save --- #

# Invert y-axis
plt.gca().invert_yaxis()
# Amount of white space in-between sub-plots
f.subplots_adjust(hspace=0.0, wspace=0.0)
# Adjust for colorbar space
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
# Add colorbar and set tick size
cbar = f.colorbar(c, cax=cbar_ax)
cbar.ax.tick_params(labelsize=colorbar_tick_size)

file_name = 'toy_operator_srg_evolution_kvnn%d_kmax%.1f_kmid%.1f_v1' % (kvnn, kmax, kmid)
    
# Save figure
chdir('Figures/SRG_operators')
f.savefig(file_name+'.pdf', bbox_inches='tight')
chdir(cwd)