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


from os import getcwd, chdir
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
import operators as op
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# Potential specifications
#kvnn = 6
kvnn = 901
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120

# SRG spefications
#generator = 'Wegner'
generator = 'Block-diag'
lamb = 1.2
lambda_bd = 2.00

# Load initial Hamiltonian, momentum, and weights
H_initial = lsp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
k_array, k_weights = lsp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
# Specify r_array
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

# Limits of axes on contours (units are fm^-1)
axes_max = 5.0

# Specifications of x and y axes
# Step-size in labeling tick marks
axes_stepsize = 1.0
# x and y axes ticks
axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
        
# Labels
x_label = "k' [fm" + r'$^{-1}$' + ']'
y_label = 'k [fm' + r'$^{-1}$' + ']'
# Label the block-diagonal Lambda
if generator == 'Block-diag':
    lambda_label = r'$\Lambda_{BD}=%.1f$' % lambda_bd + ' fm' + r'$^{-1}$'
# For band-diagonal generators, label lambda
else:
    lambda_label = r'$\lambda=%.1f$' % lamb + ' fm' + r'$^{-1}$'
generator_label = ff.generator_label_conversion(generator)
    
# Fontsize for labels and tick marks
axes_label_size = 18
lambda_label_size = 17
generator_label_size = 18
axes_tick_size = 18
colorbar_tick_size = 18
        
# Limits of colorbar (units are fm^6)
#mx = 1e6
#mn = -1e6
mx = 1e2
mn = -1e2
        
# Location of labels
generator_label_location = 'upper right'
lambda_label_location = 'lower left'
        
# Color scheme for contour plots
color_style = 'jet'
    
# Size of figure
figure_size = (4, 3.5) # (width, height) - extra width for colorbar
    
# File name of figure
file_name = 'r2_contour_kvnn%d_%s' % (kvnn, generator)

# Current working directory
cwd = getcwd()
    
    
# --- Load data and plot contours --- #

# Load unitary transformation
# SRG calls function which builds U(s) out of un-evolved and evolved eigenvectors
if generator == 'Block-diag':
    H_evolved = lsp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg', 
                                    generator, lamb, lambda_bd)
else:
    H_evolved = lsp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg', 
                                    generator, lamb)
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
    
initial_operator = op.r2_operator(k_array, k_weights, r_array, dr)[:ntot, :ntot]
evolved_operator = op.r2_operator(k_array, k_weights, r_array, dr, 
                                  U_matrix)[:ntot, :ntot]
operator_diff = (evolved_operator - initial_operator) / np.pi
        
# Interpolate the operator through 0 to axes_max for smoother looking figure 
# (the extension _int means interpolated)
#k_array_int, operator_int = ff.interpolate_matrix(k_array, operator, axes_max)
k_array_int, operator_int = ff.interpolate_matrix(k_array, operator_diff, axes_max)
        
        # Initialize figure
plt.close('all')
f, ax = plt.subplots(figsize=figure_size)
    
c = ax.pcolormesh(k_array_int, k_array_int, operator_int, cmap=color_style, 
                  vmin=mn, vmax=mx, rasterized=True)
# Set axes ticks
ax.xaxis.set_ticks(axes_ticks)
ax.yaxis.set_ticks(axes_ticks)
# Specify axes limits
ax.set_xlim( (0, axes_max) )
ax.set_ylim( (0, axes_max) )
# Position of x-axis label and tick marks
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.tick_params(labeltop=True, labelsize=axes_tick_size)
# Set axes labels
ax.set_xlabel(x_label, fontsize=axes_label_size)
ax.set_ylabel(y_label, fontsize=axes_label_size)
generator_anchored_text = AnchoredText(generator_label, 
                                       prop=dict(size=generator_label_size),
                                       loc=generator_label_location)
ax.add_artist(generator_anchored_text)
# Add lambda label as anchored text
lambda_anchored_text = AnchoredText(lambda_label, 
                                    prop=dict(size=lambda_label_size),
                                    loc=lambda_label_location)
ax.add_artist(lambda_anchored_text)
            

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
plt.show()

chdir('Figures/Operator_evolution')
f.savefig(file_name+'.pdf', bbox_inches='tight')
chdir(cwd)