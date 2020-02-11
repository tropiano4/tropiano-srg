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
#   10/29/19 --- Testing r^2 operator and RMS half-radius of deuteron.
#   01/03/20 --- Looking at block-diagonal unitary transformations.
#   02/11/20 --- Trying to understand SRG induced terms in momentum projection
#                operator by plotting some terms.
#
#------------------------------------------------------------------------------


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp
import operators as op


# --- Set-up --- #

# Potential specifications
kvnn = 111
channel = '3S1'
kmax = 10.0
kmid = 2.0
ntot = 120

# SRG specifications
generator = 'Wegner'
lamb = 1.5

# Load momentum and potential
k_array, _ = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
V_matrix = lp.load_potential(kvnn, channel, kmax, kmid, ntot, 'srg', 
                             generator, lamb)

# --- Main calculation --- #

# Set q-value and take slices of potential correspondin
#q = 0.3
q = 3.0
q_index = op.find_q_index(q, k_array)

# Calculate induced contribution
if q == 3.0:
    row, col = np.meshgrid( V_matrix[:ntot, q_index], V_matrix[q_index, :ntot] )
    induced_cont = 2 * q**4 * row * col
else:
    row, col = np.meshgrid( k_array**2 * V_matrix[:ntot, q_index], 
                            k_array**2 * V_matrix[q_index, :ntot] )
    induced_cont = 2 * row * col

# Interpolate
k_array_int, induced_cont_int = ff.interpolate_matrix(k_array, induced_cont, 
                                                      4.0)


# --- Plot figure --- #

# Size of figure
row_number = 1
col_number = 1
figure_size = (4*col_number, 3.5*row_number) # extra width for colorbar

# Limits of x and y axes
axes_max = 4.0
        
# Specifications of axes
# Step-size in labeling tick marks
axes_stepsize = 1.0
# x and y axes ticks
axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
        
# Labels and fontsize
x_label = "k' [fm" + r'$^{-1}$' + ']'
y_label = 'k [fm' + r'$^{-1}$' + ']'
axes_label_size = 18
anchored_text_label = r'$\lambda=%.1f$' % lamb + '\nq=%.1f' % q
anchored_text_size = 17
axes_tick_size = 18
colorbar_tick_size = 18
    
# Color scheme for contour plots
color_style = 'jet'

# Location of labels


# Things which depend on the q value: limits of colorbar and label location
if q < axes_max/2:
    mx = 5.0
    mn = -5.0
    anchored_text_location = 'lower right'
else:
    mx = 1.2
    mn = -1.2
    anchored_text_location = 'upper left'
    
levels = np.linspace(mn, mx, 41)
    

plt.close('all')
f, ax = plt.subplots(row_number, col_number, figsize=figure_size)

c = ax.contourf(k_array_int, k_array_int, induced_cont_int, levels, 
                cmap=color_style)
# Specify axes limits
ax.set_xlim( (0, axes_max) )
ax.set_ylim( (0, axes_max) )
# Specify axes tick marks
ax.xaxis.set_ticks(axes_ticks)
ax.xaxis.set_ticklabels(axes_ticks)
ax.yaxis.set_ticks(axes_ticks)
ax.yaxis.set_ticklabels(axes_ticks)
# Position of x-axis label and tick marks
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.tick_params(labeltop=True, labelsize=axes_tick_size)
# Set axes labels
ax.set_xlabel(x_label, fontsize=axes_label_size)
ax.set_ylabel(y_label, fontsize=axes_label_size)
# Add lambda label as anchored text
anchored_text = AnchoredText(anchored_text_label, 
                             prop=dict(size=anchored_text_size),
                             loc=anchored_text_location)
ax.add_artist(anchored_text)

# Invert y-axis
plt.gca().invert_yaxis()
# Adjust for colorbar space
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
# Add colorbar and set tick size
cbar = f.colorbar(c, cax=cbar_ax)
cbar.ax.tick_params(labelsize=colorbar_tick_size)

plt.show()