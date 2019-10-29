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
#
#------------------------------------------------------------------------------


from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp
import observables as ob
import operators as op
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


cwd = getcwd()

# SET UP
kvnn = 10
channel = '3S1'
#kmax = 30.0
kmax = 10.0
#kmid = 4.0
kmid = 2.0
ntot = 120

#generator = 'T'
generator = 'Wegner'
lambda_array = np.array( (6.0, 3.0, 2.0, 1.5) )


# Load initial Hamiltonian, momentum, and weights
H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
# The arrays below are used later to present a mesh-independent result
#factor_array = k_array * np.sqrt(k_weights)
factor_array = k_array * np.sqrt(k_weights) * np.sqrt(2/np.pi)
#row, col = np.meshgrid(factor_array, factor_array)
row, col = np.meshgrid(factor_array**2, factor_array**2)
    
# Specify r_array
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

# Limits of axes on contours (units are fm^-1)
axes_max = 0.4
# Specifications of x and y axes
# Step-size in labeling tick marks
axes_stepsize = 0.1
# x and y axes ticks
axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
        
# Labels
x_label = "k' [fm" + r'$^{-1}$' + ']'
y_label = 'k [fm' + r'$^{-1}$' + ']'
# Label the block-diagonal Lambda
if generator == 'Block-diag':
    lambda_label = r'$\Lambda=%.1f$' + ' fm' + r'$^{-1}$'
# For band-diagonal generators, label lambda
else:
    lambda_label = r'$\lambda=%.1f$' + ' fm' + r'$^{-1}$'
generator_label = ff.generator_label_conversion(generator)
    
# Fontsize for labels and tick marks
axes_label_size = 18
lambda_label_size = 17
generator_label_size = 17
axes_tick_size = 18
colorbar_tick_size = 18
        
# Limits of colorbar (units are fm^6)
mx = 1e6
mn = -1e6
        
# Location of labels
generator_label_location = 'center right'
lambda_label_location = 'lower right'
        
# Color scheme for contour plots
color_style = 'jet'
    
# Size of figure
row_number = 1
col_number = len(lambda_array)
figure_size = (4*col_number, 3.5*row_number)


# Initialize figure
plt.close('all')
f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True, figsize=figure_size)
    
# Loop over lambda values keeping track of sub-plot number i
i = 0
for lamb in lambda_array:
        
    # Load unitary transformation
    # SRG calls function which builds U(s) out of un-evolved and evolved eigenvectors
    if generator == 'Block-diag':
        H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg',
                                        generator, 1.0, lambda_bd=lamb)
    else:
        H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg',
                                        generator, lamb)
    U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
        
    # Evolved momentum projection operator
    operator = op.r2_operator(k_array, k_weights, r_array, dr, U_matrix)
    operator = operator[:ntot, :ntot] / np.pi
    # Divide by k_i * k_j * sqrt( w_i * w_j ) for mesh-independent result
    #operator = operator / row / col
    
    # Interpolate the operator
    k_array_int, operator_int = ff.interpolate_matrix(k_array, operator, axes_max)
    
    # Add sub-plot to figure (for last sub-plot, must specify colorbar c)
    if i != ( len(lambda_array) - 1):
        axs[i].pcolormesh(k_array_int, k_array_int, operator_int, 
                          cmap=color_style, vmin=mn, vmax=mx, rasterized=True)
    else:
        c = axs[i].pcolormesh(k_array_int, k_array_int, operator_int, 
                              cmap=color_style, vmin=mn, vmax=mx, rasterized=True)
    
    # Specify axes tick marks
    axs[i].xaxis.set_ticks(axes_ticks)
        
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
    axs[i].set_xlabel(x_label, fontsize=axes_label_size)
    # Only specify y axis tick marks, set label, and add generator label as anchored text for 1st sub-plot
    if i == 0:
        axs[i].yaxis.set_ticks(axes_ticks)
        axs[i].set_ylabel(y_label, fontsize=axes_label_size)
        generator_anchored_text = AnchoredText(generator_label, prop=dict(size=generator_label_size),
                                                   loc=generator_label_location)
        axs[i].add_artist(generator_anchored_text)
            
    # Add lambda label as anchored text
    lambda_anchored_text = AnchoredText(lambda_label % lambda_array[i], prop=dict(size=lambda_label_size),
                                            loc=lambda_label_location)
    axs[i].add_artist(lambda_anchored_text)
    
    i += 1
        
    
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
    
file_name = 'r2_operator_contours_TEST'
# Save figure
chdir('Figures/SRG_operators')
f.savefig(file_name+'.pdf', bbox_inches='tight')
chdir(cwd)


# --------------------------------------------------------------------------- #

# Get deuteron wave function in coordinate-space

# Grids of k (col), and r (row) values   
k_cols, r_rows = np.meshgrid(k_array, r_array)
hankel_trans_l0 = np.sqrt(2/np.pi) * k_cols**2 * k_weights * r_rows * \
                  spherical_jn(0, k_cols*r_rows)
hankel_trans_l2 = np.sqrt(2/np.pi) * k_cols**2 * k_weights * r_rows * \
                  spherical_jn(2, k_cols*r_rows)

psi_deuteron_kspace = ob.wave_function(H_initial, -2.22)
u_deuteron_kspace = psi_deuteron_kspace[:ntot] / (k_array*np.sqrt(k_weights))
w_deuteron_kspace = psi_deuteron_kspace[ntot:] / (k_array*np.sqrt(k_weights))

u_deuteron_rspace = hankel_trans_l0 @ u_deuteron_kspace
w_deuteron_rspace = hankel_trans_l2 @ w_deuteron_kspace
psi_squared_rspace = u_deuteron_rspace**2 + w_deuteron_rspace**2

normalization = np.sum(psi_squared_rspace) * dr
print('r-space normalization = %.5f' % normalization)

# RMS radius of deuteron
rms_radius = 0.5 * np.sqrt( np.sum( r_array**2 * psi_squared_rspace * dr ) )
print('rms radius = %.5f fm' % rms_radius)