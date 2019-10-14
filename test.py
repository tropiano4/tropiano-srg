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
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#
#------------------------------------------------------------------------------


from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp


# --- Variables --- #

potentials = ((10, 30.0, 4.0, 120), (111, 8.0, 2.0, 120), (222, 10.0, 2.0, 120))
lambda_array = np.array( [6.0, 3.0, 1.5, 1.0] )
channel = '3S1'
#generator = 'Wegner'
generator = 'Block-diag'

# --- Set-up --- #
    
# Limits of x and y axes (dependent on channel and potential)
xlim = [0.0, 3.0]
ylim = [-4.5, 2.0]
        
# Labels
x_label = 'k [fm' + r'$^{-1}$' + ']'
y_label = 'V(k,k) [fm]'
channel_label = ff.channel_label_conversion(channel)
generator_label = ff.generator_label_conversion(generator)
# Label for lambda and Lambda_BD
if generator == 'Wegner':
    lambda_label = r'$\lambda=%.1f$' + ' fm' + r'$^{-1}$'
else:
    lambda_label = r'$\Lambda=%.1f$' + ' fm' + r'$^{-1}$'
    
# Fontsize for labels and tick marks
x_label_size = 18
y_label_size = 20
legend_label_size = 15
channel_label_size = 22
generator_label_size = 20
lambda_label_size = 17
    
# Location of labels
legend_label_location = 'lower right'
generator_label_location = 'lower right'
channel_label_location = 'lower right'
lambda_label_location = 'upper left'
    
# Size of figure
row_number = 1
col_number = len(lambda_array)
figure_size = (4*col_number, 4*row_number) # (width, height)
    
# Initialize file name of figure
file_name = 'potential_%s_%s_kvnns' % ('diag', '3S1') # This reads 'potential_line_channel_kvnns'
kvnn_list = [10, 111, 222]

# Current working directory
cwd = getcwd()


# --- Load data and plot lines --- #
    
# Initialize figure
plt.close('all')
f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True,
                      figsize=figure_size)
    
# Loop over potential specifications
for potential in potentials:
        
    # Set kvnn, kmax, kmid, ntot, and generator
    kvnn = potential[0]
    kmax = potential[1]
    kmid = potential[2]
    ntot = potential[3]
            
    # Load momentum
    k_array, _ = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
        
    # Curve labels and styles
            
    # Curve color depends on potential
    if kvnn in [10, 900]:
        curve_color = 'xkcd:black'
    elif kvnn in [105, 106, 107, 110, 111, 112, 901]:
        curve_color = 'xkcd:red'
    elif kvnn in [222, 224, 902]:
        curve_color = 'xkcd:blue'
    curve_style = 'solid'
    potential_label = ff.kvnn_label_conversion(kvnn, full_label=False)
        
    # Loop over lambda
    i = 0 # Sub-plot number
    for lamb in lambda_array:

        # Load evolved potential
        if generator == 'Block-diag':
            V_matrix = lp.load_potential(kvnn, channel, kmax, kmid, ntot, 
                                         'srg', generator, lambda_array[-1],
                                         lambda_bd=lamb)
        else:
            V_matrix = lp.load_potential(kvnn, channel, kmax, kmid, ntot,
                                         'srg', generator, lamb)
            
        # Take a slice of the potential to plot
        # (indexing :ntot keeps the same dimension of k_array in the case of a coupled-channel potential)
        V_vector = np.diag( V_matrix[:ntot, :ntot] )
 
        # Add sub-plot to figure
        if i == 2: # 3rd sub-plot - label potential
            axs[i].plot(k_array, V_vector, color=curve_color, linestyle=curve_style, label=potential_label)
        else: # Middle sub-plots - no labels
            axs[i].plot(k_array, V_vector, color=curve_color, linestyle=curve_style)
                
        i += 1
        
        
# --- Set figure specifications and save --- #
    
# Loop over sub-plots
for j in range( len(lambda_array) ):
        
    # Specify axes limits
    axs[j].set_xlim(xlim)
    axs[j].set_ylim(ylim)
        
    # Prevent overlapping x-axis tick marks unless it's the last sub-plot
    if j != ( len(lambda_array) - 1 ):
        xticks = axs[j].xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
            
    # Set axes labels
    axs[j].set_xlabel(x_label, fontsize=x_label_size)
        
    # Only specify y label and potential labels as legend for 1st sub-plot
    if j == 0:
        # Add y label
        axs[j].set_ylabel(y_label, fontsize=y_label_size)
        
    if j == 1:
        # Add channel label as anchored text to 2nd sub-plot for diag plot only
        channel_anchored_text = AnchoredText(channel_label, 
                                         prop=dict(size=channel_label_size),
                                         loc=channel_label_location,
                                         frameon=False)
        axs[j].add_artist(channel_anchored_text)

    if j == 2:
        # Add legend for potentials
        axs[j].legend(loc=legend_label_location, frameon=False, fontsize=legend_label_size)
            
    if j == ( len(lambda_array) - 1 ):
        # Add channel label as anchored text to 2nd sub-plot for diag plot only
        generator_anchored_text = AnchoredText(generator_label, 
                                               prop=dict(size=generator_label_size),
                                               loc=generator_label_location,
                                               frameon=False)
        axs[j].add_artist(generator_anchored_text)
            
    # Add lambda's label as anchored text
    lambda_anchored_text = AnchoredText(lambda_label % lambda_array[j],
                                        prop=dict(size=lambda_label_size),
                                        loc=lambda_label_location, frameon=False)
    # Add lambda
    axs[j].add_artist(lambda_anchored_text)
    
# Amount of white space in-between sub-plots
f.subplots_adjust(hspace=0.0, wspace=0.0)
    
# Name of the file
for kvnn in kvnn_list:
    file_name += '_%d' % kvnn
# Add last value of lambda to file name
file_name += '_%s_lamb%.1f' % (generator, lambda_array[-1])
# Replace '.' with ',' in file name since LaTeX doesn't like periods
file_name = ff.replace_periods_with_commas(file_name)
    
# Save figure
chdir('Figures/SRG_potentials')
f.savefig(file_name+'.pdf', bbox_inches='tight')
chdir(cwd)