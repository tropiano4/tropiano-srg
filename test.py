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
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Testing relative strength of the r^2 operator.


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from Figures import figures_functions as ff
from Figures import register_colormap
import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Operator functions --- #


def hankel_transformation(channel, k_array, r_array, dr):
        
    # L = 0 (0th spherical Bessel function)
    if channel[1] == 'S':
        L = 0
    # L = 1
    elif channel[1] == 'P':
        L = 1
    # L = 2
    elif channel[1] == 'D':
        L = 2
        
    # r_array column vectors and k_array row vectors where both grids are
    # n x m matrices
    r_cols, k_rows = np.meshgrid(r_array, k_array)
        
    M = np.sqrt(dr) * r_cols * spherical_jn(L, k_rows * r_cols)

    return M


def r2_operator(k_array, k_weights, r_array, dr, U=np.empty(0), reg=False):
    
    # Cutoff parameter a
    a = 6.0
    if reg:
        regulator = np.exp( -r_array**2 / a**2 )
    else:
        regulator = 1
        
    # Initialize r^2 in coordinate-space first where r^2 is a diagonal matrix
    r2_coordinate_space = np.diag(r_array**2) * regulator
     
    # Transform operator to momentum-space
    s_wave_trans = hankel_transformation('3S1', k_array, r_array, dr)
    d_wave_trans = hankel_transformation('3D1', k_array, r_array, dr)
    
    # Each variable here corresponds to a sub-block of the coupled channel 
    # matrix
    ss_block = s_wave_trans @ r2_coordinate_space @ s_wave_trans.T
    dd_block = d_wave_trans @ r2_coordinate_space @ d_wave_trans.T
    
    # Grids of momenta and weights
    factor_array = np.concatenate( (np.sqrt(k_weights) * k_array, 
                                    np.sqrt(k_weights) * k_array) ) * \
                   np.sqrt(2/np.pi)
    row, col = np.meshgrid(factor_array, factor_array)
        
    # Length of k_array
    n = len(k_array)
        
    # Matrix of zeros (m x m) for coupled-channel operator
    o = np.zeros( (n, n) )
        
    # Build coupled channel operator with momenta/weights
    r2_momentum_space = np.vstack( ( np.hstack( (ss_block, o) ),
                                     np.hstack( (o, dd_block) ) ) ) * row * col
    
    # Evolve operator by applying unitary transformation U
    if U.any():
        r2_momentum_space = U @ r2_momentum_space @ U.T
        
    return r2_momentum_space


def rms_radius_from_rspace(psi, r2_operator):
    
    r2 = psi.T @ r2_operator @ psi
    
    return 0.5 * np.sqrt(r2)


# --- Plotting functions --- #


def r2_contours(kvnn, generators, lambda_array, reg=False):
   
    # --- Set-up --- #
    
    channel = '3S1'
    
    # Load momentum, weights, and initial Hamiltonian
    k_array, k_weights = lsp.load_momentum(kvnn, channel)
    # Length of k_array
    ntot = len(k_array)
    H_initial = lsp.load_hamiltonian(kvnn, channel)
    # Divide out these factors to present mesh-independent result
    #factor_array = k_array * np.sqrt(k_weights) * np.sqrt(2/np.pi)
    factor_array = np.sqrt(k_weights) * np.sqrt(2/np.pi)
    row, col = np.meshgrid(factor_array, factor_array)
    
    # Specify r_array
    r_min = 0.005
    r_max = 30.2
    dr = 0.005
    r_array = np.arange(r_min, r_max + dr, dr)
    
    # Calculate initial operator
    #initial_operator = r2_operator(k_array, k_weights, r_array, dr, reg)
    
    # Size of figure
    row_number = len(generators)
    col_number = len(lambda_array)
    figure_size = (4*col_number, 3.5*row_number) # Extra width for colorbar

    # Axes limits
    #axes_max = 0.4
    #axes_max = 1.0
    axes_max = 10.0
    axes_lim = [0.0, axes_max]
    #axes_stepsize = 0.1 # Step-size in labeling tick marks
    #axes_stepsize = 0.2 # Step-size in labeling tick marks
    axes_stepsize = 2.0 # Step-size in labeling tick marks
        
    # Axes ticks, labels, and fontsizes
    x_label = "k' [fm" + r'$^{-1}$' + ']'
    y_label = 'k [fm' + r'$^{-1}$' + ']'
    axes_label_size = 18
    axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
    axes_ticks_strings = ['%.1f' % tick for tick in axes_ticks]
    axes_tick_size = 18
    
    # Colorbar ticks, label, and fontsize
    mx = 5e2
    mn = -5e2
    # mx = 5e3
    # mn = -5e3
    #mx = 1e4
    #mn = -1e4
    #mx = 1e5
    #mn = -1e5
    #mx = 1e6
    #mn = -1e6
    levels_number = 61
    levels = np.linspace(mn, mx, levels_number)
    levels_ticks = np.linspace(mn, mx, 9)
    levels_ticks_strings = ['%.0f' % tick for tick in levels_ticks]
    colorbar_label = '[fm' + r'$^5$' + ']'
    colorbar_label_size = 22
    colorbar_tick_size = 20
    
    # Color scheme for contour plots
    color_style = 'turbo'
    
    
    # --- Load operators --- #
    
    # Initialize dictionary to store evolved potentials
    d = {}
    
    # Loop over generators
    for generator in generators:
        
        # Store momentum and operator in here
        d[generator] = {}
        
        # Loop over lambda values
        for lamb in lambda_array:
            
            # Load unitary transformation
            # SRG calls function which builds U(s) out of un-evolved and
            # evolved eigenvectors
            if generator == 'Block-diag':
                H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg',
                            generator=generator, lamb=1.0, lambda_bd=lamb)
            else:
                H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg',
                            generator=generator, lamb=lamb)
                
            U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
            
            # Evolved momentum projection operator
            evolved_operator = r2_operator(k_array, k_weights, r_array, dr,
                                           U_matrix, reg)
        
            # Calculate difference from evolved and initial operators
            #operator_diff = (evolved_operator - initial_operator)
            # Take only the upper sub-block if coupled-channel 
            evolved_operator = evolved_operator[:ntot, :ntot] / row / col
        
            # Interpolate the operator through 0 to axes_max for smoother 
            # looking figure (the extension _int means interpolated)
            #k_array_int, operator_diff_int = ff.interpolate_matrix(k_array,
            #                                 operator_diff, axes_max)
            k_array_int, evolved_operator_int = ff.interpolate_matrix(k_array, 
                                                evolved_operator, axes_max)
        
            # Store in dictionary with generator and lamb as keys
            d[generator][lamb] = evolved_operator_int
        
    
    # --- Plot data --- #
    
    # Initialize figure
    plt.close('all')
    f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True,
                          figsize=figure_size)
    
    # Loop over generators and lambda's keeping track of indices
    for i, generator in enumerate(generators):
        for j, lamb in enumerate(lambda_array):
            
            # Use contourf
            c = axs[i, j].contourf(k_array_int, k_array_int,
                                   d[generator][lamb], levels,
                                   cmap=color_style, extend='both')

                                         
            # Specify axes limits
            axs[i, j].set_xlim( axes_lim )
            axs[i, j].set_ylim( axes_lim )
                     
            # On the top row, set and label x-axis
            if i == 0:
                                         
                # Specify axes tick marks
                axs[i, j].xaxis.set_ticks(axes_ticks)
                axs[i, j].xaxis.set_ticklabels(axes_ticks_strings)
                # Switch from bottom to top
                axs[i, j].xaxis.set_label_position('top')
                axs[i, j].xaxis.tick_top()
                axs[i, j].tick_params(labeltop=True, labelsize=axes_tick_size)
                                         
                # Prevent overlapping x-axis tick marks
                if j < col_number - 1:
                    xticks = axs[i, j].xaxis.get_major_ticks()
                    xticks[-1].set_visible(False)

                # Set x-axis label
                axs[i, j].set_xlabel(x_label, fontsize=axes_label_size)
                                         
            # On the left column, set and label y-axis
            if j == 0:
                                         
                # Specify axes tick marks
                axs[i, j].yaxis.set_ticks(axes_ticks)
                axs[i, j].yaxis.set_ticklabels(axes_ticks_strings)
                axs[i, j].tick_params(labelsize=axes_tick_size)
                                      
                # Prevent overlapping y-axis tick marks
                if i < row_number - 1:
                    yticks = axs[i, j].yaxis.get_major_ticks()
                    yticks[-1].set_visible(False)
                                         
                # Set y-axis label
                axs[i, j].set_ylabel(y_label, fontsize=axes_label_size)
                                         
            # On the bottom row, switch x-axis from bottom to top
            if i == row_number - 1:
                                         
                axs[i, j].xaxis.tick_top()
                axs[i, j].tick_params(labeltop=False, labelsize=axes_tick_size)

    # Invert y-axis
    plt.gca().invert_yaxis()
                                         
    # Amount of white space in-between sub-plots
    f.subplots_adjust(hspace=0.0, wspace=0.0)
                                         
    # Set colorbar axe
    f.subplots_adjust(right=0.8) # Adjust for colorbar space
    cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
                                         
    # Set colorbar
    cbar = f.colorbar(c, cax=cbar_ax, ticks=levels_ticks)
    cbar.ax.tick_params(labelsize=colorbar_tick_size)
    cbar.ax.set_yticklabels(levels_ticks_strings)
                                         
    # Set colorbar label
    cbar.ax.set_title(colorbar_label, fontsize=colorbar_label_size)

    return f, axs


# --- Fixed variables --- #

kvnns_default = [79, 111, 222]
kvnn_default = 111
channel_default = '3S1'
generators = ['Wegner', 'Block-diag']


# --- Test regulated operators --- #


# # Contours of r^2 operator under RKE N4LO (450 MeV) transformations
# lambda_array = np.array([3.0, 2.0, 1.5])
# f, axs = r2_contours(kvnn_default, generators, lambda_array)

# # Add generator label to each subplot on the 1st column
# generator_label_size = 17
# #generator_label_location = 'center right'
# generator_label_location = 'upper right'
# for i, generator in enumerate(generators):
#     generator_label = ff.generator_label_conversion(generator)
#     anchored_text = AnchoredText(generator_label, loc=generator_label_location,
#                                   prop=dict(size=generator_label_size))
#     axs[i, 0].add_artist(anchored_text)

# # Add \lambda label to each sub-plot
# lambda_label_size = 17
# #lambda_label_location = 'lower right'
# lambda_label_location = 'lower left'
# for i, generator in enumerate(generators):
#     for j, lamb in enumerate(lambda_array):
#         if generator == 'Block-diag':
#             # Labels the block-diagonal cutoff \Lambda_BD
#             lambda_label = ff.lambda_label_conversion(lamb, 
#                                                       block_diag_bool=True)
#         else:
#             # Labels the evolution parameter \lambda
#             lambda_label = ff.lambda_label_conversion(lamb)
#         anchored_text = AnchoredText(lambda_label, loc=lambda_label_location,
#                                      prop=dict(size=lambda_label_size))
#         axs[i, j].add_artist(anchored_text)
        
# Calculate RMS radius of deuteron
#kvnn = 111
kvnn = 6
#evolve = False
evolve = True
reg = False
#reg = True
H_initial = lsp.load_hamiltonian(kvnn, channel_default)
H_evolved = lsp.load_hamiltonian(kvnn, channel_default, method='srg',
                                 generator='Wegner', lamb=2.0)
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
k_array, k_weights = lsp.load_momentum(kvnn, channel_default)
ntot = len(k_array)
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)
if evolve:
    psi = ob.wave_function(H_initial, U=U_matrix)
    r2_op = r2_operator(k_array, k_weights, r_array, dr, U=U_matrix, reg=reg)
else:
    psi = ob.wave_function(H_initial)
    r2_op = r2_operator(k_array, k_weights, r_array, dr, reg=reg)
deuteron_radius_exact = ob.rms_radius_from_rspace(psi, r2_op)
print('r = %.5f fm' % deuteron_radius_exact) # Should give 1.96574 fm


# Test strength of regions of r^2 by isolating integral to different
# regions in k and k'

k_points = np.zeros(11)
k_points[0] = k_array[0]
for i in range(1, 10):
    k_points[i] = i
k_points[10] = k_array[-1]
m = len(k_points)
rel_errors = np.zeros((m, m))

for i, k_init in enumerate(k_points):
    for j, k_final in enumerate(k_points):

        k_init_index = op.find_q_index(k_init, k_array)
        k_init_index_l2 = k_init_index+ntot
        k_final_index = op.find_q_index(k_final, k_array)
        k_final_index_l2 = k_final_index+ntot

        psi_resized = np.concatenate( (psi[k_init_index:k_final_index], 
                                       psi[k_init_index_l2:k_final_index_l2]) )
        r2_ss = r2_op[k_init_index:k_final_index, k_init_index:k_final_index]
        r2_sd = r2_op[k_init_index:k_final_index,
                      k_init_index_l2:k_final_index_l2]
        r2_ds = r2_op[k_init_index_l2:k_final_index_l2,
                      k_init_index:k_final_index]
        r2_dd = r2_op[k_init_index_l2:k_final_index_l2,
                      k_init_index_l2:k_final_index_l2]
        r2_resized = np.vstack( ( np.hstack( (r2_ss, r2_sd) ),
                                  np.hstack( (r2_ds, r2_dd) ) ) )
        deuteron_radius = ob.rms_radius_from_rspace(psi_resized, r2_resized)
        rel_errors[i, j] = abs( (deuteron_radius-deuteron_radius_exact) / \
                                 deuteron_radius_exact )
            
        # If k_init > k_final then the integral doesn't make any sense
        # Set rel_error to 0 here
        if k_init >= k_final:
            rel_errors[i, j] =  100.0
        else:
            print('k_init=%.1f, k_final=%.1f, rel_err=%.3f'%(k_init,k_final,
                                                             rel_errors[i,j]))
            
# Plot the relative errors as a contour

plt.close('all')
f, ax = plt.subplots(figsize=(4, 4))

mx = 1.0
mn = 0.0

c = ax.pcolor(k_points, k_points, rel_errors, vmin=mn, vmax=mx, cmap='Blues')

size = 16
ax.set_ylabel(r'$\Lambda_{initial}$' + ' fm' + r'$^{-1}$', fontsize=size)
ax.set_xlabel(r'$\Lambda_{final}$' + ' fm' + r'$^{-1}$', fontsize=size)
                                         
# Set colorbar axe
f.subplots_adjust(right=0.8) # Adjust for colorbar space
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
                                         
# Set colorbar
cbar = f.colorbar(c, cax=cbar_ax)

plt.show()


# Plot the relative errors as a scatter plot holding \Lambda_init = 0.0

plt.close('all')
f, ax = plt.subplots(figsize=(4, 4))

blue_label = r'$\Lambda_{initial}=0$' + ' fm' + r'$^{-1}$' + ', vary ' + \
             r'$\Lambda_{final}$'
red_label = r'$\Lambda_{final}=10$' + ' fm' + r'$^{-1}$' + ', vary ' + \
            r'$\Lambda_{initial}$'           
             
ax.scatter(k_points, rel_errors[0, :], color='xkcd:blue', marker='o',
           label=blue_label)
ax.scatter(k_points, rel_errors[:, m-1], color='xkcd:red', marker='s',
           label=red_label)
ax.axhline(y=0.0, color='xkcd:black', linestyle='dotted')
ax.legend()
ax.set_xlim([-0.1, 10.1])
ax.set_ylim([-0.1, 1.1])
ax.set_ylabel('Relative error')
ax.set_xlabel(r'$\Lambda$' + ' fm' + r'$^{-1}$')

plt.show()