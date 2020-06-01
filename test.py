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
#   Testing regulated operators.


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from Figures import figures_functions as ff
from Figures import register_colormap
import observables as ob
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Operator functions --- #

def find_q_index(q, k_array):
    
    k_difference_array = np.fabs(k_array - q)
    q_index = k_difference_array.argmin()
    
    return q_index


def momentum_projection_operator(q, k_array, k_weights, channel,
                                 U=np.empty(0)):
    
    # Length of k_array
    m = len(k_array)
        
    # Find index of q in k_array
    q_index = find_q_index(q, k_array)
        
    # Build momentum projection operator 
    operator = np.zeros( (m, m) )
    operator[q_index, q_index] = 1/4
    operator[q_index+1, q_index] = 1/8
    operator[q_index, q_index+1] = 1/8
    operator[q_index-1, q_index] = 1/8
    operator[q_index, q_index-1] = 1/8
    operator[q_index+1, q_index+1] = 1/16
    operator[q_index-1, q_index+1] = 1/16
    operator[q_index+1, q_index-1] = 1/16
    operator[q_index-1, q_index-1] = 1/16
    
    # Divide by momenta/weights
    factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
    row, col = np.meshgrid(factor_array, factor_array)
    operator = operator / row / col
    
    # Build coupled channel operator 
    if lsp.coupled_channel(channel):
    
        # Matrix of zeros (m x m) for coupled-channel operator
        o = np.zeros( (m, m) )
    
        # Build coupled channel operator
        operator = np.vstack( ( np.hstack( (operator, o) ),
                                np.hstack( (o, operator) ) ) )
            
    # Evolve operator by applying unitary transformation U
    if U.any():
        operator = U @ operator @ U.T

    return operator


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


def r2_operator(k_array, k_weights, r_array, dr, U=np.empty(0)):
    
    # Cutoff parameter a
    #a = 4.0
    a = 10.0
    regulator = np.exp( -r_array**2 / a**2 )
    #regulator = 1
        
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


# --- Plotting functions --- #

def momentum_projection_contours(q, kvnn, channel, generators, lambda_array,
                                 contour_type='contourf'):

    # --- Set-up --- #
    
    # Load momentum, weights, and initial Hamiltonian
    k_array, k_weights = lsp.load_momentum(kvnn, channel)
    # Length of k_array
    ntot = len(k_array)
    H_initial = lsp.load_hamiltonian(kvnn, channel)
    # Divide out these factors to present mesh-independent result
    factor_array = k_array * np.sqrt(k_weights) * np.sqrt(2/np.pi)
    row, col = np.meshgrid(factor_array, factor_array)
    
    # Size of figure
    row_number = len(generators)
    col_number = len(lambda_array)
    figure_size = (4*col_number, 3.5*row_number) # Extra width for colorbar

    # Axes limits
    axes_max = 4.0
    axes_lim = [0.0, axes_max]
        
    # Axes ticks, labels, and fontsizes
    x_label = "k' [fm" + r'$^{-1}$' + ']'
    y_label = 'k [fm' + r'$^{-1}$' + ']'
    axes_label_size = 18
    axes_stepsize = 1.0 # Step-size in labeling tick marks
    axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
    axes_tick_size = 18
    
    # Colorbar ticks, label, and fontsize
    if q < axes_max / 2:
        mx = 0.1
        mn = -0.1
    else:
        mx = 0.01
        mn = -0.01
    levels_number = 61
    levels = np.linspace(mn, mx, levels_number)
    levels_ticks = np.linspace(mn, mx, 9)
    if q < axes_max / 2:
        levels_ticks_strings = ['%.3f' % tick for tick in levels_ticks]
    else:
        levels_ticks_strings = ['%.4f' % tick for tick in levels_ticks]
    colorbar_label = '[fm' + r'$^6$' + ']'
    colorbar_label_size = 22
    colorbar_tick_size = 20
    
    # Color scheme for contour plots
    #color_style = 'jet'
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
            # SRG calls function which builds U(s) out of un-evolved and evolved eigenvectors
            if generator == 'Block-diag':
                H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=1.0,
                                                 lambda_bd=lamb)
            else:
                H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=lamb)
                
            U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
            
            # Evolved momentum projection operator
            operator = momentum_projection_operator(q, k_array, k_weights, channel, U_matrix)
            # Take only the upper sub-block if coupled-channel 
            if lsp.coupled_channel(channel):
                operator = operator[:ntot, :ntot]
            # Divide by k_i * k_j * sqrt( w_i * w_j ) for mesh-independent result
            operator = operator / row / col
                
            # Interpolate the operator through 0 to axes_max for smoother looking figure (the extension _int means 
            # interpolated)
            k_array_int, operator_int = ff.interpolate_matrix(k_array, operator, axes_max)
            
            # Store in dictionary with generator and lamb as keys
            d[generator][lamb] = operator_int

        
    # --- Plot data --- #
    
    # Initialize figure
    plt.close('all')
    f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True, figsize=figure_size)
    
    # Loop over generators and lambda's keeping track of indices
    for i, generator in enumerate(generators):
        for j, lamb in enumerate(lambda_array):
            
            # Use contourf
            if contour_type == 'contourf':
                c = axs[i, j].contourf(k_array_int, k_array_int, d[generator][lamb], levels, cmap=color_style, 
                                       extend='both')
            # Otherwise use pcolormesh
            else:
                c = axs[i, j].pcolormesh(k_array_int, k_array_int, d[generator][lamb], cmap=color_style, vmin=mn, 
                                         vmax=mx, rasterized=True)
                                         
            # Specify axes limits
            axs[i, j].set_xlim( axes_lim )
            axs[i, j].set_ylim( axes_lim )
                     
            # On the top row, set and label x-axis
            if i == 0:
                                         
                # Specify axes tick marks
                axs[i, j].xaxis.set_ticks(axes_ticks)
                axs[i, j].xaxis.set_ticklabels(axes_ticks)
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
                axs[i, j].yaxis.set_ticklabels(axes_ticks)
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


def momentum_projection_slices(q, channel, kvnns, generators, lambda_array):

    # --- Set-up --- #

    # Size of figure
    row_number = 2 # For diagonal (top) and far off-diagonal (bottom)
    col_number = len(lambda_array)
    figure_size = (4*col_number, 4*row_number)
    
    # Limits of x axis
    xlim = [0.0, 4.0]
    
    # Axes ticks, labels, and fontsizes
    # x-axis
    x_label = 'k [fm' + r'$^{-1}$' + ']'
    x_label_size = 18
    x_stepsize = 1.0 # Step-size in labeling tick marks
    x_ticks = np.arange(0.0, xlim[1] + x_stepsize, x_stepsize)
    # y-axis
    y_diag_label = r'$a^{\dagger}_q a_q$' + '(k,k) [fm' + r'$^6$' + ']'
    y_off_diag_label = r'$a^{\dagger}_q a_q$' + '(k,0) [fm' + r'$^6$' + ']'
    y_label_size = 20
    axes_tick_size = 16
    
    # Curve width
    curve_width = 2.0


    # --- Load operators --- #
    
    # Initialize dictionary to store evolved operators and momentum arrays
    d = {}
    
    # Loop over kvnns
    for kvnn in kvnns:
    
        d[kvnn] = {}
        
        # Load initial Hamiltonian
        H_initial = lsp.load_hamiltonian(kvnn, channel)
        
        # Load momentum and weights
        k_array, k_weights = lsp.load_momentum(kvnn, channel)
        # Length of momentum array
        ntot = len(k_array)
        # Build factor_array to divide out weights/momenta
        factor_array = k_array * np.sqrt( (2 * k_weights) / np.pi )
        row, col = np.meshgrid(factor_array, factor_array)
        
        # Store in dictionary with kvnn as key
        d[kvnn]['k_array'] = k_array
        
        for generator in generators:
            
            d[kvnn][generator] = {}
        
            # Loop over lambda values
            for lamb in lambda_array:
            
                # Split dictionary further at lamb key for diagonal and far off-diagonal elements
                d[kvnn][generator][lamb] = {}
            
                # Load evolved Hamiltonian
                if generator == 'Block-diag':
                    H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=1.0, 
                                                     lambda_bd=lamb)
                else:
                    H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=lamb)
                
                # SRG calls function which builds U(s) out of un-evolved and evolved eigenvectors
                U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
                
                # Evolved momentum projection operator
                operator = momentum_projection_operator(q, k_array, k_weights, channel, U_matrix)
                # Take only the upper sub-block if coupled-channel 
                if lsp.coupled_channel(channel):
                    operator = operator[:ntot, :ntot]
                # Divide by 2/pi * k_i * k_j * sqrt( w_i * w_j ) for mesh-independent result
                operator = operator / row / col
                
                # Save diagonal and far off-diagonal elements to dictionary
                d[kvnn][generator][lamb]['diag'] = np.diag( operator )[:ntot]
                d[kvnn][generator][lamb]['off-diag'] = operator[:ntot, 0]


    # --- Plot data --- #
    
    # Initialize figure
    plt.close('all')
    f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True, figsize=figure_size)
    
    # Loop for diagonal and off-diagonals
    for i in range(2):
        # Loop over lambda's, kvnns, and generators keeping track of lambda indices
        for j, lamb in enumerate(lambda_array):
            for m, kvnn in enumerate(kvnns):
                curve_color = ff.xkcd_colors(m) # Vary the curve color for kvnn
                for n, generator in enumerate(generators):
                
                    curve_style = ff.line_styles(n) # Vary the curve style for generator
                    # Legend labels (kvnn for top row and generator for bottom row) and slice key for dictionary
                    if i == 0:
                        if m == 0: # Don't repeat labels for each generator
                            curve_label = ff.generator_label_conversion(generator)
                        else:
                            curve_label = ''
                        slice_key = 'diag'
                    else:
                        if n == 0: # Don't repeat labels for each kvnn
                            curve_label = ff.kvnn_label_conversion(kvnn)
                        else:
                            curve_label = ''
                        slice_key = 'off-diag'

                    # Plot slice
                    axs[i, j].plot(d[kvnn]['k_array'], d[kvnn][generator][lamb][slice_key], color=curve_color, 
                                   label=curve_label, linestyle=curve_style, linewidth=curve_width)
                                         
            # Specify x-axis limits
            axs[i, j].set_xlim( xlim )
        
            # On the left column, label y-axis
            if j == 0:
                                      
                # Set y-axis label
                if i == 0:         
                    axs[i, j].set_ylabel(y_diag_label, fontsize=y_label_size)
                else:
                    axs[i, j].set_ylabel(y_off_diag_label, fontsize=y_label_size)
                                         
            # On the bottom row,  set and label x-axis
            if i == 1:
                                         
                # Specify axes tick marks
                axs[i, j].xaxis.set_ticks(x_ticks)
                axs[i, j].xaxis.set_ticklabels(x_ticks)
                                         
                # Prevent overlapping x-axis tick marks
                if j < col_number - 1:
                    xticks = axs[i, j].xaxis.get_major_ticks()
                    xticks[-1].set_visible(False)
                    
                # Set x-axis label
                axs[i, j].set_xlabel(x_label, fontsize=x_label_size)
                
            # Prevent overlapping y-axis tick marks
            yticks = axs[0, 0].yaxis.get_major_ticks()
            yticks[0].set_visible(False)     
                    
            # Enlarge axes tick marks
            axs[i, j].tick_params(labelsize=axes_tick_size)
            
    # Amount of white space in-between sub-plots
    f.subplots_adjust(hspace=0.0, wspace=0.0)
                    
    return f, axs


def r2_diff_contours(kvnn, generators, lambda_array, contour_type='contourf'):
   
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
    #initial_operator = r2_operator(k_array, k_weights, r_array, dr)
    
    # Size of figure
    row_number = len(generators)
    col_number = len(lambda_array)
    figure_size = (4*col_number, 3.5*row_number) # Extra width for colorbar

    # Axes limits
    #axes_max = 0.4
    #axes_max = 1.0
    axes_max = 10.0
    axes_lim = [0.0, axes_max]
        
    # Axes ticks, labels, and fontsizes
    x_label = "k' [fm" + r'$^{-1}$' + ']'
    y_label = 'k [fm' + r'$^{-1}$' + ']'
    axes_label_size = 18
    #axes_stepsize = 0.1 # Step-size in labeling tick marks
    #axes_stepsize = 0.2 # Step-size in labeling tick marks
    axes_stepsize = 2.0 # Step-size in labeling tick marks
    axes_ticks = np.arange(0.0, axes_max + axes_stepsize, axes_stepsize)
    axes_ticks_strings = ['%.1f' % tick for tick in axes_ticks]
    axes_tick_size = 18
    
    # Colorbar ticks, label, and fontsize
    mx = 5e2
    mn = -5e2
    # mx = 5e3
    # mn = -5e3
    #mx = 1e5
    #mx = 1e6
    #mx = 1e4
    #mn = -1e5
    #mn = -1e6
    #mn = -1e4
    levels_number = 61
    levels = np.linspace(mn, mx, levels_number)
    levels_ticks = np.linspace(mn, mx, 9)
    levels_ticks_strings = ['%.0f' % tick for tick in levels_ticks]
    colorbar_label = '[fm' + r'$^5$' + ']'
    colorbar_label_size = 22
    colorbar_tick_size = 20
    
    # Color scheme for contour plots
#     color_style = 'jet'
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
            # SRG calls function which builds U(s) out of un-evolved and evolved eigenvectors
            if generator == 'Block-diag':
#                 H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=1.0,
#                                                  lambda_bd=lamb)
                H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=1.5,
                                                 lambda_bd=lamb)
            else:
                H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg', generator=generator, lamb=lamb)
                
            U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
            
            # Evolved momentum projection operator
            evolved_operator = r2_operator(k_array, k_weights, r_array, dr, U_matrix)
        
            # Calculate difference from evolved and initial operators
            #operator_diff = (evolved_operator - initial_operator)
            # Take only the upper sub-block if coupled-channel 
            if lsp.coupled_channel(channel):
                #operator_diff = operator_diff[:ntot, :ntot]
                evolved_operator = evolved_operator[:ntot, :ntot] / row / col
            # Divide by k_i * k_j * sqrt( w_i * w_j ) for mesh-independent result
            #operator_diff = operator_diff / row / col
        
            # Interpolate the operator through 0 to axes_max for smoother looking figure (the extension _int means 
            # interpolated)
            #k_array_int, operator_diff_int = ff.interpolate_matrix(k_array, operator_diff, axes_max)
            k_array_int, evolved_operator_int = ff.interpolate_matrix(k_array, evolved_operator, axes_max)
        
            # Store in dictionary with generator and lamb as keys
            d[generator][lamb] = evolved_operator_int
        
    
    # --- Plot data --- #
    
    # Initialize figure
    plt.close('all')
    f, axs = plt.subplots(row_number, col_number, sharex=True, sharey=True, figsize=figure_size)
    
    # Loop over generators and lambda's keeping track of indices
    for i, generator in enumerate(generators):
        for j, lamb in enumerate(lambda_array):
            
            # Use contourf
            if contour_type == 'contourf':
                c = axs[i, j].contourf(k_array_int, k_array_int, d[generator][lamb], levels, cmap=color_style, 
                                       extend='both')
            # Otherwise use pcolormesh
            else:
                c = axs[i, j].pcolormesh(k_array_int, k_array_int, d[generator][lamb], cmap=color_style, vmin=mn, 
                                         vmax=mx, rasterized=True)
                                         
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
#q = 3.0
q = 0.3


# --- Test regulated operators --- #


# # Contours of evolved momentum projection operator under RKE N4LO (450 MeV) transformations where q = 3 fm^-1
# lambda_array = np.array([6.0, 3.0, 2.0, 1.5])
# f, axs = momentum_projection_contours(q, kvnn_default, channel_default, 
#                                       generators, lambda_array)

# # Add generator label to each subplot on the 1st column
# generator_label_size = 17
# generator_label_location = 'upper left'
# for i, generator in enumerate(generators):
#     generator_label = ff.generator_label_conversion(generator)
#     anchored_text = AnchoredText(generator_label, loc=generator_label_location,
#                                   prop=dict(size=generator_label_size))
#     axs[i, 0].add_artist(anchored_text)

# # Add \lambda label to each sub-plot
# lambda_label_size = 17
# lambda_label_location = 'lower left'
# for i, generator in enumerate(generators):
#     for j, lamb in enumerate(lambda_array):
#         if generator == 'Block-diag':
#             # Labels the block-diagonal cutoff \Lambda_BD
#             lambda_label = ff.lambda_label_conversion(lamb, block_diag_bool=True)
#         else:
#             # Labels the evolution parameter \lambda
#             lambda_label = ff.lambda_label_conversion(lamb)
#         anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, prop=dict(size=lambda_label_size))
#         axs[i, j].add_artist(anchored_text)
        
# plt.show()
        
        
# # Diagonal and far off-diagonal slices of momentum projection operator under EMN N4LO (500 MeV), RKE N4LO 
# # (450 MeV), Gezerlis N2LO (1 fm) transformations with q = 3.0 fm^-1
# lambda_array = np.array([6.0, 3.0, 2.0, 1.5])
# f, axs = momentum_projection_slices(q, channel_default, kvnns_default, 
#                                     generators, lambda_array)

# # Set the y-axis limit and tickmarks (this will vary based on q value)
# #ylim = [-0.003, 0.012]
# ylim = [-0.03, 0.12]
# #y_stepsize = 0.003
# y_stepsize = 0.03
# y_ticks = np.arange(ylim[0], ylim[1] + y_stepsize, y_stepsize)
# y_ticks_labels = ['%.3f' % tick for tick in y_ticks]
# for i in range(2):
#     for j in range(len(lambda_array)):
#         axs[i, j].set_ylim(ylim)
#         if j == 0:
#             axs[i, j].yaxis.set_ticks(y_ticks)
#             axs[i, j].yaxis.set_ticklabels(y_ticks_labels)

# # Add legend for generators to upper left sub-plot
# legend_size = 18
# #legend_location = 'upper left'
# legend_location = 'upper right'
# axs[0, 0].legend(loc=legend_location, frameon=False, fontsize=legend_size)

# # Add legend for kvnns to lower left sub-plot
# legend_size = 18
# #legend_location = 'upper left'
# legend_location = 'upper right'
# axs[1, 0].legend(loc=legend_location, frameon=False, fontsize=legend_size)

# # Add \lambda and \Lambda_BD labels to each sub-plot
# lambda_label = r'$\lambda$' + ', ' + r'$\Lambda_{BD}=%.1f$' + ' fm' + r'$^{-1}$'
# lambda_label_size = 16
# #lambda_label_location = 'lower left'
# lambda_label_location = 'lower right'
# for i in range(2):
#     for j, lamb in enumerate(lambda_array):
#         anchored_text = AnchoredText(lambda_label % lamb, loc=lambda_label_location, 
#                                       prop=dict(size=lambda_label_size), frameon=False)
#         axs[i, j].add_artist(anchored_text)
    
    
# Contours of r^2 operator under RKE N4LO (450 MeV) transformations
lambda_array = np.array([3.0, 2.0, 1.5])
f, axs = r2_diff_contours(kvnn_default, generators, lambda_array)

# Add generator label to each subplot on the 1st column
generator_label_size = 17
#generator_label_location = 'center right'
generator_label_location = 'upper right'
for i, generator in enumerate(generators):
    generator_label = ff.generator_label_conversion(generator)
    anchored_text = AnchoredText(generator_label, loc=generator_label_location,
                                  prop=dict(size=generator_label_size))
    axs[i, 0].add_artist(anchored_text)

# Add \lambda label to each sub-plot
lambda_label_size = 17
#lambda_label_location = 'lower right'
lambda_label_location = 'lower left'
for i, generator in enumerate(generators):
    for j, lamb in enumerate(lambda_array):
        if generator == 'Block-diag':
            # Labels the block-diagonal cutoff \Lambda_BD
            lambda_label = ff.lambda_label_conversion(lamb, block_diag_bool=True)
        else:
            # Labels the evolution parameter \lambda
            lambda_label = ff.lambda_label_conversion(lamb)
        anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, prop=dict(size=lambda_label_size))
        axs[i, j].add_artist(anchored_text)
        
# Calculate RMS radius of deuteron
H_matrix = lsp.load_hamiltonian(kvnn_default, channel_default)
k_array, k_weights = lsp.load_momentum(kvnn_default, channel_default)
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)
psi = ob.wave_function(H_matrix)
r2_op = r2_operator(k_array, k_weights, r_array, dr)
deuteron_radius = ob.rms_radius_from_rspace(psi, r2_op)
print('r = %.5f fm' % deuteron_radius) # Should give 1.96574 fm


