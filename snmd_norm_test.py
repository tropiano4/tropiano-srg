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
#
#------------------------------------------------------------------------------


# Description of this test:
#   Testing normalization and contributions of \delta U, etc. to single-
#   nucleon momentun distributions.


import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from densities import load_density
from figures import figures_functions as ff
from potentials.vsrg_macos import vnn
from snmd import single_nucleon_momentum_distributions


# Load data from the following directory
data_directory = 'Figures/SRC_physics/Data'

# Figure size
row_number = 1
col_number = 1
figure_size = (4*col_number, 4*row_number)

# Axes labels and fontsize
title_size = 16
x_label = 'q [fm' + r'$^{-1}$' + ']'
x_label_size = 16
y_label = 'proton ' + r'$n^{\lambda}_A(q)/Z$' + ' [fm' + r'$^3$' + ']'
y_label_size = 16
legend_size = 14
legend_location = 'upper right'
    
# Curve width
curve_width = 2.0

# Axes limits
xlim = (0.0, 5.0)
ylim = (1e-2, 1e3)

# Set potential and other inputs
kvnn = 6
lambda_array = np.array([1.35])
# lambda_array = np.array([6.0, 3.0, 2.0, 1.35])
# kmax, kmid, ntot = 10.0, 2.0, 120
kmax, kmid, ntot = 15.0, 3.0, 120
# kmax, kmid, ntot = 30.0, 4.0, 120
# nuclei_list = [ ['O16', 8, 8] ]
nuclei_list = [ ['Pb208', 82, 126] ]
# nuclei_list = ( ('O16', 8, 8), ('Ca40', 20, 20), ('Pb208', 82, 126) )
# nuclei_list = ( ('C12', 6, 6), ('Ca40', 20, 20) )
# nuclei_list = ( ('C12', 6, 6), ('O16', 8, 8), ('Ca40', 20, 20),
#                 ('Ca48', 20, 28), ('Fe56', 26, 30), ('Pb208', 82, 126) )
channels = ('1S0', '3S1')
# channels = ('1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3')

# Load momentum (channel argument doesn't matter here)
q_array, q_weights = vnn.load_momentum(kvnn, '3S1', kmax, kmid, ntot)

interp = True
# interp = False

# Loop over lambda
for lamb in lambda_array:
    
    # Initialize single-nucleon momentum distributions class for given
    # potential
    snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax,
                                                 kmid, ntot, 'Wegner', interp)
    
    print('_'*50)
    print( 'lambda = %s fm^-1\n' % ff.convert_number_to_string(lamb) )

    # Loop over nuclei
    for nucleus in nuclei_list:
    
        nucleus_name = nucleus[0]
        Z = nucleus[1]
        N = nucleus[2]
    
        if interp:
            
            n_array_funcs = snmd.n_lambda_interp(nucleus_name, 'proton', Z, N,
                                                 edf='SLY4')
            
            n_total_array = n_array_funcs[0](q_array)
            n_1_array = n_array_funcs[1](q_array)
            n_delU_array = n_array_funcs[2](q_array)
            n_delU2_array = n_array_funcs[3](q_array)
            
        else:
            
            # Load r values and nucleonic densities (the r_array's are the same)
            R_array, rho_p_array = load_density(nucleus_name, 'proton', Z, N)
            R_array, rho_n_array = load_density(nucleus_name, 'neutron', Z, N)
            dR = R_array[2]-R_array[1]
            
            n_1_array, n_delU_array, n_delU2_array = snmd.n_contributions(
                               q_array, R_array, dR, rho_p_array, rho_n_array )
            # # Calculate nuclear-averaged momentum distributions
            # n_array_cont = snmd.n_lambda(q_array, r_array, rho_p_array,
            #                          rho_n_array)
    
            # # Proton contributions
            # n_total_array = n_array_cont[0] 
            # n_1_array = n_array_cont[:, 1]
            # n_delU_array = n_array_cont[:, 2]
            # n_delU2_array = n_array_cont[:, 3]
            n_total_array = n_1_array + n_delU_array + n_delU2_array
    
        # Compute normalizations here
        # Proton
        factor = 4*np.pi/(2*np.pi)**3
        p_total_norm = factor * np.sum(q_array**2 * q_weights * n_total_array)
        p_1_norm = factor * np.sum(q_array**2 * q_weights * n_1_array)
        p_delU_norm = factor * np.sum(q_array**2 * q_weights * n_delU_array)
        p_delU2_norm = factor * np.sum(q_array**2 * q_weights * n_delU2_array)
        
        p_1_delU_norm = factor * np.sum( q_array**2 * q_weights * \
                                         (n_1_array + n_delU_array) )
        r_norm = p_1_delU_norm / p_1_norm
        p_1_delU_norm2 = factor * np.sum( q_array**2 * q_weights * \
                                         (n_1_array + n_delU_array/2) )
        r_norm2 = (p_1_delU_norm2 / p_1_norm)**2
    
        # Print normalizations
        print('-'*50)
        print(nucleus_name)
        print('Total proton normalization = %.5f' % p_total_norm)
        print('1 term proton normalization = %.5f' % p_1_norm)
        print('\delta U term proton normalization = %.5f' % p_delU_norm)
        print('\delta U^2 term proton normalization = %.5f' % p_delU2_norm)
        
        print('(1 term + \delta U term) norm / 1 term norm = %.5f' % r_norm)
        print('( (1 term + \delta U term/2) norm / 1 term norm )^2 = %.5f' % r_norm2)
    
        # Plot with respect to AV18 data
        plt.close('all')
        f, ax = plt.subplots(figsize=figure_size)
        
        # Add curve to figure
        ax.set_yscale('log')
        ax.plot(q_array, n_total_array/Z, color='xkcd:grey', label='total',
                linewidth=curve_width)
        ax.plot(q_array, n_1_array/Z, color='xkcd:blue', label='1',
                linestyle='dotted', linewidth=curve_width)
        ax.plot(q_array, abs(n_delU_array)/Z, color='xkcd:green',
                label=r'$|\delta U|$', linestyle='dashed', linewidth=curve_width)
        ax.plot(q_array, n_delU2_array/Z, color='xkcd:red', linestyle='dashdot',
                label=r'$\delta U \delta U^{\dagger}$', linewidth=curve_width)
        
        # Add AV18 data with error bars
        if nucleus_name in ('C12', 'O16', 'Ca40'):
            av18_data = np.loadtxt('data/qmc'+'/'+'AV18_%s_snmd.txt'
                                    % nucleus_name)
            q_array_av18 = av18_data[:, 0] # fm^-1
            n_p_array_av18 = av18_data[:, 1] / Z
            error_bars_array_av18 = av18_data[:, 2] / Z
                
            # AV18 data with error bars
            ax.errorbar(q_array_av18, n_p_array_av18, yerr=error_bars_array_av18,
                        color='xkcd:black', label='AV18', linestyle='', marker='.')

        # Specify axes limits
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
        # Title of plot
        # plot_title = ff.nuclei_label_conversion(nucleus_name) # Labels the nucleus
        plot_title = ff.nuclei_label_conversion(nucleus_name) + ' and ' + \
                     ff.lambda_label_conversion(lamb)
        ax.set_title(plot_title, fontsize=title_size)
        
        # Set axes labels and legend
        ax.legend(loc=legend_location, frameon=False, fontsize=legend_size)
        ax.set_xlabel(x_label, fontsize=x_label_size)
        ax.set_ylabel(y_label, fontsize=y_label_size)
        
        plt.show()