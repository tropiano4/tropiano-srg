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
#   Testing overall factor difference between SNMD and AV18 data.


import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
from Figures import figures_functions as ff
from lda import load_density, LDA
import operators as op
from Potentials.vsrg_macos import vnn
from snmd import single_nucleon_momentum_distributions
from SRG.srg_unitary_transformation import SRG_unitary_transformation


# Copied from operator_evolution_fig.ipynb
def snmd_tails_with_AV18(nucleus, channels, kvnn, lamb, kmax=10.0, kmid=2.0,
                         ntot=120, xlim=(2, 6), ylim=(1e-2, 1e4)):
    """
    Nuclear-averaged and SRG-evolved proton momentum distributions
    n_\lambda^A(q) / Z. This figure compares our LDA momentum distributions
    with AV18 data.
    
    Parameters
    ----------
    nucleus : tuple
        Details for various nuclei formatted as a tuple: (name (str), Z (int),
        N (int)) (e.g., ('O16', 8, 8)).
    channels : tuple
        Partial wave channels to include in the calculation
        (e.g., ('1S0', '3S1')).
    kvnn : int
        This number specifies the potential.
    lamb : float
        SRG \lambda parameter [fm^-1].
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh.
    xlim : tuple, optional
        Limits of x-axis [fm^-1].
    ylim : tuple, optional
        Limits of y-axis [fm^3].
        
    Returns
    -------
    f : Figure
        Figure object from matplotlib subplots function.
    ax : axes.Axes object
        Single Axes object from matplotlib subplots function.
    d : dict
        Dictionary for momentum distribution, momentum nodes and weights.
    
    """
    
    # Load data from the following directory
    data_directory = 'Figures/SRC_physics/Data'
    
    # Name of nucleus (e.g., 'Ca40')
    nucleus_name = nucleus[0]
    Z = nucleus[1]
    N = nucleus[2]
    
    d = {}
    
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
        
    # Initialize pair momentum distributions class for given potential
    snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb,
                                                     kmax, kmid, ntot)
    
    # Load r values and nucleonic densities (the r_array's are the same)
    r_array, rho_p_array = load_density(nucleus_name, 'proton', Z, N)
    r_array, rho_n_array = load_density(nucleus_name, 'neutron', Z, N)
    
    # Call LDA class
    lda = LDA(r_array, rho_p_array, rho_n_array)
    
    # Calculate nuclear-averaged momentum distributions
    snmd_data = lda.local_density_approximation(q_array, snmd.n_lambda,
                                                    'p', 'q_contributions')
        
    n_p_total_array = snmd_data[:, 0]
    n_p_1_array = snmd_data[:, 1]
    n_p_delU_array = snmd_data[:, 2]
    n_p_delU2_array = snmd_data[:, 3]
        
    d['snmd_data'] = snmd_data
    d['momentum'] = q_array
    d['weights'] = q_weights

    # Figure size
    row_number = 1
    col_number = 1
    figure_size = (4*col_number, 4*row_number)

    # Axes labels and fontsize
    x_label = 'q [fm' + r'$^{-1}$' + ']'
    x_label_size = 16
    y_label = 'proton ' + r'$n^{\lambda}_A(q)$' + ' [fm' + r'$^3$' + ']'
    y_label_size = 16
    
    # Curve width
    curve_width = 2.0

    # Initialize figure
    plt.close('all')
    f, ax = plt.subplots(figsize=figure_size)

    # Legend label
    # curve_label = ff.nuclei_label_conversion(nucleus_name) # Labels the nucleus
        
    # Add curve to figure
    ax.semilogy(q_array, n_p_total_array, color='xkcd:grey', label='total',
                linewidth=curve_width, linestyle='solid')
    ax.semilogy(q_array, n_p_1_array, color='xkcd:blue', label='1',
                linewidth=curve_width, linestyle='dotted')
    ax.semilogy(q_array, abs(n_p_delU_array), color='xkcd:green',
                label=r'$|\delta U|$', linewidth=curve_width,
                linestyle='dashed')
    ax.semilogy(q_array, n_p_delU2_array, color='xkcd:red',
                label=r'$\delta U \delta U^{\dagger}$', linewidth=curve_width,
                linestyle='dashdot')
        
    # Add AV18 data with error bars
    av18_data = np.loadtxt(data_directory + '/' + 'AV18_%s_snmd.txt' % nucleus_name)
    q_array_av18 = av18_data[:, 0] # fm^-1
    n_p_array_av18 = av18_data[:, 1]
    error_bars_array_av18 = av18_data[:, 2]
            
    # AV18 data with error bars
    ax.set_yscale('log')
    ax.errorbar(q_array_av18, n_p_array_av18, yerr=error_bars_array_av18,
                color='xkcd:black', label='AV18', linestyle='', marker='.')

    # Specify axes limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
        
    # Set axes labels and legend
    ax.set_xlabel(x_label, fontsize=x_label_size)
    ax.set_ylabel(y_label, fontsize=y_label_size)

    return f, ax, d


if __name__ == '__main__':
    
    # Compare O16 proton momentum distribution to AV18
    
    # Further test
    #   Mesh dependence
    #   - Very little mesh dependence
    #   Higher partial waves?
    #   Dependence on K mesh
    #   Dependence on x mesh

    nucleus =  ('O16', 8, 8)
    # channels =  ('1S0', '3S1')
    channels = ('1S0', '3S1', '3P0', '1P1', '3P1')
    # channels = ('1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3')
    # Bugs in 3P2, 3D3?
    kvnn = 6
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120 # Default mesh
    # kmax, kmid, ntot = 30.0, 4.0, 120

    # # Use data
    # xlim = (1.5, 6)
    # ylim = (1e-3, 1e2)
    # f, ax = snmd_tails_with_AV18(nucleus, channels, kvnn, lamb, kmax, kmid,
    #                               ntot, xlim, ylim)
    
    # Calculate
    xlim = (0, 6)
    ylim = (1e-3, 1e4)
    f, ax, d = snmd_tails_with_AV18(nucleus, channels, kvnn, lamb, kmax, kmid,
                                    ntot, xlim, ylim)
    
    # Add legend
    legend_size = 16
    legend_location = 'upper right'
    ax.legend(loc=legend_location, frameon=False, fontsize=legend_size)
    
    # Print normalization
    q_array = d['momentum']
    q_weights = d['weights']
    snmd_data = d['snmd_data']
    n_p_total_array = snmd_data[:, 0]
    n_p_1_array = snmd_data[:, 1]
    n_p_delU_array = snmd_data[:, 2]
    n_p_delU2_array = snmd_data[:, 3]
    
    # Overall factor
    factor = 4*np.pi * 1/(2*np.pi)**3
    
    # Total normalization
    norm_total = factor*np.sum(q_weights * q_array**2 * n_p_total_array)
    print(norm_total)
    norm_1 = factor*np.sum(q_weights * q_array**2 * n_p_1_array)
    print(norm_1)
    norm_delU = factor*np.sum(q_weights * q_array**2 * n_p_delU_array)
    print(norm_delU)
    norm_delU2 = factor*np.sum(q_weights * q_array**2 * n_p_delU2_array)
    print(norm_delU2)
    
    # # Normalization above 2 fm^-1
    # q_index = op.find_q_index(2, q_array)
    # norm_2 = np.sum( (q_weights * q_array**2 * n_p_array)[q_index:] ) * \
    #          overall_factor
    # print(norm_2)
    
    
    # --- Test deuteron contribution normalizations --- #
    
    # from dmd import deuteron_momentum_distributions
    
    # # Load momentum
    # q_array, q_weights = vnn.load_momentum(kvnn, '3S1', kmax, kmid, ntot)
    
    # # Initialize deuteron momentum distributions class for given potential
    # dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)

    # # Loop over momenta q
    # total_array = np.zeros(ntot)
    # term_1_array = np.zeros(ntot)
    # term_deltaU_array = np.zeros(ntot)
    # term_deltaU2_array = np.zeros(ntot)
    # for iq, q in enumerate(q_array):
        
    #     # Calculate each contribution to n_\lambda(q) including total
    #     total, term_1, term_deltaU, term_deltaU2 = dmd.n_lambda_pair_exact(q,
    #                                                'q_contributions')
        
    #     total_array[iq] = total
    #     term_1_array[iq] = term_1
    #     term_deltaU_array[iq] = term_deltaU
    #     term_deltaU2_array[iq] = term_deltaU2
        
    # # Print normalizations
    # integration_measure = 2/np.pi * q_array**2 * q_weights
    # total_norm = np.sum(total_array*integration_measure)
    # term_1_norm = np.sum(term_1_array*integration_measure)
    # term_deltaU_norm = np.sum(term_deltaU_array*integration_measure)
    # term_deltaU2_norm = np.sum(term_deltaU2_array*integration_measure)
    # print(total_norm)
    # print(term_1_norm)
    # print(term_deltaU_norm)
    # print(term_deltaU2_norm)
        

    # --- Test interpolation routine from SNMD --- #
    
    # # Load momentum
    # q_array, q_weights = vnn.load_momentum(kvnn, '3S1', kmax, kmid, ntot)
    # # For dividing out momenta/weights
    # factor_array = np.sqrt( (2*q_weights) / np.pi ) * q_array
    # # For coupled-channel matrices
    # factor_array_cc = np.concatenate( (factor_array, factor_array) )
    
    # # Get \delta U matrix here for 1S0
    # # channel = '1S0'
    # # channel = '3S1'
    # channel = '3P2'
    # # Load SRG transformation
    # H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    # H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
    #                                   method='srg', generator='Wegner',
    #                                   lamb=lamb)
    # # Load U(k, k') [unitless]
    # U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
    
    # # Isolate 2-body term and convert to fm^3
    # if vnn.coupled_channel(channel):
    #     I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
    #     row, col = np.meshgrid(factor_array_cc, factor_array_cc)
    # else:
    #     I_matrix_unitless = np.eye(ntot, ntot)
    #     row, col = np.meshgrid(factor_array, factor_array)

    # delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
    # delta_U_matrix = delta_U_matrix_unitless / row / col
    
    # # Calculate < k | \delta U \delta U^\dagger | k' >
    # if vnn.coupled_channel(channel):
    #     deltaU2_matrix = delta_U_matrix[:ntot, :ntot]**2 + \
    #                       delta_U_matrix[:ntot, ntot:]**2 + \
    #                       delta_U_matrix[ntot:, :ntot]**2 + \
    #                       delta_U_matrix[ntot:, ntot:]**2
    # else:
    #     deltaU2_matrix = delta_U_matrix**2
    
    # # Interpolate < k | \delta U \delta U^\dagger | k' >
    # deltaU2_func = RectBivariateSpline(q_array, q_array, deltaU2_matrix)
    
    # # Evaluate < k | \delta U \delta U^\dagger | q > for some q
    # # Note, it would be q - K/2 in the code
    # q_value = 2.5
    # q_index = op.find_q_index(q_value, q_array)
    # deltaU2_k = deltaU2_func.ev(q_array, q_value)
    # # print(deltaU2_k)
    # print(deltaU2_matrix[:ntot, q_index])
    
    # # Compare to < k | \delta U | q > < q | \delta U^\dagger | k >
    # if vnn.coupled_channel(channel):
    #     deltaU2_k_exact = delta_U_matrix[:ntot, q_index]*delta_U_matrix.T[q_index, :ntot]+\
    #                       delta_U_matrix[:ntot, ntot+q_index]*delta_U_matrix.T[ntot+q_index, :ntot]+\
    #                       delta_U_matrix[ntot:, q_index]*delta_U_matrix.T[q_index, ntot:]+\
    #                       delta_U_matrix[ntot:, ntot+q_index]*delta_U_matrix.T[ntot+q_index, ntot:]
    # else:
    #     deltaU2_k_exact = delta_U_matrix[:ntot, q_index] * delta_U_matrix.T[q_index, :ntot]
    # print(deltaU2_k_exact)
    # print(deltaU2_matrix[:ntot, q_index]-deltaU2_k_exact)