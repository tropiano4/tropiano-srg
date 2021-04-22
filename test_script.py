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
# Scripts made by A.T.
from Figures import figures_functions as ff
from lda import load_density, LDA
from Potentials.vsrg_macos import vnn
from snmd import single_nucleon_momentum_distributions


# Copied from operator_evolution_fig.ipynb
def snmd_tails_with_AV18(nucleus, channels, kvnn, lamb, kmax=10.0, kmid=2.0,
                         ntot=120, xlim=(2, 6), ylim=(1e-2, 1e4), calc=False):
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
    
    """
    
    # Load data from the following directory
    data_directory = 'Figures/SRC_physics/Data'
    
    # Name of nucleus (e.g., 'Ca40')
    nucleus_name = nucleus[0]
    Z = nucleus[1]
    N = nucleus[2]
    
    if calc:
        
        q_array, _ = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
        
        # Initialize pair momentum distributions class for given potential
        snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb,
                                                     kmax, kmid, ntot)
    
        # Load r values and nucleonic densities (the r_array's are the same)
        r_array, rho_p_array = load_density(nucleus_name, 'proton', Z, N)
        r_array, rho_n_array = load_density(nucleus_name, 'neutron', Z, N)
    
        # Call LDA class
        lda = LDA(r_array, rho_p_array, rho_n_array)
    
        # Calculate nuclear-averaged momentum distributions
        n_p_array = lda.local_density_approximation(q_array, snmd.n_lambda, 'p')
        
    else:
        # for now scale up by this factor (based on LDA / AV18 at q=4.1 fm^-1)
        # factor =  1.18700e-1 / 1.25227e-7
    
        # Get rid of 1/(2\pi)^6 (2/\pi)^2
        # factor = (2*np.pi)**6 * (np.pi/2)**2
    
        # include factor of 2^4 for four (1 - (-1)^(L+S+T)) terms
        # factor = (2*np.pi)**6 * (np.pi/2)**2 * 2**4
    
        # keep 2/\pi terms
        factor = (2*np.pi)**6 * 2**4
        
        # Data file name
        file_name = 'lda_snmd_%s_channels' % nucleus_name
        # Add each channel to file name
        for channel in channels:
            file_name += '_%s' % channel
        file_name += '_kvnn_%d_lamb_%.2f_kmax_%.1f' % (kvnn, lamb, kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        # Load momentum and single-nucleon momentum distributions
        data = np.loadtxt(data_directory + '/' + file_name)
        # Momentum in fm^-1
        q_array = data[:, 0]
        # Proton momentum distribution scaled up by overall factor
        n_p_array = data[:, 1] * factor

    # Figure size
    row_number = 1
    col_number = 1
    figure_size = (4*col_number, 4*row_number)

    # Axes labels and fontsize
    x_label = 'q [fm' + r'$^{-1}$' + ']'
    x_label_size = 16
    y_label = 'proton ' + r'$n^{\lambda}_A(q)/Z$' + ' [fm' + r'$^3$' + ']'
    y_label_size = 16
    
    # Curve width
    curve_width = 2.0

    # Initialize figure
    plt.close('all')
    f, ax = plt.subplots(figsize=figure_size)

    # Legend label
    curve_label = ff.nuclei_label_conversion(nucleus_name) # Labels the nucleus
        
    # Add curve to figure
    ax.semilogy(q_array, n_p_array, color='xkcd:red', label=curve_label,
                linewidth=curve_width)
        
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

    return f, ax


if __name__ == '__main__':
    
    # Compare O16 proton momentum distribution to AV18

    nucleus =  ('O16', 8, 8)
    channels =  ('1S0', '3S1')
    kvnn = 6
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120 # Default mesh

    # Use data
    xlim = (1.5, 6)
    ylim = (1e-3, 1e2)
    f, ax = snmd_tails_with_AV18(nucleus, channels, kvnn, lamb, kmax, kmid,
                                  ntot, xlim, ylim)
    
    # # Calculate
    # xlim = (0, 6)
    # ylim = (1e-3, 1e4)
    # f, ax = snmd_tails_with_AV18(nucleus, channels, kvnn, lamb, kmax, kmid,
    #                               ntot, xlim, ylim, calc=True)
    

    # Add legend
    legend_size = 16
    legend_location = 'upper right'
    ax.legend(loc=legend_location, frameon=False, fontsize=legend_size)