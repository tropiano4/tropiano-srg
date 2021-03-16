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
#   01/26/21 --- Renamed to test_script.py.
#
#------------------------------------------------------------------------------

# Description of this test:
#   Calculate n_{\lambda}(q) for different nuclei using LDA with simple
#   filled Fermi sea single-particle momentum distributions. Split into pp,
#   pn, nn, np contributions as in Ryckebusch_2019oya.


import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from lda import load_density, LDA
from Potentials.vsrg_macos import vnn


def n_single_particle(q, k_F):
        
    # This should just be \theta(k_F-p)
        
    if q <= k_F:
            
        return 1
        
    else:
            
        return 0


if __name__ == '__main__':

    
    # --- Set up --- #
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    # kmax, kmid, ntot = 10.0, 2.0, 120
    kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2

    # Details of example nuclei (format is [nuclei, Z, N])
    # nuclei_list = ['O16', 8, 8]
    nuclei_list = ['Ca40', 20, 20]
    # nuclei_list = ['Ca48', 20, 28]
    # nuclei_list = ['Pb208', 82, 126]
    nucleus = nuclei_list[0]
    Z = nuclei_list[1]
    N = nuclei_list[2]
    A = N + Z
    r_array, rho_p_array = load_density(nucleus, 'proton', Z, N)
    r_array, rho_n_array = load_density(nucleus, 'neutron', Z, N)
    
    # Call LDA class
    lda = LDA(r_array, rho_p_array, rho_n_array)
    
    
    # --- Calculate momentum distributions --- #

    n_infty_p_exp_array = lda.local_density_approximation( q_array, 
                          n_single_particle, 'p' )
    n_infty_n_exp_array = lda.local_density_approximation( q_array, 
                          n_single_particle, 'n' )
    
    overall_factor = (4*np.pi)**2 * 2 * 1/(2*np.pi)**3
    
    p_a_ipm_p = q_array**2 * n_infty_p_exp_array / A * overall_factor
    p_a_ipm_n = q_array**2 * n_infty_n_exp_array / A * overall_factor

    # Split into four pp, pn, nn, np pieces
    p_a_ipm_pp = p_a_ipm_p * (Z-1)/(A-1)
    p_a_ipm_pn = p_a_ipm_p * N/(A-1)
    p_a_ipm_np = p_a_ipm_n * Z/(A-1)
    p_a_ipm_nn = p_a_ipm_n * (N-1)/(A-1)
    
    p_a_ipm = p_a_ipm_pp + p_a_ipm_pn + p_a_ipm_np + p_a_ipm_nn
    
    # This is 4 \pi for d3p, 2 for spin sum, and 1/(2*\pi)^3 for converting
    # from sums to integrals
    print(np.sum(p_a_ipm*q_weights))
    
    
    # --- Plot pair momentum distributions --- #
    
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )
    plt.semilogy(q_array, p_a_ipm_pp, color='red', linestyle='dashdot', 
                  label=nucleus+' (pp)')
    plt.semilogy(q_array, p_a_ipm_nn, color='blue', linestyle='dashed', 
                  label=nucleus+' (nn)')
    plt.semilogy(q_array, p_a_ipm_pn+p_a_ipm_np, color='green', 
                 linestyle='dotted', label=nucleus+' (pn+np)') 
    plt.semilogy(q_array, p_a_ipm, color='black', label=nucleus+' (total)')
    plt.ylabel(r'$P^A(p)$')
    ax.set_xlim( [min(q_array), 4] )
    ax.set_ylim( [1e-3, 2e0])
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('p [fm' + r'$^{-1}$' + ']')
    
    plt.show()