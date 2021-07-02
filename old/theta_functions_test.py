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
# extension _test.py. For example, momentum_projection_operator_testv1.py (v1
# means 'version 1'). Use the revision history below to document when and why
# these files are created.
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
#   04/28/21 --- Testing normalization and contributions of \delta U, etc. or
#                pp/pn to single-nucleon momentun distributions. Created
#                lda_normalizations_test.py in Old_codes.
#   05/04/21 --- Testing higher partial waves of SRG transformations: 3P2-3F2
#                and 3D3-3G3 have numerical artifacts. Created
#                high_partial_waves_srg_test.py in Old_codes.
#   06/10/21 --- Verifying \theta functions averaging in snmd.py and dmd.py by
#                comparing numerical functions to analytic evaluation of 
#                \int d3K \int d3k \theta(kF-|K/2+k|) \theta(kF-|K/2-k|).
#                Created theta_functions_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Test evaluation of F_2(Q,k) in paper compared to analytic result.
#   Integration of \int d3K \int d3k \theta(kF-|K/2+k|) \theta(kF-|K/2-k|)
#   where k and K are vectors should give (4*\pi)^2 kF^6 / 9.


import numpy as np
# Scripts made by A.T.
from Misc.integration import gaussian_quadrature_mesh


def select_number_integration_points(k_max, k_min=0.0):

    # Interval of integration
    interval = k_max - k_min
    
    # Basing these numbers off expected kF values
    if interval >= 1.2:
        ntot_k = 60
    elif 1.2 > interval >= 1.0:
        ntot_k = 50
    elif 1.0 > interval >= 0.8:
        ntot_k = 40
    elif 0.8 > interval >= 0.6:
        ntot_k = 30
    elif 0.6 > interval >= 0.4:
        ntot_k = 20
    else:
        ntot_k = 10
            
    return ntot_k


def theta_deltaU2_diff(kF_1, kF_2, K, k_array, ntot_k):
    
    # Make \theta( k_F - \abs(K_vec/2 +(-) k_vec) ) the same length as
    # k_array
    theta_deltaU2 = np.zeros(ntot_k)
        
    # Loop over each momenta k and go through the four cases
    for i, k in enumerate(k_array):
            
        # Case 1: 2k+K < 2kF_1 and 2k+K < 2kF_2
        if 2*k+K <= 2*kF_1 and 2*k+K <= 2*kF_2:
            theta_deltaU2[i] = 1
                
        # Case 2: 2k+K > 2kF_1 and 2k+K > 2kF_2 and 
        # 4k^2+K^2 < 2(kF_1^2+kF_2^2)
        elif 2*k+K > 2*kF_1 and 2*k+K > 2*kF_2 and \
             4*k**2 + K**2 <= 2*(kF_1**2 + kF_2**2):
            theta_deltaU2[i] = ( 2*(kF_1**2 + kF_2**2) - 4*k**2 - K**2 ) \
                               / (4*k*K)
                            
        # Case 3: 2k+K < 2kF_2 and -4 < (4k^2 - 4kF_1^2 + K^2)/(kK) < 4
        elif 2*k+K <= 2*kF_2 and -4 < (4*k**2-4*kF_1**2+K**2)/(k*K) <= 4:
            theta_deltaU2[i] = ( 4*kF_1**2 - (K-2*k)**2 ) / (8*k*K)
                
        # Case 4: 2k+K < 2kF_1 and -4 < (4k^2 - 4kF_2^2 + K^2)/(kK) < 4
        elif 2*k+K <= 2*kF_1 and -4 < (4*k**2-4*kF_2**2+K**2)/(k*K) <= 4:
            theta_deltaU2[i] = ( 4*kF_2**2 - (K-2*k)**2 ) / (8*k*K)
                
        # Otherwise, F(K,k) = 0
            
    return theta_deltaU2


def theta_deltaU2_same(kF, K, k_array, ntot_k):
    
    # Make \theta( k_F - \abs(K_vec/2 +(-) k_vec) ) the same length as k_array
    theta_deltaU2 = np.zeros(ntot_k)
        
    # Loop over each momenta k and go through the two inequalities
    for i, k in enumerate(k_array):
                
        # Case 1: k < kF-K/2 F(K,k) = 1
        if k < kF-K/2:
            theta_deltaU2[i] = 1
                
        # Case 2: kF-K/2 < k < \sqrt(kF^2-K^2/4)
        # -> F(K,k) = ( kF^2 - k^2 - K^2/4 ) / (k*K)
        elif kF-K/2 < k < np.sqrt(kF**2-K**2/4):
            theta_deltaU2[i] = ( kF**2 - k**2 - K**2/4 ) / ( k*K )

        # Otherwise, k > \sqrt(kF^2-K^2/4) and F(K,k) = 0
                
    return theta_deltaU2


def integrand_K(kF_1, kF_2, K, k_array, k_weights, ntot_k, case='snmd'):
    
    if case == 'snmd':
        theta_k = theta_deltaU2_diff(kF_1, kF_2, K, k_array, ntot_k)
    elif case == 'dmd':
        kF = kF_1
        theta_k = theta_deltaU2_same(kF, K, k_array, ntot_k)
    integration_measure = 4*np.pi*k_weights*k_array**2
    return np.sum(theta_k*integration_measure)


def integral(kF_1, kF_2, K_array, K_weights, ntot_K, case='snmd'):
    
    integrand = np.zeros(ntot_K)
    for iK, K in enumerate(K_array):
        
        kF_min = min(kF_1, kF_2)
        kmin = max(K/2 - kF_min, 0)
        kmax = min( np.sqrt( ( kF_1**2 + kF_2**2 )/2 - K**2/4 ), kF_min + K/2 )

        # Select number of integration points based on kmax_delU2
        ntot_k = select_number_integration_points(kmax, kmin)

        # Get Gaussian quadrature mesh
        k_array, k_weights = gaussian_quadrature_mesh(kmax, ntot_k, xmin=kmin)
        
        integrand[iK] = integrand_K(kF_1, kF_2, K, k_array, k_weights, ntot_k,
                                    case)
        
    integrand *= 4*np.pi*K_weights*K_array**2
    return np.sum(integrand)
    
if __name__ == '__main__':
    
    # kF = 0.1
    # kF = 0.5
    # kF = 1.0
    kF = 1.3
    kF_1 = kF
    kF_2 = kF
    # case = 'snmd'
    case = 'dmd'
    
    # Create meshes
    Kmax = kF_1 + kF_2
    # Total number of points
    if Kmax >= 2.0:
        Ntot = 50
    elif 2.0 > Kmax >= 1.6:
        Ntot = 40
    elif 1.6 > Kmax >= 1.2:
        Ntot = 30
    elif 1.2 > Kmax >= 0.8:
        Ntot = 20
    else:
        Ntot = 10
    K_array, K_weights = gaussian_quadrature_mesh(Kmax, Ntot)
    
    value = integral(kF_1, kF_2, K_array, K_weights, Ntot, case)
    
    exact = (4*np.pi)**2 * kF**6/9
    print('Numerical = %.5f' % value )
    print('Exact = %.5f' % exact)