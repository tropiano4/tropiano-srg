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
#   03/16/21 --- Testing pair momentum distributions for different nuclei
#                using LDA with simple filled Fermi sea single-particle
#                momentum distributions. Created pmd_simple_test.py in
#                Old_codes.
#   03/16/21 --- Testing pair momentum distribution for deuteron starting from
#                second quantization derivation and using expansion of
#                U(k, k'). Created pmd_deuteron_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Calculate n_{\lambda}(q) for single-nucleons using LDA.


import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline
import time
# Scripts made by A.T.
from lda import load_density, LDA
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


# --- Fix --- #
# 1. You can speed this up by creating full integrand and doing one dx
#    integral, then dk integral, etc. (Didn't work.)
# 2. Seems to something happening in last q_array point.
# 3. Fix the normalization business. \int dq P^A(q) = A_Ca40/A_Ca48?

# --- Checks --- #
# 1. What is the sensitivity to K_mesh details?
#    K_max = 2-4 fm^-1 doesn't change.
#    15-30 points doesn't change.
# 2. Are you doing the integration over K correctly?
# 3. What is the sensitivity to x_mesh details?
#    10 -> 6 points doesn't change.
# 4. Are you doing the integration over angles correctly?
# 5. Check long expressions for typos (coupled-channel for sure!)


def construct_K_mesh(K_max=3.0, ntot=30):
    """
    Creates a momentum mesh for COM mometum.

    Parameters
    ----------
    K_max : float
        Maximum value in the momentum mesh [fm^-1].
    ntot : int
        Number of momentum points in mesh.

    Returns
    -------
    K_array : 1-D ndarray
        Momentum array [fm^-1].
    K_weights: 1-D ndarray
        Momentum weights [fm^-1].

    """
    
    # Minimum momentum value
    K_min = 0.0

    x_array, x_weights = leggauss(ntot) # Interval [-1,1]
    
    # Convert from interval [-1, 1] to [a, b] (meaning x_array -> k_array)
    K_array = 0.5 * (x_array + 1) * (K_max - K_min) + K_min
    K_weights = (K_max - K_min) / 2 * x_weights

    return K_array, K_weights


class single_nucleon_momentum_distributions(object):
    
    
    def __init__(self, kvnn, lamb, kmax=0.0, kmid=0.0, ntot=0):
        """
        Description.

        Parameters
        ----------
        kvnn : TYPE
            DESCRIPTION.
        lamb : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        channels = ['1S0', '3S1']
        
        d = {}
        
        for channel in channels:
            
            # Store relevant pieces in dictionary for each channel
            d[channel] = {}
            
            # Load momentum mesh and store in dictionary
            k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid,
                                                   ntot)
            ntot = len(k_array)
            d[channel]['k_array'] = k_array
            d[channel]['k_weights'] = k_weights
            d[channel]['ntot'] = ntot
            
            # Load SRG transformation in 1S0 and 3S1-3D1 channels
            H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
            H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
                                             method='srg', generator='Wegner', 
                                             lamb=lamb)
            # Load U(k, k') [unitless]
            U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
            
            # Divide out momentum/weights factor
            factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
            if channel == '3S1':
                factor_array = np.concatenate( (factor_array, factor_array) )
            row, col = np.meshgrid(factor_array, factor_array)
            
            # Converting to units [fm^3]
            I_matrix_unitless = np.eye( len(factor_array), len(factor_array) )
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            
            # Store both matrices in the dictionary
            # d[channel]['I'] = I_matrix_unitless
            # d[channel]['delta_U'] = delta_U_matrix_unitless
            
            # U_matrix = U_matrix_unitless / row / col
            I_matrix = I_matrix_unitless / row / col
            delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # Store both matrices in the dictionary
            d[channel]['I'] = I_matrix # [fm^3]
            d[channel]['delta_U'] = delta_U_matrix # [fm^3]
            
            # Interpolate deltaU matrices for evaluation of fourth term in
            # single-particle momentum distributions
        
            # 1S0 channel
            if channel == '1S0':
                d['deltaU_1S0_func'] = RectBivariateSpline(k_array, k_array,
                                       delta_U_matrix)
                d['deltaU_1S0_dag_func'] = RectBivariateSpline(k_array, k_array,
                                           delta_U_matrix.T)
            # 3S1-3D1 channel
            else:
                d['deltaU_3S1-3S1_func'] = RectBivariateSpline(k_array, k_array,
                                           delta_U_matrix[:ntot, :ntot] )
                d['deltaU_3S1-3S1_dag_func'] = RectBivariateSpline(k_array, k_array,
                                               delta_U_matrix.T[:ntot, :ntot] )
                d['deltaU_3S1-3D1_func'] = RectBivariateSpline(k_array, k_array,
                                           delta_U_matrix[:ntot, ntot:] )
                d['deltaU_3D1-3S1_dag_func'] = RectBivariateSpline(k_array, k_array,
                                               delta_U_matrix.T[ntot:, :ntot] )
                # d['deltaU_3D1-3D1_func'] = RectBivariateSpline(k_array, k_array,
                #                              delta_U_matrix[ntot:, ntot:] )
            
        # Create mesh for integration over COM momentum
        # K_array, K_weights = construct_K_mesh() # K_max = 3 fm^-1, 30 points
        # K_array, K_weights = construct_K_mesh(ntot=20) # K_max = 3 fm^-1, 20 points
        K_array, K_weights = construct_K_mesh(ntot=15) # K_max = 3 fm^-1, 15 points
        Ktot = len(K_array)
        d['Ktot'] = Ktot
        d['K_array'] = K_array
        d['K_weights'] = K_weights
            
        # Create mesh for integration over angles
        xtot = 10
        # xtot = 6 # No change
        x_array, x_weights = leggauss(xtot) # Interval [-1,1]
        d['xtot'] = xtot
        d['x_array'] = x_array
        d['x_weights'] = x_weights
            
        # Save dictionary
        self.d = d


    def theta_q_2k(self, k_F, q):
        # Evaluates \theta( k_F - \abs( q_vec - 2k_vec ) ) for every momentum 
        # k and angle x
        
        k_array = self.d['1S0']['k_array']
        x_array = self.d['x_array']
        x_weights = self.d['x_weights']
        
        k_grid, x_grid = np.meshgrid(k_array, x_array, indexing='ij')
        q_2k_magnitude = q**2 + 4*k_grid**2 - 4*q*k_grid*x_grid
        
        # This returns a (ntot, xtot) array of boolean values at every point
        # and converts to 1's or 0's by multiplying 1
        theta_q_2k = ( q_2k_magnitude < k_F ) * 1
        
        # # Test this first by comparing to double loop
        # # It works!
        # ntot = self.d['1S0']['ntot']
        # xtot = self.d['xtot']
        # theta_q_2k_test = np.zeros( (ntot, xtot) )
        # for i, k in enumerate(k_array):
        #     for j, x in enumerate(x_array):
        #         theta_q_2k_test[i, j] = int ( ( q**2 + 4*k**2 - 4*q*k*x ) < k_F )
                
        # if theta_q_2k.all() == theta_q_2k_test.all():
        #     print('True')
        # else:
        #     print('False')
        
        # Use this to weight the matrix with dx for integration purposes
        k_grid, x_weights_grid = np.meshgrid(k_array, x_weights, indexing='ij')
        theta_q_2k_weights = theta_q_2k * x_weights_grid

        return theta_q_2k_weights
    
    
    def theta_K_k(self, k_F, sign=1):
        # Evaluates \theta( k_F - \abs( 1/2*K_vec +/- k_vec ) ) for every
        # COM momentum K, relative momentum k, and angle x where sign
        # specifies the sign of k_vec
        
        k_array = self.d['1S0']['k_array']
        K_array = self.d['K_array']
        x_array = self.d['x_array']
        x_weights = self.d['x_weights']
        
        k_grid, K_grid, x_grid = np.meshgrid(k_array, K_array, x_array,
                                             indexing='ij')
        kK_magnitude = K_grid**2/4 + k_grid**2 + sign*k_grid*K_grid*x_grid
        # This returns a (ntot, Ktot, xtot) array of boolean values at every
        # point and converts to 1's or 0's by multiplying 1
        theta_K_k = ( kK_magnitude < k_F ) * 1
        
        # Use this to weight the matrix with dx for integration purposes
        k_grid, K_grid, x_weights_grid = np.meshgrid(k_array, K_array,
                                                     x_weights, indexing='ij')
        theta_K_k_weights = theta_K_k * x_weights_grid
        
        return theta_K_k_weights
            
    
    def n_lambda_N(self, q, kF_1, kF_2):
        # don't forget docstring
        # kF_1: Fermi momentum for nucleon N
        # kF_2 : Fermi momentum for opposite nucleon (so if you specified
        # a proton distribution then Np refers to a neutron)
        
        # Load momentum mesh and SRG transformation for 1S0 and 3S1-3D1 
        # channels
        
        # The momentum mesh is the same for 1S0 and 3S1-3D1
        k_array = self.d['1S0']['k_array']
        k_weights = self.d['1S0']['k_weights']
        ntot = self.d['1S0']['ntot']
        
        # Momentum mesh for COM momentum
        K_array = self.d['K_array']
        K_weights = self.d['K_weights']
        
        # x_array = self.d['x_array']
        
        deltaU_1S0 = self.d['1S0']['delta_U'] # fm^3
        deltaU_3S1 = self.d['3S1']['delta_U'] # fm^3
        
        # Common factors
        integration_k_measure = 2/np.pi * k_weights * k_array**2

        # Split into low- and high-q terms
        
        # Low-q terms
        # Term 1: 1 * n(q) * 1
        # Term 2: \deltaU * n(q) * 1
        # Term 3: 1 * n(q) * \deltaU^\dagger = Term 2
        
        # The first three terms have \theta(kF_1 - q)
        if q < kF_1:
            
            first_term = 2
            
            # These are (ntot, xtot) matrices of 
            # \theta( kF1 - \abs( q_vec - 2k_vec ) ) for every momentum k and 
            # angle x
            theta_kF1_k_x_matrix = self.theta_q_2k(kF_1, q)
            theta_kF2_k_x_matrix = self.theta_q_2k(kF_2, q)
            
            # Do integration over x first (angle-averaging) where the 
            # integration weights of x are already built-in
            # This sum collapses theta_k_x_matrix to a vector dependent only
            # on k: (ntot, xtot) -> (ntot, 1)
            theta_kF1_k_vector = np.sum( theta_kF1_k_x_matrix, axis=-1 ) / 2
            theta_kF2_k_vector = np.sum( theta_kF2_k_x_matrix, axis=-1 ) / 2
            
            # Build integrand for k integration
            integrand_k = integration_k_measure * ( np.diag( deltaU_1S0 ) *\
                          ( theta_kF1_k_vector + theta_kF2_k_vector/2 ) + \
                          3/2 * np.diag( deltaU_3S1[:ntot, :ntot] ) * \
                          theta_kF2_k_vector )
            
            # # Build integrand with respect to dk dx
            # deltaU_1S0_grid, _ = np.meshgrid( np.diag( deltaU_1S0 ), x_array,
            #                                   indexing='ij' )
            # deltaU_3S1_grid, _ = np.meshgrid( np.diag( deltaU_3S1 )[:ntot],
            #                                   x_array, indexing='ij' )
            # # This is a (ntot, xtot matrix)
            # integrand_k_x = deltaU_1S0_grid * ( theta_kF1_k_x_matrix + \
            #                 theta_kF2_k_x_matrix/2 ) + 3/2 * deltaU_3S1_grid * \
            #                 theta_kF2_k_x_matrix
            # # Do integration over x first (angle-averaging) where the 
            # # integration weights of x are already built-in
            # integrand_k = integration_k_measure * np.sum( integrand_k_x/2,
            #                                               axis=-1 )
            
            # Do integration over k now where the factor of 2 is for combining
            # the second and third terms
            middle_terms = 2 * np.sum(integrand_k)
            
        # q > kF_1
        else:
            
            first_term = 0
            middle_terms = 0
        
        # High-q term
        
        # Load interpolated deltaU's
        deltaU_1S0_func = self.d['deltaU_1S0_func']
        deltaU_1S0_dag_func = self.d['deltaU_1S0_dag_func']
        deltaU_3S1_3S1_func = self.d['deltaU_3S1-3S1_func']
        deltaU_3S1_3S1_dag_func = self.d['deltaU_3S1-3S1_dag_func']
        deltaU_3S1_3D1_func = self.d['deltaU_3S1-3D1_func']
        deltaU_3D1_3S1_dag_func = self.d['deltaU_3D1-3S1_dag_func']
        
        q_K_array = np.sqrt( q**2 + K_array**2/4 )
        # We can only interpolate \deltaU(k,k') up to k_max
        # Do not allow this next array to go above k_max - this may impose
        # a cutoff on the K integration
        q_K_array_cutoff = q_K_array[ q_K_array < max(k_array) ]
        K_cutoff_index = len(q_K_array_cutoff)
        
        k_grid, q_K_grid, = np.meshgrid(k_array, q_K_array_cutoff,
                                        indexing='ij')
        # k_grid, q_K_grid, _ = np.meshgrid(k_array, q_K_array_cutoff, x_array,
        #                                   indexing='ij')
        
        # Evaluate deltaU at (k, \abs(q_vec-K_vec/2)) where we approximate
        # \abs(q_vec-K_vec/2) ~ q_vec^2 + K_vec^2/4
        deltaU_squared_1S0 = deltaU_1S0_func.ev(k_grid, q_K_grid) * \
                              deltaU_1S0_dag_func.ev(q_K_grid, k_grid)
        deltaU_squared_3S1_3S1 = deltaU_3S1_3S1_func.ev(k_grid, q_K_grid) * \
                                  deltaU_3S1_3S1_dag_func.ev(q_K_grid, k_grid)
        deltaU_squared_3S1_3D1 = deltaU_3S1_3D1_func.ev(k_grid, q_K_grid) * \
                                  deltaU_3D1_3S1_dag_func.ev(q_K_grid, k_grid)
                                 
        # Build x-dependent part and integrate with resized K
        
        # These are (ntot, Ktot, xtot) arrays of
        # \theta( k_F - \abs( 1/2*K_vec +/- k_vec ) ) for every COM momentum
        # K, relative momentum k, and angle x where sign specifies the sign of
        # k_vec
        theta_kF1_K_plus_k_x = self.theta_K_k(kF_1, sign=1)
        theta_kF1_K_minus_k_x = self.theta_K_k(kF_1, sign=-1)
        theta_kF2_K_plus_k_x = self.theta_K_k(kF_2, sign=1)
        theta_kF2_K_minus_k_x = self.theta_K_k(kF_2, sign=-1)
            
        # Do integration over x first (angle-averaging) where the integration
        # weights of x are already built-in
        # This sum collapses theta_K_k_x to a matrix dependent only on K and 
        # k: (ntot, Ktot, xtot) -> (ntot, K_cutoff_index)
            
        # \int dx \theta( kF_1 - \abs( 1/2*K_vec + k_vec ) ) \times
        # \theta( kF_1 - \abs( 1/2*K_vec - k_vec ) ) / 2 
        theta_kF1_kF1_K_k = ( np.sum( 
                              theta_kF1_K_plus_k_x * theta_kF1_K_minus_k_x,
                              axis=-1 ) )[:, :K_cutoff_index] / 2
        # \int dx \theta( kF_1 - \abs( 1/2*K_vec + k_vec ) ) \times
        # \theta( kF_2 - \abs( 1/2*K_vec - k_vec ) ) / 2 
        theta_kF1_kF2_K_k = ( np.sum( 
                              theta_kF1_K_plus_k_x * theta_kF2_K_minus_k_x,
                              axis=-1 ) )[:, :K_cutoff_index] / 2
        # \int dx \theta( kF_2 - \abs( 1/2*K_vec + k_vec ) ) \times
        # \theta( kF_1 - \abs( 1/2*K_vec - k_vec ) ) / 2 
        theta_kF2_kF1_K_k = ( np.sum( 
                              theta_kF2_K_plus_k_x * theta_kF1_K_minus_k_x,
                              axis=-1 ) )[:, :K_cutoff_index] / 2
        
        # Build K integrand
        integration_K_measure = ( K_weights * K_array**2 )[:K_cutoff_index]
        integrand_k_K = integration_K_measure * (\
                          deltaU_squared_1S0 * theta_kF1_kF1_K_k + \
                          deltaU_squared_1S0/4 * (\
                            theta_kF1_kF2_K_k + theta_kF2_kF1_K_k ) + \
                          3/4*deltaU_squared_3S1_3S1 * (\
                            theta_kF1_kF2_K_k + theta_kF2_kF1_K_k ) + \
                          3/4*deltaU_squared_3S1_3D1 * (\
                            theta_kF1_kF2_K_k + theta_kF2_kF1_K_k ) )
        
        # # Build dk dK dx integrand
        # integrand_k_K_x = deltaU_squared_1S0 * theta_kF1_K_plus_k_x + \
        #                   deltaU_squared_1S0/4 * ( theta_kF1_K_plus_k_x * \
        #                   theta_kF2_K_minus_k_x + theta_kF2_K_plus_k_x * \
        #                   theta_kF1_K_minus_k_x ) + 3/4 * \
        #                   deltaU_squared_3S1_3S1 * ( theta_kF1_K_plus_k_x * \
        #                   theta_kF2_K_minus_k_x + theta_kF2_K_plus_k_x * \
        #                   theta_kF1_K_minus_k_x ) + 3/4 * \
        #                   deltaU_squared_3S1_3D1 * ( theta_kF1_K_plus_k_x * \
        #                   theta_kF2_K_minus_k_x + theta_kF2_K_plus_k_x * \
        #                   theta_kF1_K_minus_k_x )
                              
        # # Do integration over x first (angle-averaging) where the integration
        # # weights of x are already built-in
        # # This sum collapses theta_K_k_x to a matrix dependent only on K and 
        # # k: (ntot, Ktot, xtot) -> (ntot, K_cutoff_index)
        # integration_K_measure = ( K_weights * K_array**2 )[:K_cutoff_index]
        # integrand_k_K = integration_K_measure * np.sum( integrand_k_K_x,
        #                                                 axis=-1 )                             

            
        # Do integration over K
        integrand_k = np.sum( integrand_k_K, axis=-1 ) * integration_k_measure
        
        # Lastly do integration over k
        fourth_term = 1/2 * 2/np.pi * np.sum(integrand_k)

        # return first_term
        # return middle_terms
        # return fourth_term
        return first_term + middle_terms + fourth_term


if __name__ == '__main__':

    
    # --- Set up --- #
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120
    # kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    pmd = single_nucleon_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
    
    # Details of example nuclei (format is [nuclei, Z, N])
    # nuclei_list = ['C12', 6, 6]
    # nuclei_list = ['O16', 8, 8]
    # nuclei_list = ['Ca40', 20, 20]
    nuclei_list = ['Ca48', 20, 28]
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
    
    # n(q) p single-particle
    t0 = time.time()
    n_p_array = lda.local_density_approximation( q_array, pmd.n_lambda_N, 'p' )
    t1 = time.time() # End time
    mins = round( (t1 - t0) / 60.0, 2)
    print('n_lambda^p(q) done after %.2f minutes' % mins)
    n_n_array = lda.local_density_approximation( q_array, pmd.n_lambda_N, 'n' )
    
    
    # This is 4 \pi for d3p, 2 for spin sum, and 1/(2*\pi)^3 for converting
    # from sums to integrals
    overall_factor = 4*np.pi * 1/(2*np.pi)**3
    
    # BUG: this doesn't work for asymmetric nuclei for some reason
    p_a_p = q_array**2 * n_p_array / A * overall_factor
    p_a_n = q_array**2 * n_n_array / A * overall_factor

    p_a = p_a_p + p_a_n

    print( np.sum(p_a*q_weights) )
    
    
    # --- Plot --- #
    
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )
    
    # n(q) p single-particle
    plt.semilogy(q_array, p_a_p, color='red', linestyle='dashdot', 
                  label=nucleus+' (p)')
    plt.semilogy(q_array, p_a_n, color='blue', linestyle='dashed', 
                  label=nucleus+' (n)')
    plt.semilogy(q_array, p_a, color='black', label=nucleus+' (total)')
    plt.ylabel(r'$P^A(q)$')
    ax.set_xlim( [min(q_array), 4] )
    ax.set_ylim( [1e-4, 2e0])
    
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    lambda_label_location = 'center right'
    anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
                                  frameon=False)
    ax.add_artist(anchored_text)
    
    plt.show()