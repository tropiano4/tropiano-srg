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
#   Calculate n_{\lambda}(q, Q=0) for different nuclei using LDA. Split into
#   pp, pn, nn, np contributions as in Ryckebusch_2019oya.


from os import getcwd, chdir
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
from lda import load_density, LDA
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


#### TO DO ####
# 1. Just put the n_lambda_1S0, ... functions here
# 2. 1S0 (pp and nn) 1S0 (pn) 3S1-3S1 (pn) and 3S1-3D1 (pn) 1S0 (nn) 
#    contributions next


class pair_momentum_distributions(object):
    
    
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
            
            # U_matrix = U_matrix_unitless / row / col
            I_matrix = I_matrix_unitless / row / col
            delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # # BUG CHECK (unitless matrices)
            # I_matrix = np.eye( len(factor_array), len(factor_array) )
            # delta_U_matrix = U_matrix_unitless - I_matrix
            
            # Store both matrices in the dictionary
            d[channel]['I'] = I_matrix # [fm^3]
            d[channel]['delta_U'] = delta_U_matrix # [fm^3]
            
        # Save dictionary
        self.d = d
        
        
    def n_lambda_1S0(self, q, k_F):
        """
        Description.

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        k_F : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['1S0']['k_array']
        k_weights = self.d['1S0']['k_weights']
        ntot = self.d['1S0']['ntot']
        I_matrix = self.d['1S0']['I'] # [fm^3]
        delta_U_matrix = self.d['1S0']['delta_U'] # [fm^3]
        
        # Factors in front of everything
        overall_factor = 0.5 * (2/np.pi) * (4*np.pi/3*k_F**3) # [fm^-3]
        # overall_factor = 0.5 * (2/np.pi)
        
        # Spin and isospin factors appear in even powers. We can calculate
        # the squared value here.
        
        # T = 1
        # M_T = 1 or M_T = -1 (doesn't matter!)
        # Make function that takes T, M_T as argument?
        # Clebsch-Gordan coefficients for isospin are 1 here
        # < 1/2 1/2 | 1 1 > = 1 or < -1/2 -1/2 | 1 -1 > = 1
        isospin_factor_squared = 1
        
        # Select the allowed value of S and M_S
        # S = 0
        # M_S = 0
        # Make function that takes S, M_S as argument?
        # Each time we get factors of 1/\sqrt(2) or -1/\sqrt(2) but
        # appearing in even powers - thus, it's always a factor of 1/2.
        # But we sum over \sigma and \sigma' so the factor of two cancels!
        spin_factor_squared = 1/2
        
        # Select the associated partial wave channel
        # channel = '1S0'
        # This means L=0, M_L=0, J=0, M_J=0
        # < M_L M_S | J M_J> factors are always 1 in this case
        total_J_factor = 1
        # Is this ever not 1?
        
        # Find index of q in k_array
        q_index = find_q_index(q, k_array)
        
        # Bare operator: \pi / (2 w_q q^2) [fm^3]
        bare_operator = np.pi / (2*k_weights[q_index]*k_array[q_index]**2)
        # Bug checking
        # bare_operator = 1
        
        # First terms where q must be low-momentum
        if q <= 2*k_F:
            
            # matrix_element = ( I_matrix[q_index, q_index] + \
            #                     delta_U_matrix[q_index, q_index] + \
            #                     delta_U_matrix.T[q_index, q_index] ) * \
            #                   bare_operator # [fm^6]
            matrix_element = ( I_matrix[q_index, q_index] + \
                                delta_U_matrix[q_index, q_index] + \
                                delta_U_matrix.T[q_index, q_index] ) / \
                              bare_operator # [fm^6]
                             
            # Define x variable for d3Q integration
            x = q / (2*k_F)
            Q_integration_factor_low = 1 - 3/2*x + x**3/2
            
            n_lambda_low_q = overall_factor * isospin_factor_squared * \
                              (2*spin_factor_squared) * total_J_factor * \
                              matrix_element * Q_integration_factor_low # [fm^3]
            
            # n_lambda_low_q = 0          
        # q > 2*k_F
        else:
            
            n_lambda_low_q = 0

        # Now do the last term with two \delta_U's (should be an array over
        # k)
        # matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot] * \
        #                         bare_operator # [fm^9]
        
        matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
                                delta_U_matrix.T[q_index, :ntot] / \
                                bare_operator # [fm^9]
                                
        # Define y variable for d3Q integration (function of k)
        y_array = k_array / (2*k_F)
        Q_integration_factor_high = 1 - 3/2*y_array + y_array**3/2

        # n_\lambda_high_q before k-integration [fm^6]
        n_lambda_k_vector = 0.5*overall_factor * isospin_factor_squared**2 * \
                            (4*spin_factor_squared**2) * total_J_factor**2 * \
                            matrix_element_vector * Q_integration_factor_high
        
        # Do the k-integration
        integration_measure = 2/np.pi * k_weights * k_array**2 # [fm^-3]
        # Bug test
        # integration_measure = 1
        integrand = integration_measure * n_lambda_k_vector # [fm^3]
            
        # Set 2*k_F cutoff in integration limit where we prevent overshooting
        # 2*k_F (meaning \Lambda will never be more than 2*k_F)
        # DO THIS LATER
        cutoff_index = find_q_index(2*k_F, k_array)
        # while k_array[cutoff_index] > 2*k_F:
        #     cutoff_index -= 1
            
        integrand_resized = integrand[:cutoff_index+1]
        
        # Integrate to 2*k_F
        n_lambda_high_q = np.sum(integrand_resized) # [fm^3]
        # n_lambda_high_q = 0
        
        # Combine for full contribution
        n_lambda = n_lambda_low_q + n_lambda_high_q

        return n_lambda
    

# ### Here's what's in your notebook from before ###
# # Calculate SRG transformation with eigenvectors of initial and evolved Hamiltonians dividing out 
# # integration factors
# U_matrix_1s0 = SRG_unitary_transformation(H_initial_1s0, H_evolved_1s0) / row / col
# U_matrix_3s1 = SRG_unitary_transformation(H_initial_3s1, H_evolved_3s1) / row_cc / col_cc
    
# # Index of k_0 in k_array
# k_0_index = op.find_q_index(k_0, k_array)
        
# # Calculate |U(k_0, q)_{3S1}|^2 [fm^3] (divide 2/\pi k_i k_j \sqrt(w_i w_j))
# # after calculation
# numerator_array = U_matrix_3s1[k_0_index, :ntot] * U_matrix_3s1.T[:ntot, k_0_index] + \
#                   U_matrix_3s1[k_0_index, ntot:] * U_matrix_3s1.T[ntot:, k_0_index]

# NO BARE OPERATOR!
    
    
    def n_lambda_3S1_3S1(self, q, k_F):
        """
        Description.

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        k_F : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # HOW DOES THIS DEPEND ON \alpha(r)?

        return None
    
    def delta_U(self):
        
        return self.d['1S0']['delta_U']
    
    def I(self):
        
        return self.d['1S0']['I']


if __name__ == '__main__':

    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax=10.0, kmid=2.0, ntot=120)
    # q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax=30.0, kmid=4.0, ntot=120)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    pmd = pair_momentum_distributions(kvnn, lamb)
    
    # Details of example nuclei (format is [nuclei, Z, N])
    nuclei_list = ['O16', 8, 8]
    # nuclei_list = ['Ca48', 20, 28]
    nucleus = nuclei_list[0]
    Z = nuclei_list[1]
    N = nuclei_list[2]
    r_array, rho_array_p = load_density(nucleus, 'proton', Z, N)
    r_array, rho_array_n = load_density(nucleus, 'neutron', Z, N)
    
    k_F_array = np.linspace(0.001, 1.8, 100)
    
    # --- Write the relevant tables --- #
    # Comment this out when you're done
    # Load LDA class and input pair momentum distribution function
    # lda = LDA(pmd.n_lambda_1S0, 'n_lambda_1S0', kvnn, r_array, rho_array, 
    #           True, q_array, k_F_array)
    
    # No interpolation here
    lda_p = LDA(pmd.n_lambda_1S0, 'n_lambda_1S0', kvnn, r_array, rho_array_p)
    lda_n = LDA(pmd.n_lambda_1S0, 'n_lambda_1S0', kvnn, r_array, rho_array_n)
    
    
    # --- Plot n_lambda_pair_exp for pp pairs --- #
    
    # q_array_fine = np.linspace(0.1, 6.0, 100)
    q_array_fine = q_array
    # n_lambda_pair_exp_array = lda.local_density_approximation(q_array_fine)
    # BUG TESTING
    n_lambda_pp_exp_array = lda_p.local_density_approximation(q_array_fine, 
                                                              pmd.n_lambda_1S0)
    n_lambda_nn_exp_array = lda_n.local_density_approximation(q_array_fine, 
                                                              pmd.n_lambda_1S0)
    
    # # BUG TEST
    # n_lambda_pp_exp_array = lda_p.local_density_approximation(q_array_fine, 
    #                         pmd.n_lambda_1S0) / factor_array
    # n_lambda_nn_exp_array = lda_n.local_density_approximation(q_array_fine, 
    #                         pmd.n_lambda_1S0) / factor_array
    
    plt.clf()
    
    # plt.plot(q_array_fine, n_lambda_pair_exp_array, label=nucleus)
    
    plt.semilogy(q_array_fine, n_lambda_pp_exp_array, color='blue', 
                  label=nucleus+' (pp)')  
    plt.semilogy(q_array_fine, n_lambda_nn_exp_array, color='red', 
                  linestyle='dashed', label=nucleus+' (nn)')  
    plt.ylabel(r'$<n_{\lambda}^{NN}(q)>$' + ' [fm' + r'$^3$' + ']')
    
    # plt.semilogy(q_array_fine, n_lambda_pp_exp_array*factor_array, color='blue', 
    #               label=nucleus+' (pp)')
    # plt.semilogy(q_array_fine, n_lambda_nn_exp_array*factor_array, color='red', 
    #               linestyle='dashed', label=nucleus+' (nn)')
    # plt.ylabel(r'$<n_{\lambda}^{NN}(q)>q^2 dq$')
    
    # plt.semilogy(q_array_fine, n_lambda_pp_exp_array*factor_array**2, color='blue', 
    #               label=nucleus+' (pp)')  
    # plt.semilogy(q_array_fine, n_lambda_nn_exp_array*factor_array**2, color='red', 
    #               linestyle='dashed', label=nucleus+' (nn)')  
    # plt.ylabel(r'$<n_{\lambda}^{NN}(q)>(q^2 dq)^2$' + ' [fm' + r'$^{-3}$' + ']')
        
    # BUG TESTING 
    plt.xlim( [min(q_array_fine), 6] )
    plt.ylim( [1e-6, 1e1])
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('q [fm' + r'$^{-1}$' + ']')
    
    plt.title('kvnn = %d, ' % kvnn + r'$\lambda=%.2f$' % lamb + 
              ' [fm' + r'$^{-1}$' + ']')
    plt.show()