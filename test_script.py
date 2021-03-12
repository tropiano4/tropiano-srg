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
from matplotlib.offsetbox import AnchoredText
import numpy as np
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
import Figures.figures_functions as ff
from lda import load_density, LDA
import observables as ob
from operators import find_q_index, momentum_projection_operator
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
            
            # Store both matrices in the dictionary
            d[channel]['I'] = I_matrix_unitless
            d[channel]['delta_U'] = delta_U_matrix_unitless
            
            # U_matrix = U_matrix_unitless / row / col
            # I_matrix = I_matrix_unitless / row / col
            # delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # # Store both matrices in the dictionary
            # d[channel]['I'] = I_matrix # [fm^3]
            # d[channel]['delta_U'] = delta_U_matrix # [fm^3]
            
        # Save dictionary
        self.d = d
        
        
    def n_lambda_1S0(self, q, k_F, pair='pp'):
        """
        Description.

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        k_F : TYPE
            DESCRIPTION.
        pair : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['1S0']['k_array']
        k_weights = self.d['1S0']['k_weights']
        ntot = self.d['1S0']['ntot']
        I_matrix = self.d['1S0']['I']
        delta_U_matrix = self.d['1S0']['delta_U']
        
        # Create bare operator here ??
        # n(q, Q) = \pi / (2 w_q q^2) *  1 / (4 \pi w_Q Q^2) [fm^6]
        
        # Find index of q in k_array
        deltas_array = np.zeros(ntot)
        q_index = find_q_index(q, k_array)
        deltas_array[q_index] = 1 
        # factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array * \
        #                 np.sqrt(4*np.pi*k_weights) * k_array # [fm^-3]
        factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array # [fm^-3/2]
        # factor_array = 1
        deltas_array /= factor_array
        operator_row, operator_col = np.meshgrid(deltas_array, deltas_array)
        operator = operator_row * operator_col # [fm^6]
        
        # Factors in front of everything
        # overall_factor = 0.5 * (2/np.pi) * (4*np.pi/3*k_F**3) # [fm^-3]
        # Techinically you should get rid of this factor of 2 for the true
        # pp or nn pair distribution because you haven't summed over \tau and
        # \tau' in this code
        overall_factor = (2/np.pi) * (4*np.pi/3*k_F**3) # [fm^-3]
        
        # Spin and isospin factors appear in even powers. We can calculate
        # the squared value here.
        
        # T = 1
        if pair == 'pn' or pair == 'np':
            # M_T = 0
            # Make function that takes T, M_T as argument?
            # Clebsch-Gordan coefficients for isospin are 1/\sqrt(2) here
            # < 1/2 -1/2 | 1 0 > = 1/\sqrt(2) or < -1/2 1/2 | 1 0 > = 1/\sqrt(2)
            isospin_factor_squared = 1/2
        
        # pp or nn
        else:
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
        
        # First terms where q must be low-momentum
        #if q <= 2*k_F:
        if q <= k_F:
            
            # matrix_element = ( I_matrix[q_index, q_index] + \
            #                     delta_U_matrix[q_index, q_index] + \
            #                     delta_U_matrix.T[q_index, q_index] ) * \
            #                   bare_operator
            
            # matrix_element = ( I_matrix[q_index, q_index] + \
            #                     delta_U_matrix[q_index, q_index] + \
            #                     delta_U_matrix.T[q_index, q_index] ) / \
            #                   bare_operator
            
            # matrix_element = I_matrix @ operator @ I_matrix + \
            #                  delta_U_matrix @ operator @ I_matrix + \
            #                  I_matrix @ operator @ delta_U_matrix.T
                             
            matrix_element = ( I_matrix[q_index, q_index] + \
                               delta_U_matrix[q_index, q_index] + \
                               delta_U_matrix.T[q_index, q_index] ) 
                             
            # Define x variable for d3Q integration
            # x = q / (2*k_F)
            # Q_integration_factor_low = 1 - 3/2*x + x**3/2
            x = q / k_F
            Q_integration_factor_low = 8*(1 - 3/2*x + x**3/2)
            
            # n_lambda_low_q = overall_factor * isospin_factor_squared * \
            #                  (2*spin_factor_squared) * total_J_factor * \
            #                  matrix_element[q_index, q_index] * \
            #                  Q_integration_factor_low
            # 02/16/21: currently [fm^-3] tail shows mesh-dependence
            # 02/17/21: trying operator [fm^6] so this is [fm^3]
            n_lambda_low_q = overall_factor * isospin_factor_squared * \
                              (2*spin_factor_squared) * total_J_factor * \
                              matrix_element * Q_integration_factor_low
            
            # n_lambda_low_q = 0          
        # q > 2*k_F
        else:
            
            n_lambda_low_q = 0

        # Now do the last term with two \delta_U's (should be an array over
        # k)
        
        # matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot] * \
        #                         bare_operator
        # matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot] / \
        #                         bare_operator
        
        matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
                                operator[q_index, q_index] * \
                                delta_U_matrix.T[q_index, :ntot]
                                # 02/16/21: currently [unitless]
                                # 02/17/21: trying fm^6
                                
        # Define y variable for d3Q integration (function of k)
        # y_array = k_array / (2*k_F)
        # Q_integration_factor_high = 1 - 3/2*y_array + y_array**3/2
        y_array = k_array / k_F
        Q_integration_factor_high = 8*(1 - 3/2*y_array + y_array**3/2)

        # n_\lambda_high_q before k-integration
        n_lambda_k_vector = 0.5*overall_factor * isospin_factor_squared**2 * \
                            (4*spin_factor_squared**2) * total_J_factor**2 * \
                            matrix_element_vector * Q_integration_factor_high
                            # 02/16/21: currently [fm^-3]
                            # 02/17/21: trying fm^3
        
        # Do the k-integration
        integration_measure = 2/np.pi * k_weights * k_array**2 # [fm^-3]
        # Bug test
        # integration_measure = 1
        integrand = integration_measure * n_lambda_k_vector # [fm^3]
        # 02/16/21: currently [fm^-3]
        # 02/17/21: currently fm^3
            
        # Set 2*k_F cutoff in integration limit where we prevent overshooting
        # 2*k_F (meaning \Lambda will never be more than 2*k_F)
        # DO THIS LATER
        cutoff_index = find_q_index(k_F, k_array)
        # while k_array[cutoff_index] > 2*k_F:
        #     cutoff_index -= 1
            
        integrand_resized = integrand[:cutoff_index+1]
        
        # Integrate to 2*k_F
        n_lambda_high_q = np.sum(integrand_resized) # [fm^3]
        # n_lambda_high_q = 0
        
        # Combine for full contribution
        n_lambda = n_lambda_low_q + n_lambda_high_q
        # 02/16/21: currently [fm^-3]

        return n_lambda


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
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['3S1']['k_array']
        k_weights = self.d['3S1']['k_weights']
        ntot = self.d['3S1']['ntot']
        I_matrix = self.d['3S1']['I']
        delta_U_matrix = self.d['3S1']['delta_U']
        
        # Create bare operator here
        # n(q, Q) = \pi / (2 w_q q^2) *  1 / (4 \pi w_Q Q^2) [fm^6]
        
        # Find index of q in k_array
        deltas_array = np.zeros(ntot)
        q_index = find_q_index(q, k_array)
        deltas_array[q_index] = 1 
        # factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array * \
        #                np.sqrt(4*np.pi*k_weights) * k_array # [fm^-3]
        factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array # [fm^-3/2]
        # factor_array = 1
        deltas_array /= factor_array
        operator_row, operator_col = np.meshgrid(deltas_array, deltas_array)
        operator = operator_row * operator_col # [fm^6]
        
        # Matrix of zeros (m x m) for coupled-channel operator
        o = np.zeros( (ntot, ntot) )
    
        # Build coupled channel operator
        operator = np.vstack( ( np.hstack( (operator, o) ),
                                np.hstack( (o, operator) ) ) )
        
        # Factors in front of everything
        # overall_factor = 0.5 * (2/np.pi) * (4*np.pi/3*k_F**3) # [fm^-3]
        # Techinically you should get rid of this factor of 2 for the true
        # pp or nn pair distribution because you haven't summed over \tau and
        # \tau' in this code
        overall_factor = (2/np.pi) * (4*np.pi/3*k_F**3) # [fm^-3]
        
        # Spin and isospin factors appear in even powers. We can calculate
        # the squared value here.
        
        # T = 0
        # M_T = 0
        # Make function that takes T, M_T as argument?
        # Clebsch-Gordan coefficients for isospin are 1/\sqrt(2) here
        # < 1/2 -1/2 | 1 0 > = 1/\sqrt(2) or < -1/2 1/2 | 1 0 > = 1/\sqrt(2)
        # isospin_factor_squared = 1/2
        
        # Select the allowed value of S and M_S
        # S = 1
        # M_S = 0
        # Make function that takes S, M_S as argument?
        # Each time we get factors of 1/\sqrt(2) or -1/\sqrt(2) but
        # appearing in even powers - thus, it's always a factor of 1/2.
        # But we sum over \sigma and \sigma' so the factor of two cancels!
        # spin_factor_squared = 1/2
        
        # Select the associated partial wave channel
        # channel = '1S0'
        # This means L=0, M_L=0, J=0, M_J=0
        # < M_L M_S | J M_J> factors are always 1 in this case
        # total_J_factor = 1
        # Is this ever not 1?
        
        # First terms where q must be low-momentum
        #if q <= 2*k_F:
        if q <= k_F:
            
            # matrix_element = ( I_matrix[q_index, q_index] + \
            #                     delta_U_matrix[q_index, q_index] + \
            #                     delta_U_matrix.T[q_index, q_index] ) * \
            #                   bare_operator # [fm^6]
            
            # matrix_element = ( I_matrix[q_index, q_index] + \
            #                     delta_U_matrix[q_index, q_index] + \
            #                     delta_U_matrix.T[q_index, q_index] ) / \
            #                   bare_operator # [fm^6]
            
            # matrix_element = I_matrix @ operator @ I_matrix + \
            #                  delta_U_matrix @ operator @ I_matrix + \
            #                  I_matrix @ operator @ delta_U_matrix.T
            
            matrix_element = I_matrix[q_index, q_index] + \
                             delta_U_matrix[q_index, q_index] + \
                             delta_U_matrix.T[q_index, q_index]
                             
            # Define x variable for d3Q integration
            # x = q / (2*k_F)
            # Q_integration_factor_low = 1 - 3/2*x + x**3/2
            x = q / k_F
            Q_integration_factor_low = 8*(1 - 3/2*x + x**3/2)
            
            # n_lambda_low_q = overall_factor * \
            #                  matrix_element[q_index, q_index] * \
            #                  Q_integration_factor_low
            # 02/16/21: currently [fm^-3]
            n_lambda_low_q = overall_factor * matrix_element * \
                             Q_integration_factor_low
            
            # n_lambda_low_q = 0          
        # q > 2*k_F
        else:
            
            n_lambda_low_q = 0

        # Now do the last term with two \delta_U's (should be an array over
        # k)
        
        # matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot] * \
        #                         bare_operator
        # matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
        #                         delta_U_matrix.T[q_index, :ntot] / \
        #                         bare_operator
        
        matrix_element_vector = delta_U_matrix[:ntot, q_index] * \
                                operator[q_index, q_index] * \
                                delta_U_matrix.T[q_index, :ntot] + \
                                delta_U_matrix[ntot:, ntot+q_index] * \
                                operator[ntot+q_index, ntot+q_index] * \
                                delta_U_matrix.T[ntot+q_index, ntot:]
                                # 02/16/21: currently [unitless]
                                
        # Define y variable for d3Q integration (function of k)
        # y_array = k_array / (2*k_F)
        # Q_integration_factor_high = 1 - 3/2*y_array + y_array**3/2
        y_array = k_array / k_F
        Q_integration_factor_high = 8*(1 - 3/2*y_array + y_array**3/2)

        # n_\lambda_high_q before k-integration
        n_lambda_k_vector = 0.5*overall_factor * \
                            matrix_element_vector * Q_integration_factor_high
                            # 02/16/21: currently [fm^-3]
        
        # Do the k-integration
        integration_measure = 2/np.pi * k_weights * k_array**2 # [fm^-3]
        # Bug test
        # integration_measure = 1
        integrand = integration_measure * n_lambda_k_vector
        # 02/16/21: currently [fm^-3]
            
        # Set 2*k_F cutoff in integration limit where we prevent overshooting
        # 2*k_F (meaning \Lambda will never be more than 2*k_F)
        # DO THIS LATER
        cutoff_index = find_q_index(k_F, k_array)
        # while k_array[cutoff_index] > 2*k_F:
        #     cutoff_index -= 1
            
        integrand_resized = integrand[:cutoff_index+1]
        
        # Integrate to 2*k_F
        n_lambda_high_q = np.sum(integrand_resized) # [fm^3]
        # n_lambda_high_q = 0
        
        # Combine for full contribution
        n_lambda = n_lambda_low_q + n_lambda_high_q
        # 02/16/21: currently [fm^-3]

        return n_lambda
    
    
    def n_infty(self, q, k_F):
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['1S0']['k_array']
        k_weights = self.d['1S0']['k_weights']
        ntot = self.d['1S0']['ntot']
        
        # Create bare operator here
        # n(q, Q) = \pi / (2 w_q q^2) *  1 / (4 \pi w_Q Q^2) [fm^6]
        
        # Find index of q in k_array
        deltas_array = np.zeros(ntot)
        q_index = find_q_index(q, k_array)
        deltas_array[q_index] = 1 
        # factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array * \
        #                np.sqrt(4*np.pi*k_weights) * k_array # [fm^-3]
        factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array # [fm^-3/2]
        deltas_array /= factor_array
        operator_row, operator_col = np.meshgrid(deltas_array, deltas_array)
        # operator = operator_row * operator_col
        # operator_qq = operator[q_index, q_index]
        operator_qq = 1
        
    
        if q <= k_F:
            
            x = q/k_F
            
            overall_factor = 8*4*np.pi/3*k_F**3 # fm^-3
            # overall_factor = 1
            x_part = 1-3/2*x+x**3/2 # unitless
            
            return overall_factor * x_part * operator_qq
            # 02/16/21: currently [fm^-3]
        
        else:
            
            return 0
        
    def n_single_particle(self, p, k_F):
        
        # This should just be \theta(k_F-p)
        
        if p <= k_F:
        # if True:
            
            return 1
        
        else:
            
            return 0
        
    def n_deuteron(self, q, psi_vector_unitless):
        
        # Load momentum mesh, I, and delta_U from dictionary
        k_array = self.d['3S1']['k_array']
        k_weights = self.d['3S1']['k_weights']
        ntot = self.d['3S1']['ntot']
        delta_U_matrix_unitless = self.d['3S1']['delta_U']
        
        factor_array = np.sqrt(2/np.pi*k_weights) * k_array # fm^-3/2
        
        ### NEW ###
        
        factor_array_long = np.concatenate( (factor_array, factor_array) )
        row, col = np.meshgrid(factor_array_long, factor_array_long)
        
        
        # # USE MOMENTUM PROJECTION OPERATOR HERE
        # # units are fm^3
        # bare_operator = momentum_projection_operator(q, k_array, k_weights, 
        #                                              '3S1', smeared=False)
        
        # I_matrix_unitless = self.d['3S1']['I']
        # psi_vector = psi_vector_unitless
        # first_term = psi_vector.T @ I_matrix_unitless @ bare_operator @ \
        #              I_matrix_unitless @ psi_vector
        # second_term = psi_vector.T @ delta_U_matrix_unitless @ \
        #               bare_operator @ I_matrix_unitless @ psi_vector
        # third_term = second_term
        # fourth_term = psi_vector.T @ delta_U_matrix_unitless @ \
        #               bare_operator @ delta_U_matrix_unitless.T @ \
        #               psi_vector
                     
        # Make delta_U and psi have units
        # fm^3 and dimensions 2*ntot, 2*ntot
        delta_U_matrix = delta_U_matrix_unitless / row / col
        # Convert wave function to fm^3/2
        psi_vector = psi_vector_unitless / factor_array_long
        
        # Find index of q
        q_index = find_q_index(q, k_array)
        
        # First term
        # first_term = psi_vector[q_index]**2 # fm^3
        # 3S1^2 + 3D1^2
        first_term = psi_vector[q_index]**2 + psi_vector[ntot+q_index]**2
        
        # This is incorrect at the moment (@ operator not doing what you
        # think it's doing)
        # Second term
        # integration_measure = k_array**2 * k_weights # fm^-3
        # integrand_k = psi_vector.T @ delta_U_matrix[:ntot, q_index] * \
        #               psi_vector[q_index] # fm^6
        # second_term = 2/np.pi * np.sum( integration_measure * integrand_k )
        
        # Third term
        # integrand_kp = psi_vector.T[q_index] * \
        #                delta_U_matrix.T[q_index, :ntot] @ psi_vector
        # third_term = 2/np.pi * np.sum( integration_measure * integrand_kp )
        
        # Fourth term
        # fourth_term = (2/np.pi**2) * ( ( integration_measure * psi_vector.T ) @\
        #               ( delta_U_matrix[:ntot, q_index] * \
        #                 delta_U_matrix.T[q_index, :ntot] + \
        #                 delta_U_matrix[:ntot, ntot+q_index] * \
        #                 delta_U_matrix.T[ntot+q_index, :ntot] ) @ \
        #               ( integration_measure * psi_vector ) )
        
        # DO A DOUBLE FOUR LOOP TO MAKE SURE YOU'RE DOING THIS CORRECTLY
        second_term = 0
        third_term = 0
        fourth_term = 0
        for i, k in enumerate(k_array):
            
            # 3S1 3S1-3S1 3S1
            # 3S1 3S1-3D1 3D1
            # 3D1 3D1-3S1 3S1
            # 3D1 3D1-3D1 3D1
            second_term += 2/np.pi * k**2 * k_weights[i] * ( \
                            psi_vector.T[i] * delta_U_matrix[i, q_index] * \
                            psi_vector[q_index] + \
                                
                            psi_vector.T[i] * \
                            delta_U_matrix[i, ntot+q_index] * \
                            psi_vector[ntot+q_index] + \
                            
                            psi_vector.T[ntot+i] * \
                            delta_U_matrix[ntot+i, q_index] * \
                            psi_vector[q_index] + \
                            
                            psi_vector.T[ntot+i] * \
                            delta_U_matrix[ntot+i, ntot+q_index] * \
                            psi_vector[ntot+q_index] )

            for j, kp in enumerate(k_array):
                
                # \psi.T \delta_U \delta_U.T \psi
                # 3S1 3S1-3S1 3S1-3S1 3S1
                # 3S1 3S1-3D1 3D1-3S1 3S1
                
                # 3S1 3S1-3S1 3S1-3D1 3D1
                # 3S1 3S1-3D1 3D1-3D1 3D1
                
                # 3D1 3D1-3S1 3S1-3S1 3S1
                # 3D1 3D1-3D1 3D1-3S1 3S1
                
                # 3D1 3D1-3S1 3S1-3D1 3D1
                # 3D1 3D1-3D1 3D1-3D1 3D1
                fourth_term +=  (2/np.pi)**2 * k**2 * k_weights[i] *\
                                kp**2 * k_weights[j] * (\
                                psi_vector.T[i] * ( 
                                    delta_U_matrix[i, q_index] * \
                                    delta_U_matrix.T[q_index, j] +
                                    delta_U_matrix[i, ntot+q_index] * \
                                    delta_U_matrix.T[ntot+q_index, j] ) * \
                                psi_vector[j] + \
                                    
                                psi_vector.T[i] * ( 
                                    delta_U_matrix[i, q_index] * \
                                    delta_U_matrix.T[q_index, j+ntot] +
                                    delta_U_matrix[i, ntot+q_index] * \
                                    delta_U_matrix.T[ntot+q_index, j+ntot] ) * \
                                psi_vector[j+ntot] + \
                                    
                                psi_vector.T[i+ntot] * ( 
                                    delta_U_matrix[i+ntot, q_index] * \
                                    delta_U_matrix.T[q_index, j] +
                                    delta_U_matrix[i+ntot, ntot+q_index] * \
                                    delta_U_matrix.T[ntot+q_index, j] ) * \
                                psi_vector[j] + \
                                
                                psi_vector.T[i+ntot] * ( 
                                    delta_U_matrix[i+ntot, q_index] * \
                                    delta_U_matrix.T[q_index, j+ntot] +
                                    delta_U_matrix[i+ntot, ntot+q_index] * \
                                    delta_U_matrix.T[ntot+q_index, j+ntot] ) * \
                                psi_vector[j+ntot] )
                                    
        third_term = second_term
                                    
        # low_q_part = first_term + second_term + third_term
        # high_q_part = fourth_term
        
        # return low_q_part, high_q_part
        
        return first_term, second_term + third_term, fourth_term


if __name__ == '__main__':

    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120
    # kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    pmd = pair_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
    
    # Details of example nuclei (format is [nuclei, Z, N])
    # nuclei_list = ['O16', 8, 8]
    # nuclei_list = ['Ca40', 20, 20]
    # nuclei_list = ['Ca48', 20, 28]
    # nuclei_list = ['Pb208', 82, 126]
    # nucleus = nuclei_list[0]
    # Z = nuclei_list[1]
    # N = nuclei_list[2]
    # A = N + Z
    # r_array, rho_array_p = load_density(nucleus, 'proton', Z, N)
    # r_array, rho_array_n = load_density(nucleus, 'neutron', Z, N)
    
    # k_F_array = np.linspace(0.001, 1.8, 100)
    
    # --- Write the relevant tables --- #
    # Comment this out when you're done
    # Load LDA class and input pair momentum distribution function
    # lda = LDA(pmd.n_lambda_1S0, 'n_lambda_1S0', kvnn, r_array, rho_array, 
    #           True, q_array, k_F_array)
    # GET RID OF INTERPOLATION
    
    # No interpolation here
    # THIS WILL LIKELY CAUSE BUGS LATER: FIX!
    # lda_p = LDA(pmd.n_lambda_1S0, 'n_lambda_1S0', kvnn, r_array, rho_array_p)
    # lda_n = LDA(pmd.n_lambda_1S0, 'n_lambda_1S0', kvnn, r_array, rho_array_n)
    
    
    # --- Plot n_lambda_pair_exp for pp pairs --- #
    
    # q_array_fine = np.linspace(0.01, 6.0, 1000)
    q_array_fine = q_array

    # BUG TESTING
    # n_lambda_pp_exp_array = lda_p.local_density_approximation(q_array_fine, 
    #                                                           pmd.n_lambda_1S0)
    # n_lambda_nn_exp_array = lda_n.local_density_approximation(q_array_fine, 
    #                                                           pmd.n_lambda_1S0)
    # n_lambda_pn_exp_array = lda_p.local_density_approximation(q_array_fine, 
    #                                                           pmd.n_lambda_3S1_3S1)
    # n_infty_pp_exp_array = lda_p.local_density_approximation(q_array_fine, 
    #                                                          pmd.n_infty)
    # n_infty_p_exp_array = lda_p.local_density_approximation(q_array_fine, 
    #                                                         pmd.n_single_particle)
    # n_infty_n_exp_array = lda_n.local_density_approximation(q_array_fine, 
    #                                                         pmd.n_single_particle)
    # overall_factor = 4*np.pi * 2 * 1/(2*np.pi)**3
    # p_a_ipm_p = q_array_fine**2 * n_infty_p_exp_array / A * overall_factor
    # p_a_ipm_n = q_array_fine**2 * n_infty_n_exp_array / A * overall_factor
    
    # p_a_ipm = p_a_ipm_p + p_a_ipm_n
    
    # Split into four pp, pn, nn, np pieces
    # p_a_ipm_pp = p_a_ipm_p * (Z-1)/(A-1)
    # p_a_ipm_pn = p_a_ipm_p * N/(A-1)
    # p_a_ipm_np = p_a_ipm_n * Z/(A-1)
    # p_a_ipm_nn = p_a_ipm_n * (N-1)/(A-1)
    
    # p_a_ipm = p_a_ipm_pp + p_a_ipm_pn + p_a_ipm_np + p_a_ipm_nn
    
    
    # This is 4 \pi for d3p, 2 for spin sum, and 1/(2*\pi)^3 for converting
    # from sums to integrals
    # print(np.sum(p_a_ipm*q_weights))
    
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )
    
    # plt.plot(q_array_fine, n_lambda_pair_exp_array, label=nucleus)
    
    # plt.semilogy(q_array_fine, n_lambda_pp_exp_array, color='blue', 
    #               linestyle='dashdot', label=nucleus+' (pp)')  
    # plt.semilogy(q_array_fine, n_lambda_nn_exp_array, color='red', 
    #               linestyle='dashed', label=nucleus+' (nn)') 
    # plt.semilogy(q_array_fine, n_infty_pp_exp_array, color='green', 
    #               linestyle='dotted', label=nucleus+' (pp), unevolved')  
    # plt.semilogy(q_array_fine, n_lambda_pn_exp_array, color='black', 
    #               linestyle='solid', label=nucleus+' (pn)')  
    # # plt.ylabel(r'$<n_{\lambda}^{NN}(q)>$' + ' [fm' + r'$^3$' + ']')
    # plt.ylabel(r'$<n_{\lambda}^{NN}(q)>$' + ' [fm' + r'$^{-3}$' + ']')
    
    # plt.semilogy(q_array_fine, n_lambda_pp_exp_array*factor_array, 
    #              color='blue', linestyle='dashdot', label=nucleus+' (pp)')  
    # plt.semilogy(q_array_fine, n_lambda_nn_exp_array*factor_array,
    #              color='red', linestyle='dashed', label=nucleus+' (nn)') 
    # plt.semilogy(q_array_fine, n_infty_pp_exp_array*factor_array, 
    #              color='green', linestyle='dotted', label=nucleus+' (pp), unevolved')  
    # plt.semilogy(q_array_fine, n_lambda_pn_exp_array*factor_array,
    #              color='black', linestyle='solid', label=nucleus+' (pn)') 
    # plt.ylabel(r'$<n_{\lambda}^{NN}(q)>q^2 dq$' + ' [fm' + r'$^{-6}$' + ']')
    
    # plt.semilogy(q_array_fine, n_lambda_pp_exp_array/factor_array, 
    #               color='blue', linestyle='dashdot', label=nucleus+' (pp)')  
    # plt.semilogy(q_array_fine, n_lambda_nn_exp_array/factor_array, 
    #               color='red', linestyle='dashed', label=nucleus+' (nn)') 
    # plt.semilogy(q_array_fine, n_infty_pp_exp_array/factor_array, 
    #               color='green', linestyle='dotted', label=nucleus+' (pp), unevolved')  
    # plt.semilogy(q_array_fine, n_lambda_pn_exp_array/factor_array,
    #               color='black', linestyle='solid', label=nucleus+' (pn)') 
    # plt.ylabel(r'$<n_{\lambda}^{NN}(q)>/(q^2 dq)$' + ' [unitless]')
    
    # plt.semilogy(q_array_fine, p_a_ipm_p, color='red', linestyle='dashdot', 
    #              label=nucleus+' (p)')
    # plt.semilogy(q_array_fine, p_a_ipm_n, color='blue', linestyle='dashed', 
    #              label=nucleus+' (n)') 
    # plt.semilogy(q_array_fine, p_a_ipm_pp, color='red', linestyle='dashdot', 
    #              label=nucleus+' (pp)')
    # plt.semilogy(q_array_fine, p_a_ipm_nn, color='blue', linestyle='dashed', 
    #              label=nucleus+' (nn)')
    # plt.semilogy(q_array_fine, p_a_ipm_pn+p_a_ipm_np, color='green', linestyle='dotted', 
    #              label=nucleus+' (pn+np)') 
    # plt.semilogy(q_array_fine, p_a_ipm, color='black', label=nucleus+' (total)')
    # plt.ylabel(r'$P^A(p)$')
    
    
    # --- use deuteron to check pair momentum distribution --- #
    
    # Load evolved wave function here (unitless)
    H_initial = vnn.load_hamiltonian(kvnn, '3S1', kmax, kmid, ntot)
    H_evolved = vnn.load_hamiltonian(kvnn, '3S1', kmax, kmid, ntot, 
                                     method='srg', generator='Wegner', 
                                     lamb=lamb)
    U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
    
    psi_deuteron = ob.wave_function(H_initial, U=U_matrix_unitless)
    
    # n_d_low_array = np.zeros(ntot)
    # n_d_high_array = np.zeros(ntot)
    # for iq, q in enumerate(q_array):
    #     low_q, high_q = pmd.n_deuteron(q, psi_deuteron)
    #     n_d_low_array[iq] = low_q
    #     n_d_high_array[iq] = high_q
    # n_d_total_array = n_d_low_array + n_d_high_array
    
    n_d_1 = np.zeros(ntot)
    n_d_23 = np.zeros(ntot)
    n_d_4 = np.zeros(ntot)
    for iq, q in enumerate(q_array):
        first_term, second_third_term, fourth_term = pmd.n_deuteron(q, psi_deuteron)
        n_d_1[iq] = first_term
        n_d_23[iq] = second_third_term
        n_d_4[iq] = fourth_term
    n_d_total_array = n_d_1 + n_d_23 + n_d_4
    
    # Correct momentum distribution for comparison
    n_d_exact = ( ob.wave_function(H_initial)[:ntot]**2 + \
                  ob.wave_function(H_initial)[ntot:]**2 ) / factor_array
        
    ax.semilogy(q_array_fine, n_d_total_array, color='black',
                 linestyle='solid', label='total')
    ax.semilogy(q_array_fine, n_d_1, color='blue', linestyle='dotted', 
                 label='First term') 
    ax.semilogy(q_array_fine, abs(n_d_23), color='purple', linestyle='dotted', 
                 label='|Second+Third Term|') 
    ax.semilogy(q_array_fine, n_d_4, color='red', linestyle='dotted', 
                 label='Fourth Term') 
    ax.semilogy(q_array_fine, n_d_exact, color='green', linestyle='dashed', 
                 label='exact') 
    ax.set_ylabel(r'$<n^{\lambda}_d(q)>$' + ' [fm' + r'$^3$' + ']')
        
        
        
    # BUG TESTING 
    ax.set_xlim( [min(q_array_fine), 4] )
    ax.set_ylim( [1e-5, 1e3])
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    lambda_label_location = 'center right'
    
    anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
                                 frameon=False)
    ax.add_artist(anchored_text)

    plt.show()