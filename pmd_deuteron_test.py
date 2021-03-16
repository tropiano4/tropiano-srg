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
#   Calculate n_{\lambda}(q) for deuteron starting from second quantization
#   derivation and using expansion of U(k, k'):
#   U(k, k') = 1 + \sum_{1234} < 12 | \delta U | 34 > a^t_1 a^t_2 a_4 a_3


import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
# Scripts made by A.T.
import observables as ob
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


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
            
            # Store in dictionary
            d[channel]['delta_U'] = delta_U_matrix_unitless
            
        # Save dictionary
        self.d = d

        
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

        return first_term, second_term + third_term, fourth_term


if __name__ == '__main__':

    
    # --- Set up --- #
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    # kmax, kmid, ntot = 10.0, 2.0, 120
    kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    pmd = pair_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
    
    
    # --- Use deuteron to check pair momentum distribution --- #
    
    # Load evolved wave function here (unitless)
    H_initial = vnn.load_hamiltonian(kvnn, '3S1', kmax, kmid, ntot)
    H_evolved = vnn.load_hamiltonian(kvnn, '3S1', kmax, kmid, ntot, 
                                     method='srg', generator='Wegner', 
                                     lamb=lamb)
    U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)
    
    psi_deuteron = ob.wave_function(H_initial, U=U_matrix_unitless)
    
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


    # --- Plot n_lambda^d(q) --- #
    
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )

    ax.semilogy(q_array, n_d_total_array, color='black', linestyle='solid',
                label='total')
    ax.semilogy(q_array, n_d_1, color='blue', linestyle='dotted', 
                label='First term') 
    ax.semilogy(q_array, abs(n_d_23), color='purple', linestyle='dotted', 
                label='|Second+Third Term|') 
    ax.semilogy(q_array, n_d_4, color='red', linestyle='dotted', 
                label='Fourth Term') 
    ax.semilogy(q_array, n_d_exact, color='green', linestyle='dashed', 
                label='exact') 
    ax.set_ylabel(r'$<n^{\lambda}_d(q)>$' + ' [fm' + r'$^3$' + ']')
    ax.set_xlim( [min(q_array), 4] )
    ax.set_ylim( [1e-5, 1e3])
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    lambda_label_location = 'center right'
    anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
                                  frameon=False)
    ax.add_artist(anchored_text)

    plt.show()