#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: pmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 31, 2021
# 
# Calculates SRG-evolved pair momentum distributions for nuclei assuming the
# evolved wave function is given by a free Fermi gas. (Note, 'pmd' stands for
# pair momentum distribution.) Combine these functions with lda.py for
# nuclear-averaged momentum distributions. This script is based off a
# combination of old and current scripts: snmd.py, pmd_simple_test.py, and
# pmd_deuteron_test.py with the latter two in Old_codes.
#
# Revision history:
#   xx/xx/xx --- ...
#
#------------------------------------------------------------------------------


import numpy as np
# Scripts made by A.T.
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class pair_momentum_distributions(object):
    
    
    def __init__(self, kvnn, channels, lamb, kmax=0.0, kmid=0.0, ntot=0,
                 Q_dependence='Zero'):
        """
        Saves momentum arrays, grids, and matrix elements of \delta U and
        \delta U^{\dagger} given the input potential and SRG \lambda.
        Evaluates and saves the pp and pn matrix elements.
        
        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        channels : tuple
            Partial wave channels to include in the calculation.
        lamb : float
            SRG evolution parameter lambda [fm^-1].
        kmax : float, optional
            Maximum value in the momentum mesh [fm^-1].
        kmid : float, optional
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int, optional
            Number of momentum points in mesh.
        Q_dependence : str, optional
            How to deal with Q-dependence. Either set equal to zero ('Zero')
            or integrate out the Q-dependence ('Integrate').
            
        Notes
        -----
        Currently we only have Q_dependence = 'Zero' as an option. Still need
        to do the 'Integrate' option.
            
        """
        
        # --- Set up --- #
        
        # Save highest allowed L based on input channels
        highest_L = 0
        for channel in channels:
            next_L = vnn.channel_L_value(channel)
            if next_L > highest_L:
                highest_L = next_L
        
        # Load and save momentum arrays for integration
        k_array, k_weights = vnn.load_momentum(kvnn, '1S0')
        ntot = len(k_array)
        self.k_array, self.k_weights, self.ntot = k_array, k_weights, ntot
        self.k_integration_measure = k_weights * k_array**2
        
        # For dividing out momenta/weights
        factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
        # For coupled-channel matrices
        factor_array_cc = np.concatenate( (factor_array, factor_array) )
        
        
        # --- Evaluate matrix elements for each channel --- #
        
        # Initialize pp and pn matrix elements
        deltaU_pp = np.zeros( (ntot, ntot) ) # \delta U linear term
        deltaU_pn = np.zeros( (ntot, ntot) )
        deltaU2_pp = np.zeros( (ntot, ntot) ) # \delta U \delta U^\dagger term
        deltaU2_pn = np.zeros( (ntot, ntot) )
        
        # Allowed channels for pp (and nn) up through the D-waves
        pp_channels = ('1S0', '3P0', '3P1', '3P2', '1D2')
        
        # Loop over channels and evaluate matrix elements
        for channel in channels:

            # Load SRG transformation
            H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
            H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
                                             method='srg', generator='Wegner', 
                                             lamb=lamb)
            # Load U(k, k') [unitless]
            U_matrix_unitless = SRG_unitary_transformation(H_initial,
                                                           H_evolved)

            # Isolate 2-body term and convert to fm^3
            if vnn.coupled_channel(channel):
                I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
                row, col = np.meshgrid(factor_array_cc, factor_array_cc)
            else:
                I_matrix_unitless = np.eye(ntot, ntot)
                row, col = np.meshgrid(factor_array, factor_array)
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # 2J+1 factor
            J = int( channel[-1] )
            
            # Add to the pp and pn terms
            # Coupled-channel
            if vnn.coupled_channel(channel):
                    
                # First L of coupled-channel
                # Isospin CG's=1/\sqrt(2) for pn
                deltaU_pn += (2*J+1)/2 * delta_U_matrix[:ntot, :ntot]
                deltaU2_pn += (2*J+1)/4 * ( delta_U_matrix[:ntot, :ntot]**2 + \
                                            delta_U_matrix[:ntot, ntot:]**2 )

                # Isospin CG's=1 for pp
                if channel in pp_channels:
                    deltaU_pp += (2*J+1) * delta_U_matrix[:ntot, :ntot]
                    deltaU2_pp += (2*J+1) * ( delta_U_matrix[:ntot, :ntot]**2 \
                                            + delta_U_matrix[:ntot, ntot:]**2 )
                    
                # Decide whether to add second L based on highest allowed
                # L value (e.g., 0 + 2 <= 2 meaning we include the 3D1-3D1
                # part of the coupled 3S1-3D1 channel if we input D-waves
                # in channels)
                if vnn.channel_L_value(channel) + 2 <= highest_L:
                    deltaU_pn += (2*J+1)/2 * delta_U_matrix[ntot:, ntot:]
                    deltaU2_pn += (2*J+1)/4 * ( \
                                  delta_U_matrix[ntot:, :ntot]**2 + \
                                  delta_U_matrix[ntot:, ntot:]**2 )
                        
                    if channel in pp_channels:
                        deltaU_pp += (2*J+1) * delta_U_matrix[ntot:, ntot:]
                        deltaU2_pp += (2*J+1) * ( \
                                            delta_U_matrix[ntot:, :ntot]**2 + \
                                            delta_U_matrix[ntot:, ntot:]**2 )
            
            else:
                
                # Isospin CG's=1/\sqrt(2) for pn
                deltaU_pn += (2*J+1)/2 * delta_U_matrix
                deltaU2_pn += (2*J+1)/4 * delta_U_matrix**2
                
                # Isospin CG's=1 for pp
                if channel in pp_channels:
                    deltaU_pp += (2*J+1) * delta_U_matrix
                    deltaU2_pp += (2*J+1) * delta_U_matrix**2

        # Save \delta U matrix elements (these are (ntot, ntot) arrays)
        # < k | \delta U | k' >
        self.deltaU_pp = deltaU_pp
        self.deltaU_pn = deltaU_pn
        # < k | \delta U \delta U^\dagger | k' >
        self.deltaU2_pp = deltaU2_pp
        self.deltaU2_pn = deltaU2_pn


    def n_lambda_pp(self, q, kF):
        """
        Pair momentum distribution for proton-proton only where kF specifies
        the proton Fermi momentum.
        
        Parameters
        ----------
        q : float
            Relative momentum value [fm^-1].
        kF : float
            Fermi momentum [fm^-1].
            
        Return
        ------
        output : float
            Pair momentum distribution [unitless] evaluated at momentum q
            [fm^-1].
            
        """
        
        # Find index of q in k_array
        q_index = find_q_index(q, self.k_array)
        
        # Split into low- and high-q terms
        
        # Low-q terms
        # Term 1: 1 * n(q) * 1
        # Term 2: \deltaU * n(q) * 1
        # Term 3: 1 * n(q) * \deltaU^\dagger = Term 2
        # \delta U term : Term 2 + Term 3
        if q < kF: # \theta(kF - q)*\theta(kF - q)
            
            term_1 = 2 # \sum_{\sigma, \sigma'} 1/2 = 2
            
            deltaU_factor = 1/(4*np.pi) * 2/np.pi
            # delta U evaluated at q
            term_deltaU = deltaU_factor * self.deltaU_pp[q_index, q_index]
            
        # q > kF_1
        else:
            
            term_1 = 0
            term_deltaU = 0
            
            
        ## LEFT OFF HERE ##
        
        # High-q term: \deltaU * n(q) * \deltaU^\dagger
        
        # Approximate abs(q_vec - K_vec/2) as \sqrt(q^2 + K^2/4)
        q_K_array = np.sqrt( q**2 + self.K_array**2/4 ) # (Ntot, 1)
        
        # We can only interpolate \deltaU(k,k') up to k_max
        # Do not allow this next array to go above k_max - this may impose
        # a cutoff on the K integration
        q_K_array_cutoff = q_K_array[ q_K_array < max(self.k_array) ]
        K_cutoff_index = len(q_K_array_cutoff)
        
        # Create a grid for evaluation of \delta( k, abs(q_vec - K_vec/2) ) *
        # \delta U^\dagger( abs(q_vec - K_vec/2), k )
        k_grid, q_K_grid, = np.meshgrid(self.k_array, q_K_array_cutoff,
                                        indexing='ij')
        
        # Evaluate \delta( k, abs(q_vec - K_vec/2) ) *
        # \delta U^\dagger( abs(q_vec - K_vec/2), k ) for pp and pn (or nn
        # and np if kF_1 corresponds to a neutron)
        deltaU2_pp_array = self.deltaU2_pp_func.ev(k_grid, q_K_grid)
        deltaU2_pn_array = self.deltaU2_pn_func.ev(k_grid, q_K_grid)
        
        # Build x-dependent part and integrate with resized K
        
        # These are (ntot, Ktot, xtot) arrays of
        # \theta( k_F - \abs( 1/2*K_vec +/- k_vec ) ) for every total momentum
        # K, relative momentum k, and angle x where sign specifies the sign of
        # k_vec
        theta_kF1_K_plus_k_x = self.theta_K_k(kF_1, sign=1)
        theta_kF1_K_minus_k_x = self.theta_K_k(kF_1, sign=-1)
        theta_kF2_K_plus_k_x = self.theta_K_k(kF_2, sign=1)
        theta_kF2_K_minus_k_x = self.theta_K_k(kF_2, sign=-1)
            
        # Integrate over x first
        # This sum collapses theta_K_k_x to a matrix dependent only on K and 
        # k: (ntot, Ktot, xtot) -> (ntot, K_cutoff_index)
            
        # \int dx/2 \theta_pp (or \theta_nn)
        theta_kF1_kF1_K_k = ( np.sum( 
                              theta_kF1_K_plus_k_x * theta_kF1_K_minus_k_x,
                              axis=-1 ) )[:, :K_cutoff_index] / 2
        # \int dx/2 \theta_pn (or \theta_np)
        theta_kF1_kF2_K_k = ( np.sum( 
                              theta_kF1_K_plus_k_x * theta_kF2_K_minus_k_x +
                              theta_kF2_K_plus_k_x * theta_kF1_K_minus_k_x,
                              axis=-1 ) )[:, :K_cutoff_index] / 2

        # Build K integrand (ntot, K_cutoff_index) array spliting pp and np
        # contributions (or nn and np)
        integrand_k_K = self.K_integration_measure[:, :K_cutoff_index] * (\
                            deltaU2_pp_array * theta_kF1_kF1_K_k + \
                            deltaU2_pn_array * theta_kF1_kF2_K_k )

        # Integrate over K and build k integrand
        integrand_k = np.sum( integrand_k_K, axis=-1 ) * \
                      self.k_integration_measure
                      
        # Integrate over k
        deltaU2_factor = 1/2 * 1/(2*np.pi)**6 * (2/np.pi)**2
        term_deltaU2 = deltaU2_factor * np.sum(integrand_k)

        # Return total
        return term_1 + term_deltaU + term_deltaU2
        
    
    def n_lambda_pn(self):
        
        return None
    
    
    def n_lambda_total(self, q, kF_p, kF_n):
        
        # total = self.n_lambda_pp(q, kF_p) + self.n_lambda_pn(q, kF_p, kF_n) + \
        #         self.n_lambda_pp(q, kF_n) + self.n_lambda_np(q, kF_n, kF_p)
            
        # Add up each contribution
        return None