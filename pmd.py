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
# Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
# screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!
#
# Revision history:
#   04/26/21 --- Correcting overall factors in front of \delta U and
#                \delta U^2 terms.
#
#------------------------------------------------------------------------------


import numpy as np
# Scripts made by A.T.
from Misc.integration import gaussian_quadrature_mesh
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from scipy.interpolate import RectBivariateSpline
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
        k_array, k_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
        if ntot == 0:
            ntot = len(k_array) # Make sure ntot is the number of k-points
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
                deltaU2_pn += (2*J+1)/2 * ( delta_U_matrix[:ntot, :ntot]**2 + \
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
                    deltaU2_pn += (2*J+1)/2 * ( \
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
                deltaU2_pn += (2*J+1)/2 * delta_U_matrix**2
                
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
        
        # Interpolate < k | \delta U \delta U^\dagger | k' > for integration
        # over \int dk k^2 from 0 to kF(r)
        self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pp)
        self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pn)
        

    def n_lambda_pp(self, q, kF, contributions='total'):
        """
        Pair momentum distribution for proton-proton only where kF specifies
        the proton Fermi momentum. For neutron-neutron use kF corresponding to
        a neutron.
        
        Parameters
        ----------
        q : float
            Relative momentum value [fm^-1].
        kF : float
            Fermi momentum [fm^-1].
        contributions : str, optional
            Option to return different contributions to the momentum
            distribution.
            1. Default is 'total' which only returns the total momentum
               distribution.
            2. Specify 'q_contributions' for 1, \delta U, and 
               \delta U \delta U^\dagger along with the total.
            
        Return
        ------
        output : float or tuple
            Pair momentum distribution [unitless] evaluated at momentum q
            [fm^-1]. (Note, this function will return a tuple of floats if
            contributions is not 'total'.)
            
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
            
            # deltaU_factor =  2/np.pi * 1/(4*np.pi) * 2**2
            deltaU_factor =  1 * 2 * 2/np.pi * 1/(4*np.pi)
            # Need to include 1/\sqrt(2) factors from Dickhoff expression
            # Matches what is in the text of appendix
            # delta U evaluated at q
            term_deltaU = deltaU_factor * self.deltaU_pp[q_index, q_index]
            
        # q > kF_1
        else:
            
            term_1 = 0
            term_deltaU = 0

        # High-q term: \deltaU * n(q) * \deltaU^\dagger
        
        # # Evaluate \delta(k, q) * \delta U^\dagger(q, k) and multiply by
        # # dk k^2
        # integrand_k = self.deltaU2_pp[:, q_index] * self.k_integration_measure

        # # Index of where to stop the integral in the momentum mesh k_array
        # # (This is from \theta(kFp-k))
        # kF_cutoff = find_q_index(kF, self.k_array)
        # # kF_cutoff = find_q_index(kF, self.k_array) + 1
                      
        # # Integrate over k
        # deltaU2_factor = 1/4 * (2/np.pi)**2 * 1/(4*np.pi) * 2**4
        # # deltaU2_factor = 1/4 * (2/np.pi)**2 * 1/(4*np.pi) * 2**2
        # term_deltaU2 = deltaU2_factor * np.sum( integrand_k[:kF_cutoff] )
        
        # TESTING
        if kF >= 1.2:
            ntot_k = 60
        elif 1.2 > kF >= 1.0:
            ntot_k = 50
        elif 1.0 > kF >= 0.8:
            ntot_k = 40
        elif 0.8 > kF >= 0.6:
            ntot_k = 30
        elif 0.6 > kF >= 0.4:
            ntot_k = 20
        else:
            ntot_k = 10
        k_array, k_weights = gaussian_quadrature_mesh(kF, ntot_k)
        k_integration_measure = k_weights * k_array**2
        
        # Evaluate \delta(k, q) * \delta U^\dagger(q, k) and multiply by
        # dk k^2
        integrand_k = self.deltaU2_pp_func.ev(k_array, q) * \
                      k_integration_measure
        
        # Integrate over k
        # deltaU2_factor = 1/4 * (2/np.pi)**2 * 1/(4*np.pi) * 2**4
        # See comment by deltaU factor
        deltaU2_factor = 1/4 * 2**2 * (2/np.pi)**2 * 1/(4*np.pi)
        term_deltaU2 = deltaU2_factor * np.sum(integrand_k)
        
        # Add up each term for total
        total = term_1 + term_deltaU + term_deltaU2

        # Return contributions and total or just total
        if contributions == 'q_contributions':
            return total, term_1, term_deltaU, term_deltaU2
        else: # Default
            return total
        
    
    def n_lambda_pn(self, q, kF_1, kF_2, contributions='total'):
        """
        Pair momentum distribution for proton-neutron only where kF_1 (kF_2)
        specifies the proton (neutron) Fermi momentum. Flip the kF arguments
        for neutron-proton.
        
        Parameters
        ----------
        q : float
            Relative momentum value [fm^-1].
        kF_1 : float
            Fermi momentum [fm^-1] for the corresponding nucleon momentum
            distribution.
        kF_2 : float
            Fermi momentum [fm^-1] for the correlated nucleon. If kF_1
            corresponds to a proton, then kF_2 corresponds to a neutron (and
            vice versa).
        contributions : str, optional
            Option to return different contributions to the momentum
            distribution.
            1. Default is 'total' which only returns the total momentum
               distribution.
            2. Specify 'q_contributions' for 1, \delta U, and 
               \delta U \delta U^\dagger along with the total.
            
        Return
        ------
        output : float or tuple
            Pair momentum distribution [unitless] evaluated at momentum q
            [fm^-1]. (Note, this function will return a tuple of floats if
            contributions is not 'total'.)
            
        """
        
        # In the Q=0 case, one of the \theta's is more restrictive
        kF = min(kF_1, kF_2) # Take smaller kF and evaluate
        
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
            
            # deltaU_factor =  2/np.pi * 1/(4*np.pi) * 2**2
            deltaU_factor =  1 * 2 * 2/np.pi * 1/(4*np.pi)
            # Need to include 1/\sqrt(2) factors from Dickhoff expression
            # Matches what is in the text of appendix
            # delta U evaluated at q
            term_deltaU = deltaU_factor * self.deltaU_pn[q_index, q_index]
            
        # q > kF_1
        else:
            
            term_1 = 0
            term_deltaU = 0

        # High-q term: \deltaU * n(q) * \deltaU^\dagger
        
        # # Evaluate \delta(k, q) * \delta U^\dagger(q, k) and multiply by
        # # dk k^2
        # integrand_k = self.deltaU2_pn[:, q_index] * self.k_integration_measure

        # # Index of where to stop the integral in the momentum mesh k_array
        # # (This is from \theta(kFp-k))
        # kF_cutoff = find_q_index(kF, self.k_array)
                      
        # # Integrate over k
        # deltaU2_factor = 1/4 * (2/np.pi)**2 * 1/(4*np.pi) * 2**4 * 2
        # # Last factor of 2 is for \theta^p \theta^n + \theta^p \theta^n
        # term_deltaU2 = deltaU2_factor * np.sum( integrand_k[:kF_cutoff] )
        
        # TESTING
        if kF >= 1.2:
            ntot_k = 60
        elif 1.2 > kF >= 1.0:
            ntot_k = 50
        elif 1.0 > kF >= 0.8:
            ntot_k = 40
        elif 0.8 > kF >= 0.6:
            ntot_k = 30
        elif 0.6 > kF >= 0.4:
            ntot_k = 20
        else:
            ntot_k = 10
        k_array, k_weights = gaussian_quadrature_mesh(kF, ntot_k)
        k_integration_measure = k_weights * k_array**2
        
        # Evaluate \delta(k, q) * \delta U^\dagger(q, k) and multiply by
        # dk k^2
        integrand_k = self.deltaU2_pn_func.ev(k_array, q) * \
                      k_integration_measure
        
        # Integrate over k
        # deltaU2_factor = 1/4 * (2/np.pi)**2 * 1/(4*np.pi) * 2**4
        # See comment above near \delta U factor
        deltaU2_factor = 1/4 * 2**2 * (2/np.pi)**2 * 1/(4*np.pi)
        term_deltaU2 = deltaU2_factor * np.sum(integrand_k)
        
        # Add up each term for total
        total = term_1 + term_deltaU + term_deltaU2

        # Return contributions and total or just total
        if contributions == 'q_contributions':
            return total, term_1, term_deltaU, term_deltaU2
        else: # Default
            return total
        
        
    def select_number_integration_points(self, k_max, k_min=0.0):
        """
        Select the number of integration points over momenta k given the upper
        limit of integration. Assumes Gaussian quadrature integration.
        
        Parameters
        ----------
        k_max : float
            Upper limit of integration over momentum [fm^-1].
        k_min : float, optional
            Lower limit of integration over momentum [fm^-1]. Default is zero.
            
        Returns
        -------
        ntot_k : int
            Number of integration points in Gaussian quadrature mesh.
        
        """
        
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
    
    
    def n_1(self, q, kF):
        """
        Evaluates one of the \theta functions in U n(q) U^\dagger ~ I n(q) I
        which gives 2 \theta( kF(r) - q ). The full first term is given by
            2 \int d3r \theta( kF_1(r) - q ) \int d3r' \theta( kF_2(r) - q ),
        where the factor of two comes from \sum_{\sigma, \sigma'} 1/2 = 2.
        
        Parameters
        ----------
        q : float
            Single-nucleon momentum [fm^-1].
        kF : float
            Fermi momentum [fm^-1].
            
        Returns
        -------
        output : float
            Momentum distribution from I term before integration over
            \int dr r^2.
        
        """
        
        # This is simply a \theta function
        if q < kF:
            
            return 1
        
        else:
            
            return 0
        
        
    def n_deltaU(self, q, kF_1, kF_2=0.0):
        """
        Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
        we are combining \delta U and \delta U^\dagger, hence the factor of 2.
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        kF_1 : float
            Fermi momentum [fm^-1] for the corresponding nucleon momentum
            distribution.
        kF_2 : float, optional
            Fermi momentum [fm^-1] for the correlated nucleon. If kF_1
            corresponds to a proton, then kF_2 corresponds to a neutron (and
            vice versa).
            
        Returns
        -------
        output : float
            Momentum distribution from \delta U term before integration over
            \int dr r^2.
            
        """
        
        # LEFT OFF HERE
        
        # In the Q=0 case, one of the \theta's is more restrictive (if no kF_2
        # is provided then take kF = kF_1 -> pp or nn distribution)
        if kF_2 == 0:
            kF = kF_1
        else:
            kF = min(kF_1, kF_2) # Take smaller kF and evaluate
        
        # Find index of q in k_array
        q_index = find_q_index(q, self.k_array)

        # Restricted by \theta function with smaller kF value
        if q < kF: # \theta(kF - q) * \theta(kF - q)
        
            # Evaluate pp or pn < k | \delta U | k > matrix elements (or nn
            # and np if kF_1 corresponds to a neutron)
            
            deltaU_factor =  1 * 2 * 2/np.pi * 1/(4*np.pi)
            # Need to include 1/\sqrt(2) factors from Dickhoff expression
            # Matches what is in the text of appendix
            # delta U evaluated at q
            term_deltaU = deltaU_factor * self.deltaU_pn[q_index, q_index]

        else:
            
            return 0
    
    
    def n_lambda(self, q_array, r_array, rho_1_array, rho_2_array=np.empty(0)):
        """
        Pair momentum distribution where the nucleonic densities specify the
        nucleus.
        
        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        r_array : 1-D ndarray
            Coordinates array [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of r [fm^-3] for the first nucleon
            corresponding to \tau.
        rho_1_array : 1-D ndarray, optional
            Densities as a function of r [fm^-3] for the second nucleon
            corresponding to \tau'. If empty array, function assumes a proton-
            proton (or neutron-neutron) pair momentum distribution.
            
        Returns
        -------
        n_lambda : 2-D ndarray
            Array of contributions to the pair momentum distribution ordered
            according to total, 1, \delta U, and \delta^2 [fm^3].
            
        Notes
        -----
        Currently we only have Q_dependence = 'Zero' as an option. Still need
        to do the 'Integrate' option.
            
        """
        
        # Number of columns for output array (total, 1, \delta U, \delta U^2)
        axes = 4
            
        # Number of r_array points
        mtot = len(r_array)
        r2_grid, _ = np.meshgrid(r_array**2, np.zeros(axes), indexing='ij')
        dr = 0.1 # Spacing between r-points
            
        # Length of q_array
        ntot = len(q_array)
        
        # Evaluate n_\lambda(q, kF) contributions at each point in q_array,
        # rho_1_array, and rho_2_array
        n_lambda = np.zeros( (ntot, axes) )
        # Temporary n_lambda_array to store different r integrals in '1' term
        n_lambda_temp = np.zeros( (ntot, axes) )

        # Evaluate kF values at each point in r_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        # Calculate kF_2 array if pn or np distribution
        if rho_2_array.any():
            kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
        # Otherwise set kF_2 to zero where previous functions will correspond
        # to pp or nn distributions (depending on kF_1)
        else:
            kF2_array = np.zeros(mtot)
            
        # Loop over q
        for i, q in enumerate(q_array):
    
            # n_lambda contributions before integrating over r
            n_lambda_r = np.zeros( (mtot, axes) )
            
            # Loop over kF values
            for j, (kF_1, kF_2) in enumerate( zip(kF1_array, kF2_array) ):
                
                # Two r integrations for 1 * n(q) * 1 term
                term_1_tau = self.n_1(q, kF_1) # \tau nucleon
                if kF_2 == 0: # pp or nn distribution
                    term_1_taup = term_1_tau # Same kF value
                else: # Different kF value
                    term_1_taup = self.n_1(q, kF_2) # \tau' nucleon
                term_deltaU = self.n_deltaU(q, kF_1, kF_2)
                term_deltaU2 = self.n_deltaU2(q, kF_1, kF_2)
                
                n_lambda_r[j, :] = term_1_tau, term_1_taup, term_deltaU, \
                                   term_deltaU2

            # Integrate over 4\pi \int dr r^2 for each contribution (summing
            # over axis=0)
            n_lambda_temp[i, :] = 4*np.pi*dr * np.sum(r2_grid * n_lambda_r,
                                                      axis=0)
            
            # Construct n_lambda array for contributions according to total,
            # 1, \delta U, and \delta U^2 where the 1 term includes two r-
            # integrations
            n_lambda[i, 1] = 2 * n_lambda_temp[i, 0] * n_lambda_temp[i, 1]
            n_lambda[i, 2:] = n_lambda_temp[i, 2:]
            n_lambda[i, 0] = n_lambda[i, 1] + n_lambda[i, 2] + n_lambda[i, 3]
        
        # Return (ntot, 4) array of contributions to n_\lambda^\tau(q)
        return n_lambda