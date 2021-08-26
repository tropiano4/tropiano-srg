#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: pmd_test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 31, 2021
# 
# Calculates SRG-evolved pair momentum distributions for nuclei assuming the
# evolved wave function is given by HF treated in LDA. (Note, 'pmd' stands for
# pair momentum distribution.)
#
# Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
# screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!
#
# Revision history:
#   04/26/21 --- Correcting overall factors in front of \delta U and
#                \delta U^2 terms.
#   06/08/21 --- Corrected "I" term to have additional integration over R'.
#   06/22/21 --- Generalizing distribution to n(q, Q) from n(q, 0). Saved old
#                version as pmd_v1.py in Old_codes.
#   06/23/21 --- Speeding up code by switching from loops to np.sum() to do
#                integrations. Saved old version as pmd_v2.py in Old_codes. 
#   07/02/21 --- Added interpolation option for n_{\lambda}(q, Q).
#   07/15/21 --- Switched interpolation functions interp1d and 
#                RectBivariateSpline from cubic to linear since \delta U^2(k,q)
#                was returning negative values. Also, setting Q_max to 2 fm^-1
#                instead of max(kF1)+max(kF2) in file writing function.
#   08/18/21 --- Testing changes to the sharp cutoff in \theta(kF-q) and trying
#                block-diagonal SRG transformations.
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline
# Scripts made by A.T.
from densities import load_density
from figures import figures_functions as ff
from misc.integration import gaussian_quadrature_mesh
from potentials.vsrg_macos import vnn
from srg.srg_unitary_transformation import SRG_unitary_transformation


class pair_momentum_distributions(object):
    
    
    def __init__(self, kvnn, channels, lamb, kmax, kmid, ntot,
                 generator='Wegner', beta=0.0, interp=False):
        """
        Evaluates and saves the pp and pn matrix elements of \delta U and
        \delta U^{\dagger} given the input potential and SRG \lambda.
        
        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        channels : tuple
            Partial wave channels to include in the calculation.
        lamb : float
            SRG evolution parameter lambda [fm^-1].
        kmax : float
            Maximum value in the momentum mesh [fm^-1].
        kmid : float
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int
            Number of momentum points in mesh.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        beta : bool, optional
            \beta parameter [fm^-1] in Fermi function. If left zero, uses
            \theta functions as default to evaluate
              < F | n_\alpha | F > ~ \theta( kF(R) - k_\alpha ).
            Nonzero entries will replace the \theta function with a Fermi
            function 1 / ( 1 + exp^(-(kF-q)/\beta ) ). Low \beta approaches
            the \theta function.
        interp : bool, optional
            Option to use interpolated n_\lambda(q, Q) functions.
            
        """
        
        # Get relevant info for file and directory names
        # Part of data directory name
        self.kvnn = kvnn 
        # Part of file name
        self.channels = channels
        self.lamb = lamb
        self.kmax = kmax
        self.ntot = ntot
        self.generator = generator
        self.beta = beta

        if interp == False:
            
            # Save highest allowed L based on input channels
            highest_L = 0
            for channel in channels:
                next_L = vnn.channel_L_value(channel)
                if next_L > highest_L:
                    highest_L = next_L
        
            # Load and save momentum arrays
            k_array, k_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid,
                                                   ntot)
            # Save k_array for writing files
            self.k_array = k_array
        
            # For dividing out momenta/weights
            factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
            # For coupled-channel matrices
            factor_array_cc = np.concatenate( (factor_array, factor_array) )

            # Initialize pp and pn matrix elements
            # \delta U linear term
            deltaU_pp = np.zeros( (ntot, ntot) )
            deltaU_pn = np.zeros( (ntot, ntot) )
            # \delta U \delta U^\dagger term
            deltaU2_pp = np.zeros( (ntot, ntot) )
            deltaU2_pn = np.zeros( (ntot, ntot) )
        
            # Allowed channels for pp (and nn) up through the D-waves
            pp_channels = ('1S0', '3P0', '3P1', '3P2', '1D2')
        
            # Loop over channels and evaluate matrix elements
            for channel in channels:

                # Load SRG transformation
                H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid,
                                                 ntot)
                if generator == 'Block-diag':
                    # Take \lambda = 1 fm^-1 and set \Lambda_BD = input \lambda
                    H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid,
                                                     ntot, 'srg', generator,
                                                     1.0, lamb)
                else:
                    H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid,
                                                     ntot, 'srg', generator,
                                                     lamb)
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
                delta_U_matrix = delta_U_matrix_unitless / row / col # fm^3
            
                # 2J+1 factor
                J = int( channel[-1] )
            
                # Add to the pp and pn terms
                # Coupled-channel
                if vnn.coupled_channel(channel):
                    
                    # First L of coupled-channel
                    # Isospin CG's=1/\sqrt(2) for pn
                    deltaU_pn += (2*J+1)/2 * delta_U_matrix[:ntot, :ntot]
                    deltaU2_pn += (2*J+1)/2 * ( \
                                      delta_U_matrix[:ntot, :ntot]**2 + \
                                      delta_U_matrix[:ntot, ntot:]**2 )

                    # Isospin CG's=1 for pp
                    if channel in pp_channels:
                        deltaU_pp += (2*J+1) * delta_U_matrix[:ntot, :ntot]
                        deltaU2_pp += (2*J+1) * ( \
                                          delta_U_matrix[:ntot, :ntot]**2 \
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
                                              delta_U_matrix[ntot:, :ntot]**2 \
                                            + delta_U_matrix[ntot:, ntot:]**2 )
            
                else:
                
                    # Isospin CG's=1/\sqrt(2) for pn
                    deltaU_pn += (2*J+1)/2 * delta_U_matrix
                    deltaU2_pn += (2*J+1)/2 * delta_U_matrix**2
                
                    # Isospin CG's=1 for pp
                    if channel in pp_channels:
                        deltaU_pp += (2*J+1) * delta_U_matrix
                        deltaU2_pp += (2*J+1) * delta_U_matrix**2

            # Interpolate pp and pn < k | \delta U | k >
            self.deltaU_pp_func = interp1d( k_array, np.diag(deltaU_pp),
                                            kind='linear', bounds_error=False,
                                            fill_value='extrapolate' )
            self.deltaU_pn_func = interp1d( k_array, np.diag(deltaU_pn),
                                            kind='linear', bounds_error=False,
                                            fill_value='extrapolate' )
        
            # Interpolate pp and pn < k | \delta U \delta U^{\dagger} | k' > 
            self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                       deltaU2_pp, kx=1, ky=1)
            self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                       deltaU2_pn, kx=1, ky=1)


    def theta_I(self, q_mesh, Q_mesh, kF1_mesh, kF2_mesh):
        """
        Evaluates angle-average of \theta( kF1(R) - \abs(Q/2 + q) ) x
        \theta( kF2(R') - \abs(Q/2 - q) ). This function appears in the I term.

        Parameters
        ----------
        q_mesh : 4-D ndarray
            Relative momentum values [fm^-1].
        Q_mesh : 4-D ndarray
            C.o.M. momentum values [fm^-1].
        kF1_mesh : 4-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_mesh : 4-D ndarray
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R'.

        Returns
        -------
        theta_mesh : 4-D ndarray
            \theta function [unitless] evaluated for each q, Q, kF1(R), and
            kF2(R').

        """
        
        # Initialize 4-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_Q, self.ntot_R,
                                self.ntot_R) )
        
        # Evaluate each boolean case and use these to fill in the theta_mesh

        # Case 1: 2q+Q < 2kF1 and 2q+Q < 2kF2
        case_1 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
                 ( 2*q_mesh + Q_mesh <= 2*kF2_mesh )
        theta_mesh[case_1] = 1
        
        # Q = 0 case simplifies \theta functions
        # If Q = 0, skip the following cases
        if self.ntot_Q == 1:
            return theta_mesh
                
        # Case 2: 2q+Q > 2kF1 and 2q+Q > 2kF2 and 4q^2+Q^2 < 2(kF1^2+kF2^2)
        case_2 = ( 2*q_mesh + Q_mesh > 2*kF1_mesh ) * \
                 ( 2*q_mesh + Q_mesh > 2*kF2_mesh ) * \
                 ( 4*q_mesh**2 + Q_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
        theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*q_mesh**2 - \
                               Q_mesh**2 ) / (4*q_mesh*Q_mesh) )[case_2]
                            
        # Case 3: 2q+Q < 2kF2 and -4 < (4q^2 - 4kF1^2 + Q^2)/(qQ) < 4
        case_3 = ( 2*q_mesh + Q_mesh <= 2*kF2_mesh ) * \
                 ( -4 < (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) ) * \
                 ( (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) <= 4 )
        theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
                               (8*q_mesh*Q_mesh) )[case_3]
                
        # Case 4: 2q+Q < 2kF1 and -4 < (4q^2 - 4kF2^2 + Q^2)/(qQ) < 4
        case_4 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
                 ( -4 < (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) ) * \
                 ( (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) <= 4 )
        theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
                               (8*q_mesh*Q_mesh) )[case_4]
        
        # This is a (ntot_q, ntot_Q, ntot_R, ntot_R) size array
        return theta_mesh
    
    
    def theta_deltaU(self, q_mesh, Q_mesh, kF1_mesh, kF2_mesh):
        """
        Evaluates angle-average of \theta( kF1(R) - \abs(Q/2 + q) ) x
        \theta( kF2(R) - \abs(Q/2 - q) ). This function appears in the
        \delta U term.

        Parameters
        ----------
        q_mesh : 3-D ndarray
            Relative momentum values [fm^-1].
        Q_mesh : 3-D ndarray
            C.o.M. momentum values [fm^-1].
        kF1_mesh : 3-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_mesh : 3-D ndarray
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R.

        Returns
        -------
        theta_mesh : 3-D ndarray
            \theta function [unitless] evaluated for each q, Q, kF1(R), and
            kF2(R).

        """
        
        # Initialize 3-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_Q, self.ntot_R) )
        
        # Evaluate each boolean case and use these to fill in the theta_mesh

        # Case 1: 2q+Q < 2kF1 and 2q+Q < 2kF2
        case_1 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
                 ( 2*q_mesh + Q_mesh <= 2*kF2_mesh )
        theta_mesh[case_1] = 1
        
        # Q = 0 case simplifies \theta functions
        # If Q = 0, skip the following cases
        if self.ntot_Q == 1:
            return theta_mesh
                
        # Case 2: 2q+Q > 2kF1 and 2q+Q > 2kF2 and 4q^2+Q^2 < 2(kF1^2+kF2^2)
        case_2 = ( 2*q_mesh + Q_mesh > 2*kF1_mesh ) * \
                 ( 2*q_mesh + Q_mesh > 2*kF2_mesh ) * \
                 ( 4*q_mesh**2 + Q_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
        theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*q_mesh**2 - \
                               Q_mesh**2 ) / (4*q_mesh*Q_mesh) )[case_2]
                            
        # Case 3: 2q+Q < 2kF2 and -4 < (4q^2 - 4kF1^2 + Q^2)/(qQ) < 4
        case_3 = ( 2*q_mesh + Q_mesh <= 2*kF2_mesh ) * \
                 ( -4 < (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) ) * \
                 ( (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) <= 4 )
        theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
                               (8*q_mesh*Q_mesh) )[case_3]
                
        # Case 4: 2q+Q < 2kF1 and -4 < (4q^2 - 4kF2^2 + Q^2)/(qQ) < 4
        case_4 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
                 ( -4 < (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) ) * \
                 ( (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
                 (q_mesh*Q_mesh) <= 4 )
        theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
                               (8*q_mesh*Q_mesh) )[case_4]
        
        # This is a (ntot_q, ntot_Q, ntot_R) size array
        return theta_mesh
    
    
    def theta_deltaU2(self, Q_mesh, kF1_mesh, kF2_mesh, k_mesh):
        """
        Evaluates angle-average of \theta( kF1(R) - \abs(Q/2 + k) ) x
        \theta( kF2(R) - \abs(Q/2 - k) ). This function appears in the
        \delta U \delta U^\dagger term.

        Parameters
        ----------
        Q_mesh : 4-D ndarray
            C.o.M. momentum values [fm^-1].
        kF1_mesh : 4-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_mesh : 4-D ndarray
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R.
        k_mesh : 4-D ndarray
            Relative momentum values [fm^-1].

        Returns
        -------
        theta_mesh : 4-D ndarray
            \theta function [unitless] evaluated for each q, Q, kF1(R), kF2(R),
            and k. Note, this function does not depend on q but we need the
            array to match the size of \delta U \delta U^\dagger matrix
            elements.

        """
        
        # Initialize 4-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_Q, self.ntot_R,
                                self.ntot_k) )
        
        # Evaluate each boolean case and use these to fill in the theta_mesh

        # Case 1: 2k+Q < 2kF1 and 2k+Q < 2kF2
        case_1 = ( 2*k_mesh + Q_mesh <= 2*kF1_mesh ) * \
                 ( 2*k_mesh + Q_mesh <= 2*kF2_mesh )
        theta_mesh[case_1] = 1
        
        # Q = 0 case simplifies \theta functions
        # If Q = 0, skip the following cases
        if self.ntot_Q == 1:
            return theta_mesh
                
        # Case 2: 2k+Q > 2kF1 and 2k+Q > 2kF2 and 4k^2+Q^2 < 2(kF1^2+kF2^2)
        case_2 = ( 2*k_mesh + Q_mesh > 2*kF1_mesh ) * \
                 ( 2*k_mesh + Q_mesh > 2*kF2_mesh ) * \
                 ( 4*k_mesh**2 + Q_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
        theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*k_mesh**2 - \
                               Q_mesh**2 ) / (4*k_mesh*Q_mesh) )[case_2]
                            
        # Case 3: 2k+Q < 2kF2 and -4 < (4k^2 - 4kF1^2 + Q^2)/(kQ) < 4
        case_3 = ( 2*k_mesh + Q_mesh <= 2*kF2_mesh ) * \
                 ( -4 < (4*k_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
                 (k_mesh*Q_mesh) ) * \
                 ( (4*k_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
                 (k_mesh*Q_mesh) <= 4 )
        theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (Q_mesh-2*k_mesh)**2 ) / \
                               (8*k_mesh*Q_mesh) )[case_3]
                
        # Case 4: 2k+Q < 2kF1 and -4 < (4k^2 - 4kF2^2 + Q^2)/(kQ) < 4
        case_4 = ( 2*k_mesh + Q_mesh <= 2*kF1_mesh ) * \
                 ( -4 < (4*k_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
                 (k_mesh*Q_mesh) ) * \
                 ( (4*k_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
                 (k_mesh*Q_mesh) <= 4 )
        theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (Q_mesh-2*k_mesh)**2 ) / \
                               (8*k_mesh*Q_mesh) )[case_4]
        
        # This is a (ntot_q, ntot_Q, ntot_R, ntot_k) size array
        return theta_mesh
    
    
    def fermi_function_low_q(self, q_mesh, kF1_mesh, kF2_mesh):
        # Fermi function for I and \delta U terms
        
        # ff1_mesh = 1 / ( 1 + np.exp( -self.beta * (kF1_mesh**2 - q_mesh**2) ) )
        # ff2_mesh = 1 / ( 1 + np.exp( -self.beta * (kF2_mesh**2 - q_mesh**2) ) )
        ff1_mesh = 1 / ( 1 + np.exp( -(kF1_mesh - q_mesh) / self.beta ) )
        ff2_mesh = 1 / ( 1 + np.exp( -(kF2_mesh - q_mesh) / self.beta ) )
        
        # Dimensions will match dimensions of input mesh grids
        return ff1_mesh * ff2_mesh
    
    
    def fermi_function_high_q(self, kF1_mesh, kF2_mesh, k_mesh):
        # Fermi function for I and \delta U terms
        
        # ff1_mesh = 1 / ( 1 + np.exp( -self.beta * (kF1_mesh**2 - k_mesh**2) ) )
        # ff2_mesh = 1 / ( 1 + np.exp( -self.beta * (kF2_mesh**2 - k_mesh**2) ) )
        ff1_mesh = 1 / ( 1 + np.exp( -(kF1_mesh - k_mesh) / self.beta ) )
        ff2_mesh = 1 / ( 1 + np.exp( -(kF2_mesh - k_mesh) / self.beta ) )
        
        # Dimensions will match dimensions of input mesh grids
        return ff1_mesh * ff2_mesh
    
    
    def n_I(self, q_array, Q_array, R_array, dR, kF1_array, kF2_array):
        """
        Evaluates the I term in U n(q) U^\dagger ~ I n(q) I.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R'.

        Returns
        -------
        output : 2-D ndarray
            Momentum distribution [fm^6] from I term as a function of q and Q.

        """
        
        # Initialize 4-D meshgrids (q, Q, R, R')
        q_mesh, Q_mesh, R_mesh, Rp_mesh = np.meshgrid(q_array, Q_array,
                                                      R_array, R_array,
                                                      indexing='ij')
        # Get 4-D kF1(R) and kF2(R') meshes
        _, _, kF1_mesh, kF2_mesh = np.meshgrid(q_array, Q_array, kF1_array,
                                               kF2_array, indexing='ij')
        
        # Evaluate angle-average of \theta-functions in I term
        if self.beta == 0.0:
            theta_mesh = self.theta_I(q_mesh, Q_mesh, kF1_mesh, kF2_mesh)
        else:
            theta_mesh = self.fermi_function_low_q(q_mesh, kF1_mesh, kF2_mesh)
        
        # Calculate R' integrand (ntot_q, ntot_Q, ntot_R, ntot_R)
        integrand_Rp = theta_mesh * Rp_mesh**2 * dR * R_mesh**2 * dR
        
        # Integrate over R' leaving R integrand (ntot_q, ntot_Q, ntot_R)
        integrand_R = 4*np.pi * np.sum(integrand_Rp, axis=-1)
        
        # Integrate over R
        # This is a (ntot_q, ntot_Q) size array
        # Factor of 2 is overall factor
        return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    
    def n_deltaU(self, q_array, Q_array, R_array, dR, kF1_array, kF2_array,
                 pair):
        """
        Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
        we are combining \delta U and \delta U^\dagger, hence the factor of 2.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R.
        pair : str
            Type of pair momentum distribution ('pp' or 'pn'). This will tell
            the function to either set isospin CGs to 1 or 1/\sqrt(2) (no need
            to specify 'nn' or 'np').

        Returns
        -------
        output : 2-D ndarray
            Momentum distribution [fm^6] from \delta U term as a function of q
            and Q.

        """
        
        # Initialize 3-D meshgrids (q, Q, R)
        q_mesh, Q_mesh, R_mesh = np.meshgrid(q_array, Q_array, R_array,
                                             indexing='ij')
        
        # Get 3-D kF1(R) and kF2(R) meshes
        _, _, kF1_mesh = np.meshgrid(q_array, Q_array, kF1_array,
                                     indexing='ij')
        _, _, kF2_mesh = np.meshgrid(q_array, Q_array, kF2_array,
                                     indexing='ij')
        
        # Evaluate angle-average of \theta-functions in \delta U term
        if self.beta == 0.0:
            theta_mesh = self.theta_deltaU(q_mesh, Q_mesh, kF1_mesh, kF2_mesh)
        else:
            theta_mesh = self.fermi_function_low_q(q_mesh, kF1_mesh, kF2_mesh)
        
        # Evaluate < q | \delta U | q >
        if pair == 'pp':
            deltaU_mesh = self.deltaU_pp_func(q_mesh)
        elif pair == 'pn':
            deltaU_mesh = self.deltaU_pn_func(q_mesh)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 1
        # Factor of 2 from \delta U + \delta U^\dagger
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        # 1/(4\pi) for averaging over \int d\Omega_q
        deltaU_factor =  2 * 2/np.pi * (2*np.pi)**3/(4*np.pi)
        
        # Calculate R integrand (ntot_q, ntot_Q, ntot_R)
        integrand_R = deltaU_factor * deltaU_mesh * theta_mesh * R_mesh**2 * dR
        
        # Integrate over R
        # This is a (ntot_q, ntot_Q) size array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
        
    
    def n_deltaU2(self, q_array, Q_array, R_array, dR, kF1_array, kF2_array,
                  pair):
        """
        Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R.
        pair : str
            Type of pair momentum distribution ('pp' or 'pn'). This will tell
            the function to either set isospin CGs to 1 or 1/\sqrt(2) (no need
            to specify 'nn' or 'np').

        Returns
        -------
        output : 2-D ndarray
            Momentum distribution [fm^6] from \delta U \delta U^\dagger term as
            a function of q and Q.

        """
        
        # Set number of k points for integration over k
        self.ntot_k = 40
        
        # Initialize 4-D meshgrids (q, Q, R, k) with k values equal to 0
        q_mesh, Q_mesh, R_mesh, k_mesh = np.meshgrid(q_array, Q_array, R_array,
                                                     np.zeros(self.ntot_k),
                                                     indexing='ij')
        
        # Get 4-D kF1(R) and kF2(R) meshes and initialize k weights mesh
        _, _, kF1_mesh, dk_mesh = np.meshgrid(q_array, Q_array, kF1_array,
                                              np.zeros(self.ntot_k),
                                              indexing='ij')
        _, _, kF2_mesh, _ = np.meshgrid(q_array, Q_array, kF2_array,
                                        np.zeros(self.ntot_k), indexing='ij')

        # Loop over q, Q, and R to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iQ, Q in enumerate(Q_array):
            for iR, R in enumerate(R_array):
                
                kF1, kF2 = kF1_array[iR], kF2_array[iR]
                
                # Minimum kF value
                kF_min = min(kF1, kF2)
                
                # Lower limit of integration
                k_min = max(Q/2 - kF_min, 0)
                
                # Upper limit of integration
                if Q**2/4 < (kF1**2 + kF2**2)/2:
                    k_max = min( np.sqrt( ( kF1**2 + kF2**2 )/2 - Q**2/4 ),
                                 kF_min + Q/2 )
                else:
                    k_max = kF_min + Q/2
                    
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max,
                                                              self.ntot_k,
                                                              xmin=k_min)
                
                # Fill in k_mesh and dk_mesh given the specific k_array
                for iq in range(self.ntot_q):
                    k_mesh[iq, iQ, iR, :] = k_array
                    dk_mesh[iq, iQ, iR, :] = k_weights
                     
        # Evaluate angle-average of \theta-functions in \delta U^2 term
        if self.beta == 0.0:
            theta_mesh = self.theta_deltaU2(Q_mesh, kF1_mesh, kF2_mesh, k_mesh)
        else:
            theta_mesh = self.fermi_function_high_q(kF1_mesh, kF2_mesh, k_mesh)
        
        # Evaluate < k | \delta U | q > < q | \delta U^\dagger | k >
        if pair == 'pp':
            deltaU2_mesh = self.deltaU2_pp_func.ev(k_mesh, q_mesh)
        elif pair == 'pn':
            deltaU2_mesh = self.deltaU2_pn_func.ev(k_mesh, q_mesh)

        # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
        # combine to give 1
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        # 1/(4\pi) for averaging over \int d\Omega_q
        deltaU2_factor =  (2/np.pi)**2 * (2*np.pi)**3/(4*np.pi)
        
        # Calculate k integrand (ntot_q, ntot_Q, ntot_R, ntot_k)
        integrand_k = deltaU2_factor * deltaU2_mesh * theta_mesh * \
                      k_mesh**2 * dk_mesh * R_mesh**2 * dR
        
        # Integrate over k leaving R integrand (ntot_q, ntot_Q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R
        # This is a (ntot_q, ntot_Q) size array
        return 4*np.pi * np.sum(integrand_R, axis=-1)


    def n_total(self, q_array, Q_array, R_array, dR, rho_1_array,
                rho_2_array=np.empty(0)):
        """
        Pair momentum distribution where the nucleonic densities specify the
        nucleus and distribution type (e.g., O16 and pn).

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the first nucleon
            corresponding to \tau.
        rho_2_array : 1-D ndarray, optional
            Densities as a function of R [fm^-3] for the second nucleon
            corresponding to \tau'. If an empty array is input, the function
            assumes a proton-proton (or neutron-neutron) pair momentum
            distribution relying only on rho_1_array.

        Returns
        -------
        n_total : 2-D ndarray
            Pair momentum distribution [fm^6] for each q and Q.

        """
        
        # Save lengths of q_array, Q_array, and R_array
        self.ntot_q = len(q_array)
        self.ntot_Q = len(Q_array)
        self.ntot_R = len(R_array)

        # Evaluate kF values at each point in R_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        # Calculate kF2 array if pn or np distribution
        if rho_2_array.any():
            kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
            pair = 'pn'
        # Otherwise set kF2=kF1 where previous functions will give pp or nn
        # distributions (depending on kF_1)
        else:
            kF2_array = kF1_array
            pair = 'pp'
            
        # Get each contribution with respect to q and Q
        n_I = self.n_I(q_array, Q_array, R_array, dR, kF1_array, kF2_array)
        n_deltaU = self.n_deltaU(q_array, Q_array, R_array, dR, kF1_array,
                                 kF2_array, pair)
        n_deltaU2 = self.n_deltaU2(q_array, Q_array, R_array, dR, kF1_array,
                                   kF2_array, pair)
        
        # Return total (ntot_q, ntot_Q)
        return n_I + n_deltaU + n_deltaU2
    
    
    def n_contributions(self, q_array, Q_array, R_array, dR, rho_1_array,
                        rho_2_array=np.empty(0)):
        """
        Contributions to the pair momentum distribution where the nucleonic
        densities specify the nucleus and distribution type (e.g., O16 and pn).
        This function isolates the I, \delta U, and \delta U \delta U^\dagger
        terms.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the first nucleon
            corresponding to \tau.
        rho_2_array : 1-D ndarray, optional
            Densities as a function of R [fm^-3] for the second nucleon
            corresponding to \tau'. If an empty array is input, the function
            assumes a proton-proton (or neutron-neutron) pair momentum
            distribution relying only on rho_1_array.

        Returns
        -------
        n_contributions : tuple
            Tuple of 2-D ndarrays corresponding to the I, \delta U, and
            \delta U^\dagger terms of the pair momentum distribution [fm^6] for
            each q and Q.

        """
        
        # Save lengths of q_array, Q_array, and R_array
        self.ntot_q = len(q_array)
        self.ntot_Q = len(Q_array)
        self.ntot_R = len(R_array)

        # Evaluate kF values at each point in R_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        # Calculate kF2 array if pn or np distribution
        if rho_2_array.any():
            kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
            pair = 'pn'
        # Otherwise set kF2=kF1 where previous functions will give pp or nn
        # distributions (depending on kF_1)
        else:
            kF2_array = kF1_array
            pair = 'pp'
            
        # Get each contribution with respect to q and Q
        n_I = self.n_I(q_array, Q_array, R_array, dR, kF1_array, kF2_array)
        n_deltaU = self.n_deltaU(q_array, Q_array, R_array, dR, kF1_array,
                                 kF2_array, pair)
        n_deltaU2 = self.n_deltaU2(q_array, Q_array, R_array, dR, kF1_array,
                                   kF2_array, pair)
        
        # Return tuple of contributions ( (ntot_q, ntot_Q), ... )
        return n_I, n_deltaU, n_deltaU2
    
    
    def write_file(self, nucleus, pair, Z, N, edf='SLY4'):
        """
        Write pair momentum distribution file for interpolation purposes.
        Split things into total, I, \delta U, and \delta U^2 contributions.

        Parameters
        ----------
        nucleus : str
            Specify nucleus (e.g., 'O16', 'Ca40', etc.)
        pair : str
            Specify 'pp', 'pn', or 'nn'.
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
        edf : str, optional
            Name of EDF (e.g., 'SLY4').
            
        Notes
        -----
        'pn' means \tau = +1/2 and \tau' = -1/2 and does NOT account for the
        opposite case (np). Many references combine these two. Multiply by 2 to
        match those references.

        """
        
        # Get relative momentum values
        q_array = self.k_array
        
        # Directory for distributions data
        data_directory = 'data/pmd/kvnn_%d/%s' % (self.kvnn, edf)
        
        # Create file name
        file_name = '%s_%s_channels' % (nucleus, pair)
        # Add each channel to file name
        for channel in self.channels:
            file_name += '_%s' % channel
        if self.generator == 'Block-diag':
            file_name += '_LambdaBD_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        else:
            file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        # Load R values and nucleonic densities (the R_array's are the same)
        R_array, rho_p_array = load_density(nucleus, 'proton', Z, N, edf)
        if edf == 'AV18' and Z == N: # e.g., AV18 He4 densities
            rho_n_array = rho_p_array
        else: # e.g., AV18 He8 densities
            R_array, rho_n_array = load_density(nucleus, 'neutron', Z, N, edf)
        dR = R_array[2] - R_array[1] # Assuming linear spacing
        
        # Set C.o.M. momentum values
        Q_max = 2.0 # Starts to get wonky at Q_max > 2.3 fm^-1
        ntot_Q = 40
        Q_array, Q_weights = gaussian_quadrature_mesh(Q_max, ntot_Q)

        # Calculate n_\lambda^\tau(q, Q) for each q in q_array and Q in Q_array
        if pair == 'pp':
            n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
                                    q_array, Q_array, R_array, dR, rho_p_array)
        elif pair == 'nn':
            n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
                                    q_array, Q_array, R_array, dR, rho_n_array)
        else: # pn pair
            n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
                       q_array, Q_array, R_array, dR, rho_p_array, rho_n_array)
            
        # Total momentum distribution
        n_total_array = n_I_array + n_delU_array + n_delU2_array
    
        # Open file and write header where we allocate roughly 18 centered
        # spaces for each label
        f = open(data_directory + '/' + file_name, 'w')
        header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}{:^18s}'.format('q',
                 'Q', 'total', '1', '\delta U', '\delta U^2')
        f.write(header + '\n')
    
        # Loop over momenta q and Q
        for iq, q in enumerate(q_array):
            for iQ, Q in enumerate(Q_array):

                # Write to data file following the format from the header
                line = '{:^18.6f}{:^18.6f}{:^18.6e}{:^18.6e}{:^18.6e}{:^18.6e}' \
                       .format( q, Q, n_total_array[iq, iQ], n_I_array[iq, iQ],
                                n_delU_array[iq, iQ], n_delU2_array[iq, iQ] )
                f.write('\n' + line)

        # Close file
        f.close()
        
        
    def n_lambda_interp(self, nucleus, pair, Z, N, edf='SLY4'):
        """
        Interpolate the pair momentum distribution for the specified file.

        Parameters
        ----------
        nucleus : str
            Specify nucleus (e.g., 'O16', 'Ca40', etc.)
        pair : str
            Specify 'pp', 'pn', or 'nn'.
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
        edf : str, optional
            Name of EDF (e.g., 'SLY4').
            
        Notes
        -----
        'pn' means \tau = +1/2 and \tau' = -1/2 and does NOT account for the
        opposite case (np). Many references combine these two. Multiply by 2 to
        match those references.

        """
        
        # Directory for distributions data
        data_directory = 'data/pmd/kvnn_%d/%s' % (self.kvnn, edf)
        
        # Get file name
        file_name = '%s_%s_channels' % (nucleus, pair)
        # Add each channel to file name
        for channel in self.channels:
            file_name += '_%s' % channel
        if self.generator == 'Block-diag':
            file_name += '_LambdaBD_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        else:
            file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        # Load data which includes all contributions to n_\lambda(q)
        data = np.loadtxt(data_directory + '/' + file_name)

        # Get C.o.M. and relative momentum values
        ntot_Q = 40
        ntot_q = self.ntot
        q_array = np.reshape( data[:, 0], (ntot_q, ntot_Q) )[:, 0]
        Q_array = np.reshape( data[:, 1], (ntot_q, ntot_Q) )[0, :]
        
        # Split data into 1-D arrays for each column
        # Total distribution
        n_total_array = np.reshape( data[:, 2], (ntot_q, ntot_Q) )
        # 1 term
        n_I_array = np.reshape( data[:, 3], (ntot_q, ntot_Q) )
        # \delta U term
        n_delU_array = np.reshape( data[:, 4], (ntot_q, ntot_Q) )
        # \delta U^2 term
        n_delU2_array = np.reshape( data[:, 5], (ntot_q, ntot_Q) )
        
        # Interpolate each array
        n_total_func = RectBivariateSpline(q_array, Q_array, n_total_array,
                                           kx=1, ky=1)
        n_I_func = RectBivariateSpline(q_array, Q_array, n_I_array, kx=1, ky=1)
        n_delU_func = RectBivariateSpline(q_array, Q_array, n_delU_array, kx=1,
                                          ky=1)
        n_delU2_func = RectBivariateSpline(q_array, Q_array, n_delU2_array,
                                           kx=1, ky=1)
        
        # Return all contributions with total first
        # Note, these are functions of q and Q
        return n_total_func, n_I_func, n_delU_func, n_delU2_func
    
    
    def n_Q0(self, q_array, R_array, dR, rho_1_array, rho_2_array=np.empty(0)):
        """
        Pair momentum distribution evaluated at Q = 0 where the nucleonic
        densities specify the nucleus and distribution type (e.g., O16 and pn).

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the first nucleon
            corresponding to \tau.
        rho_2_array : 1-D ndarray, optional
            Densities as a function of R [fm^-3] for the second nucleon
            corresponding to \tau'. If an empty array is input, the function
            assumes a proton-proton (or neutron-neutron) pair momentum
            distribution relying only on rho_1_array.

        Returns
        -------
        n_total : 1-D ndarray
            Pair momentum distribution [fm^6] for each q.

        """
        
        # Create Q_array with a single zero entry
        # \theta functions simplify to give 1 or 0 in this case
        Q_array = np.array([0.0])
        
        # Save lengths of q_array, Q_array, and R_array
        self.ntot_q = len(q_array)
        self.ntot_Q = len(Q_array)
        self.ntot_R = len(R_array)

        # Evaluate kF values at each point in R_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        # Calculate kF2 array if pn or np distribution
        if rho_2_array.any():
            kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
            distribution = 'pn'
        # Otherwise set kF2=kF1 where previous functions will give pp or nn
        # distributions (depending on kF_1)
        else:
            kF2_array = kF1_array
            distribution = 'pp'
            
        # Get each contribution with respect to q and Q
        n_I = self.n_I(q_array, Q_array, R_array, dR, kF1_array, kF2_array)
        n_deltaU = self.n_deltaU(q_array, Q_array, R_array, dR, kF1_array,
                                 kF2_array, distribution)
        n_deltaU2 = self.n_deltaU2(q_array, Q_array, R_array, dR, kF1_array,
                                   kF2_array, distribution)
        
        # Return total (ntot_q, 1)
        return ( n_I + n_deltaU + n_deltaU2 )[:, 0]