#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: dmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 30, 2021
# 
# Calculates SRG-evolved deuteron momentum distributions. (Note, 'dmd' stands
# for deuteron momentum distribution.)
#
# Notes on normalizations:
#   1. The deuteron wave function describing relative position or momentum is
#      normalized according to
#        \int dr r^2 ( |\psi_3S1(r)|^2 + |\psi_3D1(r)|^2 ) = 1,
#        2/\pi * \int dk k^2 ( |\psi_{3S1}(k)|^2 + |\psi_{3D1}(k)|^2 ) = 1.
#   2. Under HF+LDA, we adopt the normalization
#        4\pi / (2\pi)^3 \int dk k^2 < n_d^N(k) > = 1,
#        4\pi / (2\pi)^3 \int dk k^2 < n_d^{pn}(k) > = 1,
#      where angled-brackets indicate nuclear-averaging.
#
# Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
# screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!
#
# Revision history:
#   03/31/21 --- Moved hankel_transformation function to fourier_transform.py.
#   04/06/21 --- Added pair momentum distribution function called
#                'n_lambda_pair' and renamed original single-nucleon momentum
#                distribution from 'n_lambda' to 'n_lambda_single'. Also
#                merged pmd_deuteron_test.py with this script.
#   04/26/21 --- Correcting overall factors in front of \delta U and
#                \delta U^2 terms.
#   04/27/21 --- Fixed issue with evaluating at maximum q-point in q_array for
#                single-nucleon distribution.
#   05/26/21 --- Saving older version as dmd_v1.py with new updates: angle-
#                averaging, etc.
#   06/10/21 --- Including interpolating option to speed up code. Replacing
#                UnivariateSpline with interp1d for better accuracy, though
#                RectBivariateSpline works well for 2-D interpolations.
#   06/30/21 --- Speeding up code by switching from loops to np.sum() to do
#                integrations.
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d, RectBivariateSpline
# Scripts made by A.T.
from Figures import figures_functions as ff
from Misc.fourier_transform import hankel_transformation_k2r
from Misc.integration import gaussian_quadrature_mesh
import observables as ob
from operators import momentum_projection_operator
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class deuteron_momentum_distributions(object):
    
    
    def __init__(self, kvnn, lamb, kmax=0.0, kmid=0.0, ntot=0, interp=False):
        """
        Evaluates and saves 3S1-3D1 matrix elements of \delta U and
        \delta U^{\dagger} given the input potential and SRG \lambda.
        
        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        lamb : float
            SRG evolution parameter lambda [fm^-1].
        kmax : float, optional
            Maximum value in the momentum mesh [fm^-1]. (Default of zero
            automatically selects default mesh based on kvnn.)
        kmid : float, optional
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int, optional
            Number of momentum points in mesh.
        interp : bool, optional
            Option to use interpolated n_\lambda(q) functions.
            
        """
        
        # Get relevant info for file and directory names
        # Part of data directory name
        self.kvnn = kvnn 
        # Part of file name
        self.lamb = lamb
        self.kmax = kmax
        
        # Set-up for calculation
        if interp == False:
            
            # Channel is 3S1-3D1 for deuteron
            channel = '3S1'

            # Load and save momentum arrays
            k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid,
                                                   ntot)
            # Make sure you get actual size of momentum array
            if ntot == 0:
                ntot = len(k_array)
            # Save k_array, k_weights, ntot
            self.k_array, self.k_weights, self.ntot = k_array, k_weights, ntot
            
            # For dividing out momenta/weights
            factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
            # For coupled-channel matrices
            factor_array_cc = np.concatenate( (factor_array, factor_array) )

            # Load SRG transformation
            H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
            H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot,
                                             method='srg', generator='Wegner',
                                             lamb=lamb)
            # Load U(k, k') [unitless]
            U_matrix_unitless = SRG_unitary_transformation(H_initial, 
                                                           H_evolved)
            
            # Unitless wave function in momentum space
            psi_k_unitless = ob.wave_function(H_initial, U=U_matrix_unitless)
            self.psi = psi_k_unitless
            # Divide out momenta/weights and save
            self.psi_k = psi_k_unitless / factor_array_cc # [fm^3/2]

            # Isolate 2-body term and convert to fm^3
            I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
            row, col = np.meshgrid(factor_array_cc, factor_array_cc)
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            # Save for exact unitary calculation
            self.delta_U_matrix = delta_U_matrix_unitless
            delta_U_matrix = delta_U_matrix_unitless / row / col # fm^3
            
            # No 2*J+1 factor since we're fixing M_J for deuteron
            
            # Evaluate matrix elements
            deltaU = 1/2 * ( delta_U_matrix[:ntot, :ntot] + \
                             delta_U_matrix[ntot:, ntot:] )
            deltaU2 = 1/2 * ( delta_U_matrix[:ntot, :ntot]**2 + \
                              delta_U_matrix[:ntot, ntot:]**2 + \
                              delta_U_matrix[ntot:, :ntot]**2 + \
                              delta_U_matrix[ntot:, ntot:]**2 )

            # Interpolate < k | \delta U | k >
            self.deltaU_func = interp1d( k_array, np.diag(deltaU),
                                         kind='cubic', bounds_error=False,
                                         fill_value='extrapolate' )

            # Interpolate < k | \delta U \delta U^{\dagger} | k' >
            self.deltaU2_func = RectBivariateSpline(k_array, k_array, deltaU2)


    def theta_I(self, q_mesh, kF_mesh):
        """
        Evaluates \theta( kF(R) - q ). This function appears in the I term.

        Parameters
        ----------
        q_mesh : 4-D ndarray
            Momentum values [fm^-1].
        kF_mesh : 4-D ndarray
            Deuteron Fermi momentum [fm^-1] with respect to R.

        Returns
        -------
        theta_mesh : 2-D ndarray
            \theta function [unitless] evaluated for each q and kF(R).

        """
        
        # Initialize 2-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_R) )
        
        # Gives 1 if q < kF(R)
        theta_mesh[ q_mesh < kF_mesh ] = 1
        
        # This is a (ntot_q, ntot_R) size array
        return theta_mesh
        

    def theta_deltaU(self, q_mesh, kF_mesh, k_mesh):
        """
        Evaluates angle-average of \theta( kF(R) - \abs(q - 2k) ). This
        function appears in the \delta U term.

        Parameters
        ----------
        q_mesh : 3-D ndarray
            Momentum values [fm^-1].
        kF_mesh : 4-D ndarray
            Deuteron Fermi momentum [fm^-1] with respect to R.
        k_mesh : 3-D ndarray
            Relative momentum values [fm^-1].

        Returns
        -------
        theta_mesh : 3-D ndarray
            \theta function [unitless] evaluated for each q, kF(R), and k.
            
        Notes
        -----
        Not sure why the cases had to be computed in reverse order. Does that
        mean there is an overlap of truth values in case 3 and case 1 when
        there shouldn't be?

        """
        
        # Initialize 3-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_R, self.ntot_k) )
        
        # Evaluate each boolean case and use these to fill in the theta_mesh
        
        # This applies to each case: q < kF
        case_all = q_mesh < kF_mesh
        
        # Case 3: q-kF < 2k < kF+q
        case_3 = case_all * ( q_mesh - kF_mesh < 2*k_mesh ) * \
                 ( 2*k_mesh < kF_mesh + q_mesh )
        theta_mesh[case_3] = ( ( kF_mesh**2 - ( q_mesh - 2*k_mesh )**2 ) / \
                               ( 8*k_mesh*q_mesh ) )[case_3]
            
        # Case 2: kF-q <= 2k < kF+q
        case_2 = case_all * ( kF_mesh - q_mesh <= 2*k_mesh ) * \
                 ( 2*k_mesh < kF_mesh + q_mesh )
        theta_mesh[case_2] = ( ( kF_mesh**2 - ( q_mesh - 2*k_mesh )**2 ) / \
                               ( 8*k_mesh*q_mesh ) )[case_2]
        
        # Case 1: 2k < kF-q
        case_1 = case_all * ( 2*k_mesh < kF_mesh - q_mesh )
        theta_mesh[case_1] = 1

        # This is a (ntot_q, ntot_R, ntot_k) size array
        return theta_mesh


    def theta_deltaU2(self, kF_mesh, K_mesh, k_mesh):
        """
        Evaluates angle-average of \theta( kF(R) - \abs(K/2 + k) ) x
        \theta( kF(R) - \abs(K/2 - k) ). This function appears in the
        \delta U \delta U^\dagger term.

        Parameters
        ----------
        kF_mesh : 4-D ndarray
            Deuteron Fermi momentum [fm^-1] with respect to R.
        K_mesh : 4-D ndarray
            C.o.M. momentum values [fm^-1].
        k_mesh : 4-D ndarray
            Relative momentum values [fm^-1].

        Returns
        -------
        theta_mesh : 4-D ndarray
            \theta function [unitless] evaluated for each q, kF(R), K, and k.
            Note, this function does not depend on q but we need the array to
            match the size of \delta U \delta U^\dagger matrix elements.

        """
        
        # Initialize 4-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_R, self.ntot_K,
                                self.ntot_k) )

        # Evaluate each boolean case and use these to fill in the theta_mesh
        
        # Case 2: 2k+K > 2kF and 4k^2+K^2 < 4kF^2
        case_2 = ( 2*k_mesh + K_mesh > 2*kF_mesh ) * \
                 ( 4*k_mesh**2 + K_mesh**2 <= 4*kF_mesh**2 )    
        theta_mesh[case_2] = ( ( 4*kF_mesh**2 - 4*k_mesh**2 - K_mesh**2 ) / \
                               (4*k_mesh*K_mesh) )[case_2]

        # Case 1: 2k+K < 2kF
        case_1 = ( 2*k_mesh + K_mesh <= 2*kF_mesh )
        theta_mesh[case_1] = 1
        
        # This is a (ntot_q, ntot_R, ntot_K, ntot_k) size array
        return theta_mesh
    
    
    def n_I(self, q_array, R_array, dR, kF_array):
        """
        Evaluates the I term in U n(q) U^\dagger ~ I n(q) I.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF_array : 1-D ndarray
            Deuteron Fermi momentum [fm^-1] with respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from I term as a function of q.

        """
        
        # Initialize 2-D meshgrids (q, R)
        q_mesh, R_mesh = np.meshgrid(q_array, R_array, indexing='ij')
        
        # Get 2-D kF(R) mesh
        _, kF_mesh = np.meshgrid(q_array, kF_array, indexing='ij')
        
        # Evaluate the \theta-function in I term
        theta_mesh = self.theta_I(q_mesh, kF_mesh)
        
        # Calculate R integrand (ntot_q, ntot_R)
        integrand_R = theta_mesh * R_mesh**2 * dR

        # Integrate over R
        # This is a (ntot_q, 1) size array
        # Factor of 2 is overall factor
        return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    
    def n_deltaU(self, q_array, R_array, dR, kF_array):
        """
        Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
        we are combining \delta U and \delta U^\dagger, hence the factor of 2.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF_array : 1-D ndarray
            Deuteron Fermi momentum [fm^-1] with respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U term as a function of q.

        """
        
        # Set number of k points for integration over k
        self.ntot_k = 50
        
        # Initialize 3-D meshgrids (q, R, k) with k values equal to 0
        q_mesh, R_mesh, k_mesh = np.meshgrid(q_array, R_array,
                                             np.zeros(self.ntot_k),
                                             indexing='ij')
        
        # Get 3-D kF(R) mesh and initialize k weights mesh
        _, kF_mesh, dk_mesh = np.meshgrid(q_array, kF_array,
                                           np.zeros(self.ntot_k),
                                           indexing='ij')

        # Loop over q and R to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iq, q in enumerate(q_array):
            for iR, R in enumerate(R_array):
                
                kF = kF_array[iR]

                # Create integration mesh k_array up to (kF + q)/2 which
                # corresponds to the upper limit of \theta( kF2(R) - |q-2k| )
                k_max = (kF + q)/2
 
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max,
                                                              self.ntot_k)
                
                # Fill in k_mesh and dk_mesh given the specific k_array
                k_mesh[iq, iR, :] = k_array
                dk_mesh[iq, iR, :] = k_weights
        
        # Evaluate angle-average of \theta-functions in \delta U term
        theta_mesh = self.theta_deltaU(q_mesh, kF_mesh, k_mesh)
        
        # Evaluate < k | \delta U | k >
        deltaU_mesh = self.deltaU_func(k_mesh)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # Factor of 2 from \delta U + \delta U^\dagger
        # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        deltaU_factor = 2 * 8 * 2/np.pi * 2
        
        # Calculate the k integrand where we split terms according to pp and
        # pn (or nn and np if kF_1 corresponds to a neutron)
        # (ntot_q, ntot_R, ntot_k)
        integrand_k = deltaU_factor * k_mesh**2 * dk_mesh * R_mesh**2 * dR * \
                      deltaU_mesh * theta_mesh
        
        # Integrate over k leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R
        # This is a (ntot_q, 1) size array
        return 4*np.pi * np.sum(integrand_R, axis=-1)


    def n_deltaU2(self, q_array, R_array, dR, kF_array):
        """
        Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF_array : 1-D ndarray
            Deuteron Fermi momentum [fm^-1] with respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U \delta U^\dagger term as
            a function of q.

        """
        
        # Set number of K and k points for integration over K and k
        self.ntot_K = 50
        self.ntot_k = 50
        
        # y = cos(\theta) angles for averaging integration
        ntot_y = 7
        y_array, y_weights = leggauss(ntot_y)
        
        # Initialize 4-D meshgrids (q, R, K, k) with k and K equal to 0
        q_mesh, R_mesh, K_mesh, k_mesh = np.meshgrid(q_array, R_array,
                                                     np.zeros(self.ntot_K),
                                                     np.zeros(self.ntot_k),
                                                     indexing='ij')
        
        # Get 4-D kF(R) mesh and initialize K and k weights mesh
        _, kF_mesh, dK_mesh, dk_mesh = np.meshgrid(q_array, kF_array,
                                                    np.zeros(self.ntot_K),
                                                    np.zeros(self.ntot_k),
                                                    indexing='ij')

        # Loop over q, R, and k to find limits of K integration and then create
        # K_array using Gaussian quadrature
        for iR, R in enumerate(R_array):
                
            kF = kF_array[iR]
    
            # K integration goes from 0 to 2*kF
            K_max = 2*kF
                
            # Get Gaussian quadrature mesh for K integration
            K_array, K_weights = gaussian_quadrature_mesh(K_max, self.ntot_K)
              
            # Loop over remaining variables and fill in 5-D array
            for iq in range(self.ntot_q):
                for ik in range(self.ntot_k):
                    
                    # Fill in K_mesh and dK_mesh given the specific K_array
                    K_mesh[iq, iR, :, ik] = K_array
                    dK_mesh[iq, iR, :, ik] = K_weights

        # Loop over q, R, and K to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iR, R in enumerate(R_array):
        
            kF = kF_array[iR]
            
            # K_array only depends on R, so loop over K_mesh[0, iR, :, 0]
            for iK, K in enumerate( K_mesh[0, iR, :, 0] ):
                
                # Lower limit of k integration
                k_min = max(K/2 - kF, 0)
                
                # Upper limit of k integration
                if K**2/4 < kF**2:
                    k_max = min( np.sqrt(kF**2 - K**2/4), kF + K/2 )
                else:
                    k_max = kF + K/2

                # Get Gaussian quadrature mesh for k integration
                k_array, k_weights = gaussian_quadrature_mesh(k_max,
                                                              self.ntot_k,
                                                              xmin=k_min)
                
                # Loop over remaining variables and fill in 4-D array    
                for iq in range(self.ntot_q):
                    
                    # Fill in k_mesh and dk_mesh given the specific k_array
                    k_mesh[iq, iR, iK, :] = k_array
                    dk_mesh[iq, iR, iK, :] = k_weights
                        
        # Evaluate angle-average of \theta-functions in \delta U^2 term
        theta_mesh = self.theta_deltaU2(kF_mesh, K_mesh, k_mesh)

        # Evaluate 4-D < k | \delta U | |q_vec-K_vec/2| >^2 while also
        # averaging over angle y
        deltaU2_mesh = np.zeros_like(theta_mesh) # (q, R, K, k)
        
        # Integrate over y
        for y, dy in zip(y_array, y_weights):
            
            # Evaluate |q-K/2| mesh
            q_K_mesh = np.sqrt( q_mesh**2 + K_mesh**2/4 - q_mesh*K_mesh*y )
            
            deltaU2_mesh += self.deltaU2_func.ev(k_mesh, q_K_mesh) * dy/2
        
        # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2

        # Calculate the k integrand
        # (ntot_q, ntot_R, ntot_K, ntot_k)
        integrand_k = deltaU2_factor * k_mesh**2 * dk_mesh * K_mesh**2 * \
                      dK_mesh * R_mesh**2 * dR * deltaU2_mesh * theta_mesh

        # Integrate over k leaving K integrand (ntot_q, ntot_R, ntot_K)
        integrand_K = np.sum(integrand_k, axis=-1)
        
        # Integrate over K leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_K, axis=-1)
        
        # Integrate over R
        # This is a (ntot_q, 1) size array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
            
    
    def n_total(self, q_array, R_array, dR):
        """
        Deuteron momentum distribution under HF+LDA.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].

        Returns
        -------
        n_total : 1-D ndarray
            Deuteron momentum distribution [fm^3] for each q.
            
        Notes
        -----
        Should R_array be called r_array for relative coordinate?

        """
        
        # Save lengths of q_array and R_array
        self.ntot_q = len(q_array)
        self.ntot_R = len(R_array)
        
        # Get deuteron density
        k_array, k_weights, ntot_k = self.k_array, self.k_weights, self.ntot
        
        # Transform wave function to coordinate space
        hank_trans_3S1 = hankel_transformation_k2r('3S1', k_array, k_weights,
                                                   R_array)
        hank_trans_3D1 = hankel_transformation_k2r('3D1', k_array, k_weights,
                                                   R_array)
    
        # Get 3S1 and 3D1 waves [fm^-3/2] in coordinate-space
        psi_R_3S1 = hank_trans_3S1 @ self.psi_k[:ntot_k]
        psi_R_3D1 = hank_trans_3D1 @ self.psi_k[ntot_k:]
    
        # Calculate the deuteron nucleonic density [fm^-3]
        rho_array = psi_R_3S1**2 + psi_R_3D1**2

        # Evaluate kF values at each point in R_array
        kF_array = (3*np.pi**2 * rho_array)**(1/3)
            
        # Get each contribution with respect to q
        n_I = self.n_I(q_array, R_array, dR, kF_array)
        n_deltaU = self.n_deltaU(q_array, R_array, dR, kF_array)
        n_deltaU2 = self.n_deltaU2(q_array, R_array, dR, kF_array)

        # Return total (ntot_q, 1)
        return n_I + n_deltaU + n_deltaU2
    
    
    def n_contributions(self, q_array, R_array, dR):
        """
        Contributions to the deuteron momentum distribution. This function
        isolates the I, \delta U, and \delta U \delta U^\dagger terms.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].

        Returns
        -------
        n_contributions : tuple
            Tuple of 1-D ndarrays corresponding to the I, \delta U, and
            \delta U^\dagger terms of the pair momentum distribution [fm^3] for
            each q.

        """
        
        # Save lengths of q_array and R_array
        self.ntot_q = len(q_array)
        self.ntot_R = len(R_array)
        
        # Get deuteron density
        k_array, k_weights, ntot_k = self.k_array, self.k_weights, self.ntot
        
        # Transform wave function to coordinate space
        hank_trans_3S1 = hankel_transformation_k2r('3S1', k_array, k_weights,
                                                   R_array)
        hank_trans_3D1 = hankel_transformation_k2r('3D1', k_array, k_weights,
                                                   R_array)
    
        # Get 3S1 and 3D1 waves [fm^-3/2] in coordinate-space
        psi_R_3S1 = hank_trans_3S1 @ self.psi_k[:ntot_k]
        psi_R_3D1 = hank_trans_3D1 @ self.psi_k[ntot_k:]
    
        # Calculate the deuteron nucleonic density [fm^-3]
        rho_array = psi_R_3S1**2 + psi_R_3D1**2

        # Evaluate kF values at each point in R_array
        kF_array = (3*np.pi**2 * rho_array)**(1/3)
            
        # Get each contribution with respect to q
        n_I = self.n_I(q_array, R_array, dR, kF_array)
        n_deltaU = self.n_deltaU(q_array, R_array, dR, kF_array)
        n_deltaU2 = self.n_deltaU2(q_array, R_array, dR, kF_array)
        
        # Return tuple of contributions ( (ntot_q, ntot_Q), ... )
        return n_I, n_deltaU, n_deltaU2
    
    
    def write_file(self):
        """
        Write deuteron momentum distribution file for interpolation purposes.
        Split things into total, I, \delta U, and \delta U^2 contributions.

        """
        
        # Get momentum values
        q_array = self.k_array
        
        # Directory for distributions data
        data_directory = 'Data/dmd/kvnn_%d' % self.kvnn
        
        # Create file name
        file_name = 'deuteron_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        
        # Use R values in accordance with densities code from densities.py
        R_array = np.linspace(0.1, 20.0, 200)
        dR = R_array[1] - R_array[0] # Assuming linear spacing

        # Calculate n_\lambda^d(k) for each k in k_array
        n_I_array, n_delU_array, n_delU2_array = self.n_contributions(q_array,
                                                                      R_array,
                                                                      dR)
            
        # Total momentum distribution
        n_total_array = n_I_array + n_delU_array + n_delU2_array
    
        # Open file and write header where we allocate roughly 18 centered
        # spaces for each label
        f = open(data_directory + '/' + file_name, 'w')
        header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}'.format('q',
                  'total', '1', '\delta U', '\delta U^2')
        f.write(header + '\n')
    
        # Loop over momenta k
        for ik, k in enumerate(q_array):

            # Write to data file following the format from the header
            line = '{:^18.6f}{:^18.6e}{:^18.6e}{:^18.6e}{:^18.6e}'.format( k,
                            n_total_array[ik], n_I_array[ik], n_delU_array[ik],
                            n_delU2_array[ik] )
            f.write('\n' + line)

        # Close file
        f.close()
        
        
    def n_lambda_interp(self):
        """
        Interpolate the deuteron momentum distribution for the specified file.
    
        Returns
        -------
        output : tuple
            Tuple of functions that depend only on momentum q [fm^-1] where
            each function corresponds to contributions to n_\lambda(q): total,
            I, \delta U, and \delta U^2.

        """
        
        # Directory for distributions data
        data_directory = 'Data/dmd/kvnn_%d' % self.kvnn
        
        # Get file name
        file_name = 'deuteron_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        # Load data which includes all contributions to n_\lambda(q)
        data = np.loadtxt(data_directory + '/' + file_name)
        
        # Split data into 1-D arrays for each column
        q_array = data[:, 0] # Momentum in fm^-1
        n_total_array = data[:, 1] # Total distribution
        n_1_array = data[:, 2] # 1 term
        n_delU_array = data[:, 3] # \delta U term
        n_delU2_array = data[:, 4] # \delta U^2 term
        
        # Interpolate each array (UnivariateSpline is for smoothing whereas
        # interp1d gives closer value to the actual calculation)
        n_total_func = interp1d(q_array, n_total_array, bounds_error=False,
                                kind='cubic', fill_value='extrapolate')
        n_1_func = interp1d(q_array, n_1_array, bounds_error=False,
                            kind='cubic', fill_value='extrapolate')
        n_delU_func = interp1d(q_array, n_delU_array, bounds_error=False,
                               kind='cubic', fill_value='extrapolate')
        n_delU2_func = interp1d(q_array, n_delU2_array, bounds_error=False,
                                kind='cubic', fill_value='extrapolate')
        
        # Return all contributions with total first
        # Note, these are functions of q
        return n_total_func, n_1_func, n_delU_func, n_delU2_func
    
    
    def n_lambda_exact(self, q, contributions='total'):
        """
        Computes the exact deuteron momentum distribution using the wave
        function from the input Hamiltonian. Note, set interp = False to use
        this function.

        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        contributions : str, optional
            Option to return different contributions to the momentum
            distribution.
            1. Default is 'total' which only returns the total momentum
               distribution.
            2. Specify 'q_contributions' for I, \delta U, and 
               \delta U \delta U^\dagger along with the total.
            3. Specify 'partial_wave_ratio' for the ratio of 3S1-3S1 to full
               3S1-3D1 along with the absolute total.

        Returns
        -------
        output : float or tuple
            Pair momentum distribution [fm^3] evaluated at momentum q
            [fm^-1]. (Note, this function will return a tuple of floats if
            contributions is not 'total'.)

        """

        # Load deuteron wave function and delta_U_matrix (both unitless)
        psi_vector = self.psi
        delta_U_matrix = self.delta_U_matrix
        
        # Load bare momentum projection operator [fm^3]
        bare_operator = momentum_projection_operator(q, self.k_array,
                        self.k_weights, '3S1', smeared=False)
        
        # Term 1: 1 * n(q) * 1
        term_1 = psi_vector.T @ bare_operator @ psi_vector
        
        # Term 2: \deltaU * n(q) * 1
        term_deltaU = psi_vector.T @ delta_U_matrix @ bare_operator @ \
                      psi_vector
        
        # Term 3: 1 * n(q) * \deltaU^\dagger = Term 2
        term_deltaU *= 2
        
        # High-q term: \deltaU * n(q) * \deltaU^\dagger
        term_deltaU2 = psi_vector.T @ delta_U_matrix @ bare_operator @ \
                        delta_U_matrix.T @ psi_vector
                       
        if contributions == 'partial_wave_ratio':
            
            # Length of q_array
            ntot = self.ntot
            
            # 3S1-3S1 only
            numerator = abs( psi_vector.T[:ntot] @ ( \
                                 delta_U_matrix[:ntot, :ntot] @ \
                                 bare_operator[:ntot, :ntot] @ \
                                 delta_U_matrix.T[:ntot, :ntot] +
                                 delta_U_matrix[:ntot, ntot:] @ \
                                 bare_operator[ntot:, ntot:] @ \
                                 delta_U_matrix.T[ntot:, :ntot] ) @ \
                             psi_vector[:ntot] )
                        
            # Full 3S1-3D1 taking absolute values 
            denominator = numerator + \
                          abs( psi_vector.T[:ntot] @ ( \
                                   delta_U_matrix[:ntot, :ntot] @ \
                                   bare_operator[:ntot, :ntot] @ \
                                   delta_U_matrix.T[:ntot, ntot:] +
                                   delta_U_matrix[:ntot, ntot:] @ \
                                   bare_operator[ntot:, ntot:] @ \
                                   delta_U_matrix.T[ntot:, ntot:] ) @ \
                               psi_vector[ntot:] ) + \
                          abs( psi_vector.T[ntot:] @ ( \
                                   delta_U_matrix[ntot:, :ntot] @ \
                                   bare_operator[:ntot, :ntot] @ \
                                   delta_U_matrix.T[:ntot, :ntot] +
                                   delta_U_matrix[ntot:, ntot:] @ \
                                   bare_operator[ntot:, ntot:] @ \
                                   delta_U_matrix.T[ntot:, :ntot] ) @ \
                               psi_vector[:ntot] ) + \
                          abs( psi_vector.T[ntot:] @ ( \
                                   delta_U_matrix[ntot:, :ntot] @ \
                                   bare_operator[:ntot, :ntot] @ \
                                   delta_U_matrix.T[:ntot, ntot:] +
                                   delta_U_matrix[ntot:, ntot:] @ \
                                   bare_operator[ntot:, ntot:] @ \
                                   delta_U_matrix.T[ntot:, ntot:] ) @ \
                               psi_vector[ntot:] )
                                   
        # Add up each term for total
        total = term_1 + term_deltaU + term_deltaU2

        # Return contributions and total or just total
        if contributions == 'q_contributions':
            return total, term_1, term_deltaU, term_deltaU2
        elif contributions == 'partial_wave_ratio':
            return total, numerator/denominator
        else: # Default
            return total