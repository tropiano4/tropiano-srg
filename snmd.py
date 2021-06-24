#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: snmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 25, 2021
# 
# Calculates SRG-evolved single-nucleon momentum distributions for nuclei
# assuming the evolved wave function is given by a free Fermi gas. (Note,
# 'snmd' stands for single-nucleon momentum distribution.)
#
# Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
# screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!
#
# Revision history:
#   03/29/21 --- Finalized and produces same results (up to factors of 2*\pi)
#                as single_particle_momentum_dist.py. Also, fixed a bug
#                involving the integration over total momentum K.
#   03/31/21 --- Moved channel_L_value function to vnn.py. See vnn.py for
#                details of the function.
#   04/02/21 --- Added option to return pp and pn (or nn and np) contributions
#                to \delta U \delta U^\dagger term along with total.
#   04/22/21 --- Correcting overall factors in front of \delta U and
#                \delta U^2 terms.
#   04/27/21 --- Fixed issue with evaluating at maximum q-point in q_array.
#   05/25/21 --- Saving older version as snmd_v1.py with new updates: angle-
#                averaging, etc.
#   06/04/21 --- Including interpolating option to speed up code.
#   06/10/21 --- Replacing UnivariateSpline with interp1d for better accuracy,
#                though RectBivariateSpline works well for 2-D interpolations.
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d, RectBivariateSpline
# Scripts made by A.T.
from densities import load_density
from Figures import figures_functions as ff
from Misc.integration import gaussian_quadrature_mesh
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class single_nucleon_momentum_distributions(object):
    
    
    def __init__(self, kvnn, channels, lamb, kmax=0.0, kmid=0.0, ntot=0,
                 interp=True):
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
            Maximum value in the momentum mesh [fm^-1]. (Default of zero
            automatically selects default mesh based on kvnn.)
        kmid : float, optional
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int, optional
            Number of momentum points in mesh.
        interp : bool, optional
            Option to use interpolated n_\lambda(q) functions.
            
        """
        
        # Part of data directory name
        self.kvnn = kvnn 
        # Part of file name
        self.channels = channels
        self.lamb = lamb
        self.kmax = kmax
        
        # Calculate directly and/or generate data files
        if interp == False:
            
            # Save highest allowed L based on input channels
            highest_L = 0
            for channel in channels:
                next_L = vnn.channel_L_value(channel)
                if next_L > highest_L:
                    highest_L = next_L

            # Relative momentum k [fm^-1] (channel doesn't matter here)
            k_array, k_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid,
                                                   ntot)
            self.k_array = k_array # Save for generating files
            if ntot == 0:
                ntot = len(k_array) # Make sure ntot is the number of k-points
        
            # cos(\theta) angles for averaging integration
            xtot = 9
            x_array, x_weights = leggauss(xtot)
            self.x_array, self.x_weights, self.xtot = x_array, x_weights, xtot
        
            # For dividing out momenta/weights
            factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
            # For coupled-channel matrices
            factor_array_cc = np.concatenate( (factor_array, factor_array) )
            
            
            # --- Evaluate matrix elements for each channel --- #
        
            # Initialize pp and pn matrix elements
            deltaU_pp = np.zeros( (ntot, ntot) ) # \delta U linear term
            deltaU_pn = np.zeros( (ntot, ntot) )
            deltaU2_pp = np.zeros( (ntot, ntot) ) # \delta U \delta U^\dagger
            deltaU2_pn = np.zeros( (ntot, ntot) )
        
            # Allowed channels for pp (and nn) up through the D-waves
            pp_channels = ('1S0', '3P0', '3P1', '3P2', '1D2')
        
            # Loop over channels and evaluate matrix elements
            for channel in channels:

                # Load SRG transformation
                H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid,
                                                 ntot)
                H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid,
                                                 ntot, method='srg',
                                                 generator='Wegner', lamb=lamb)
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
                    deltaU2_pn += (2*J+1)/2 * (
                                      delta_U_matrix[:ntot, :ntot]**2 + \
                                      delta_U_matrix[:ntot, ntot:]**2 )

                    # Isospin CG's=1 for pp
                    if channel in pp_channels:
                        deltaU_pp += (2*J+1) * delta_U_matrix[:ntot, :ntot]
                        deltaU2_pp += (2*J+1) * (
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

            # Interpolate pp and pn < k | \delta U | k > matrix elements
            self.deltaU_pp_func = interp1d( k_array, np.diag(deltaU_pp),
                                            bounds_error=False,
                                            fill_value='extrapolate' )
            self.deltaU_pn_func = interp1d( k_array, np.diag(deltaU_pn),
                                            bounds_error=False,
                                            fill_value='extrapolate' )
        
            # Interpolate pp and pn < k | \delta U \delta U^{\dagger} | k' > 
            # matrix elements
            self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                       deltaU2_pp)
            self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                       deltaU2_pn)
            
            # print(deltaU_pp)
            # print(deltaU_pn)
            # print(deltaU2_pp)
            # print(deltaU2_pn)
            print(self.deltaU_pp_func(2.0))
            print(self.deltaU2_pp_func.ev(2.0, 2.0))
            print(self.deltaU_pn_func(2.0))
            print(self.deltaU2_pn_func.ev(2.0, 2.0))
            

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


    def theta_deltaU(self, kF, q, k_array, ntot_k):
        """
        Evaluates \theta( k_F - \abs(q_vec - 2k_vec) ) for every momentum k.
        This function appears in the \delta U term.

        Parameters
        ----------
        kF : float
            Fermi momentum [fm^-1].
        q : float
            Single-nucleon momentum [fm^-1].
        k_array : 1-D ndarray
            Momentum values [fm^-1] of \int dk k^2 integral in \delta U term.
        ntot_k : int
            Number of momentum points in momentum mesh k_array.

        Returns
        -------
        theta_deltaU : 1-D ndarray
            \theta function evaluated at each momenta k. This is a unitless
            array of length ntot where ntot = len(k_array).
            
        Notes
        -----
        This function assumes that the k_array does not exceed (kF + q)/2
        which is the maximum limit given by the \theta function, and that
        q < kF_1 (with the kF in this function corresponding to kF_2).

        """
        
        # Make \theta( k_F - \abs(q_vec - 2k_vec) ) the same length as k_array
        theta_deltaU = np.zeros(ntot_k)
        
        # Loop over each momenta k and go through the three inequalities
        for i, k in enumerate(k_array):
            
            # Case 1: q < kF and 2k < kF-q -> h(q,k) = 1
            if q < kF and 2*k <= kF-q:
                theta_deltaU[i] = 1
            
            # Case 2: q < kF and kF-q < 2k < kF+q
            # -> h(q,k) = ( kF^2 - (q-2k)^2 ) / (8*q*k)
            elif q < kF and kF-q < 2*k < kF+q:
                theta_deltaU[i] = ( kF**2 - ( q - 2*k )**2 ) / ( 8*k*q )
                
            # Case 3: q-kF < 2k < kF+q
            # -> h(q,k) = ( kF^2 - (q-2k)^2 ) / (8*q*k)
            elif q-kF < 2*k < kF+q:
                theta_deltaU[i] = ( kF**2 - ( q - 2*k )**2 ) / ( 8*k*q )
                
            # Otherwise, h(q,k) = 0
                
        return theta_deltaU


    def theta_deltaU2(self, kF_1, kF_2, K, k_array, ntot_k):
        """
        Evaluates \theta( kF_1 - \abs(K_vec/2 + k_vec) ) \times
        \theta( kF_2 - \abs(K_vec/2 - k_vec) ) for every momentum k.
        This function appears in the \delta U \delta U^\dagger term.

        Parameters
        ----------
        kF_1 : float
            Fermi momentum [fm^-1] for the corresponding nucleon momentum
            distribution.
        kF_2 : float
            Fermi momentum [fm^-1] for the correlated nucleon. If kF_1
            corresponds to a proton, then kF_2 corresponds to a neutron (and
            vice versa).
        K : float
            C.o.M. momentum [fm^-1].
        k_array : 1-D ndarray
            Momentum values [fm^-1] of \int dk k^2 integral in \delta U term.
        ntot_k : int
            Number of momentum points in momentum mesh k_array.

        Returns
        -------
        theta_deltaU2 : 1-D ndarray
            \theta function evaluated at each momenta k. This is a unitless
            array of length ntot where ntot = len(k_array).
            
        Notes
        -----
        This function assumes that the k_array is within the interval
        [k_min, k_max] where
            k_min = max(K/2-kF_min, 0),
            k_max = min( \sqrt( kF_avg^2 - K^2/4 ), kF_min + K/2 ),
            kF_avg = \sqrt( (kF_1^2 + kF_2^2)/2 ),
        which are the limits given by the \theta functions.

        """
    
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
    
    
    def n_1(self, q, kF_1):
        """
        Evaluates first term in U n(q) U^\dagger ~ I n(q) I which gives
        2 \theta( kF(r) - q ).
        
        Parameters
        ----------
        q : float
            Single-nucleon momentum [fm^-1].
        kF_1 : float
            Fermi momentum [fm^-1] for the corresponding nucleon momentum
            distribution.
            
        Returns
        -------
        output : float
            Momentum distribution from I term before integration over
            \int dr r^2.
        
        """
        
        # This is simply a \theta function where the factor of two comes from
        # summing over particle spin projection \sigma
        if q < kF_1:
            
            return 2
        
        else:
            
            return 0
        
        
    def n_deltaU(self, q, kF_1, kF_2):
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
        kF_2 : float
            Fermi momentum [fm^-1] for the correlated nucleon. If kF_1
            corresponds to a proton, then kF_2 corresponds to a neutron (and
            vice versa).
            
        Returns
        -------
        output : float
            Momentum distribution from \delta U term before integration over
            \int dr r^2.
            
        """

        # Initial \theta( kF_1(r) - q )
        if q < kF_1:

            # Create integration mesh k_array up to (kF_2 + q)/2 where
            # (kF_2 + q)/2 corresponds to the upper limit of
            # \theta(kF_2(r)-|q_vec-2k_vec|)
            kmax_delU = (kF_2 + q)/2

            # Select number of integration points based on kmax_delU
            ntot_delU = self.select_number_integration_points(kmax_delU)

            # Get Gaussian quadrature mesh
            k_array_delU, k_weights_delU = gaussian_quadrature_mesh(kmax_delU,
                                                                    ntot_delU)
            
            # Average \theta( kF_2 - |q_vec-2k_vec| ) for kF_1 and kF_2
            # (meaning \tau=\tau' and \tau=-\tau')
            theta_pp_array = self.theta_deltaU(kF_1, q, k_array_delU,
                                               ntot_delU)
            theta_pn_array = self.theta_deltaU(kF_2, q, k_array_delU,
                                               ntot_delU)
            
            # Evaluate pp and pn < k | \delta U | k > matrix elements (or nn
            # and np if kF_1 corresponds to a neutron)
            deltaU_pp_k = self.deltaU_pp_func(k_array_delU)
            deltaU_pn_k = self.deltaU_pn_func(k_array_delU)
            
            # Build integrand for k integration where we split terms according
            # to pp and pn (or nn and np if kF_1 corresponds to a neutron)
            integrand_k = k_array_delU**2 * k_weights_delU * ( \
                              deltaU_pp_k * theta_pp_array + \
                              deltaU_pn_k * theta_pn_array )
                
            # Factor of 2 from \delta U + \delta U^\dagger
            # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
            # 2/\pi for two | k_vec > -> | k J L S ... > changes
            # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
            # combine to give 2
            deltaU_factor = 2 * 8 * 2/np.pi * 2
            
            # Integrate over \int dk k^2
            return deltaU_factor * np.sum(integrand_k)
        
        else:
            
            return 0


    def n_deltaU2_K(self, q, K, kF_1, kF_2):
        """
        Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger
        up to integration over \int dK K^2.
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        K : float
            C.o.M. momentum [fm^-1].
        kF_1 : float
            Fermi momentum [fm^-1] for the corresponding nucleon momentum
            distribution.
        kF_2 : float
            Fermi momentum [fm^-1] for the correlated nucleon. If kF_1
            corresponds to a proton, then kF_2 corresponds to a neutron (and
            vice versa).

        Returns
        -------
        output : float
            Momentum distribution from \delta U \delta U^\dagger term before
            integration over \int dr r^2 \int dK K^2.
            
        """
        
        # Create integration mesh k_array from minimum and maximum k values
        # corresponding to the limits of \theta(kF(r)-|K_vec/2 +(-) k_vec|)
        kF_min = min(kF_1, kF_2)
        kmin_delU2 = max(K/2 - kF_min, 0)
        kmax_delU2 = min( np.sqrt( ( kF_1**2 + kF_2**2 )/2 - K**2/4 ),
                          kF_min + K/2 )

        # Select number of integration points based on kmax_delU2
        ntot_delU2 = self.select_number_integration_points( kmax_delU2,
                                                            kmin_delU2 )

        # Get Gaussian quadrature mesh
        k_array_delU2, k_weights_delU2 = gaussian_quadrature_mesh(kmax_delU2,
                                         ntot_delU2, xmin=kmin_delU2)
        
        # Create a grid for evaluation of < k | \delta U | |q_vec-K_vec/2| >^2
        k_grid, x_grid = np.meshgrid(k_array_delU2, self.x_array,
                                     indexing='ij')
        # Get angle weights too
        _, x_weights_grid = np.meshgrid(k_array_delU2, self.x_weights,
                                        indexing='ij')
        
        # Create grid for |q_vec-K_vec/2| values
        q_K_grid = np.sqrt( q**2 + K**2/4 - q*K*x_grid )
        
        # Evaluate 2-D < k | \delta U | |q_vec-K_vec/2| >^2
        # These are (ntot_delU2, xtot) arrays
        deltaU2_pp_k_x = self.deltaU2_pp_func.ev(k_grid, q_K_grid) * \
                         x_weights_grid
        deltaU2_pn_k_x = self.deltaU2_pn_func.ev(k_grid, q_K_grid) * \
                         x_weights_grid
        
        # Integrate over \int dx/2 -> (ntot_k, 1)
        deltaU2_pp_k = np.sum(deltaU2_pp_k_x, axis=-1)/2
        deltaU2_pn_k = np.sum(deltaU2_pn_k_x, axis=-1)/2
        
        # Average \theta(kF(r)-|K_vec/2 +(-) k_vec|) functions
        # This is a (ntot_delU2, 1) array
        theta_array = self.theta_deltaU2(kF_1, kF_2, K, k_array_delU2,
                                         ntot_delU2)
        
        # Build integrand for k integration where we split terms according
        # to pp and pn (or nn and np if kF_1 corresponds to a neutron)
        integrand_k = k_array_delU2**2 * k_weights_delU2 * ( \
                      deltaU2_pp_k + deltaU2_pn_k ) * theta_array

        # Integrate over \int dk k^2 leaving K-dependent part and keeping the
        # overall factor in the n_deltaU2 function below
        return np.sum(integrand_k)


    def n_deltaU2(self, q, kF_1, kF_2):
        """
        Evaluates integration over \int dK K^2 in fourth term in
        U n(q) U^\dagger ~ \delta U \delta U^\dagger.
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        kF_1 : float
            Fermi momentum [fm^-1] for the corresponding nucleon momentum
            distribution.
        kF_2 : float
            Fermi momentum [fm^-1] for the correlated nucleon. If kF_1
            corresponds to a proton, then kF_2 corresponds to a neutron (and
            vice versa).
            
        Returns
        -------
        output : float
            Momentum distribution from \delta U \delta U^\dagger term before
            integration over \int dr r^2.
            
        """
        
        # Create K_mesh from 0 - (kF_1 + kF_2)
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
        
        # Calculate K-integrand which is a (Ntot, 1) array
        integrand_K = np.zeros(Ntot)
        
        # Loop over K
        for iK, K in enumerate(K_array):
            
            integrand_K[iK] = self.n_deltaU2_K(q, K, kF_1, kF_2)
            
        # Attach integration measure
        integrand_K *= K_array**2 * K_weights
        
        # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2
        
        # Integrate over \int dK K^2
        return deltaU2_factor * np.sum(integrand_K)
        
        
    def n_lambda(self, q_array, r_array, rho_1_array, rho_2_array):
        """
        Single-nucleon momentum distribution where the nucleonic densities
        specify the nucleus.
        
        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        r_array : 1-D ndarray
            Coordinates array [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of r [fm^-3] for the corresponding nucleon
            momentum distribution. For instance, if this array corresponds to
            proton densities, then this function will give a proton momentum
            distribution.
        rho_2_array : 1-D ndarray
            Densities as a function of r [fm^-3] for the correlated nucleon. If
            rho_1_array corresponds to a proton, then rho_2_array corresponds
            to a neutron (and vice versa).
            
        Returns
        -------
        n_lambda : 2-D ndarray
            Array of contributions to the single-nucleon momentum distribution
            ordered according to total, 1, \delta U, and \delta U^2 [fm^3].
            
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

        # Evaluate kF values at each point in r_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
            
        # Loop over q
        for i, q in enumerate(q_array):
    
            # n_lambda contributions before integrating over r
            n_lambda_r = np.zeros( (mtot, axes) )
            
            # Loop over kF values
            for j, (kF_1, kF_2) in enumerate( zip(kF1_array, kF2_array) ):
                
                term_1 = self.n_1(q, kF_1)
                term_deltaU = self.n_deltaU(q, kF_1, kF_2)
                term_deltaU2 = self.n_deltaU2(q, kF_1, kF_2)
                total = term_1 + term_deltaU + term_deltaU2
                
                n_lambda_r[j, :] = total, term_1, term_deltaU, term_deltaU2

            # Integrate over 4\pi \int dr r^2 for each contribution (summing
            # over axis=0)
            n_lambda[i, :] = 4*np.pi*dr * np.sum(r2_grid * n_lambda_r, axis=0)
        
        # Return (ntot, 4) array of contributions to n_\lambda^\tau(q)
        return n_lambda
    
    
    def write_files(self, nucleus, nucleon, Z, N, edf='SLY4'):
        """
        Write single-nucleon momentum distribution files for interpolation
        purposes. Split things into total, 1, \delta U, and \delta U^2
        contributions.

        Parameters
        ----------
        nucleus : str
            Specify nucleus (e.g., 'O16', 'Ca40', etc.)
        nucleon : str
            Specify 'proton' or 'neutron'.
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
        edf : str, optional
            Name of EDF (e.g., 'SLY4').

        """
        
        # Directory for distributions data
        data_directory = 'Data/snmd/kvnn_%d' % self.kvnn
        
        # Create file name
        file_name = '%s_%s_channels' % (nucleus, nucleon)
        # Add each channel to file name
        for channel in self.channels:
            file_name += '_%s' % channel
        file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        # Load r values and nucleonic densities (the r_array's are the same)
        r_array, rho_p_array = load_density(nucleus, 'proton', Z, N, edf)
        if edf == 'AV18':
            rho_n_array = rho_p_array
        else:
            r_array, rho_n_array = load_density(nucleus, 'neutron', Z, N, edf)

        # Calculate n_\lambda^\tau(k) for each k in k_array
        if nucleon == 'proton':
            n_array = self.n_lambda(self.k_array, r_array, rho_p_array,
                                    rho_n_array)
        elif nucleon == 'neutron':
            n_array = self.n_lambda(self.k_array, r_array, rho_n_array,
                                    rho_p_array)
    
        # Open file and write header where we allocate roughly 18 centered
        # spaces for each label
        f = open(data_directory + '/' + file_name, 'w')
        header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}'.format('q',
                 'total', '1', '\delta U', '\delta U^2')
        f.write(header + '\n')
    
        # Loop over momenta k
        for ik, k in enumerate(self.k_array):

            # Write to data file following the format from the header
            line = '{:^18.6f}{:^18.6e}{:^18.6e}{:^18.6e}{:^18.6e}'.format( k,
                       n_array[ik, 0], n_array[ik, 1], n_array[ik, 2],
                       n_array[ik, 3] )
            f.write('\n' + line)

        # Close file
        f.close()
        
        
    def n_lambda_interp(self, nucleus, nucleon, Z, N):
        """
        Interpolate the single-nucleon momentum distribution for the specified
        file.

        Parameters
        ----------
        nucleus : str
            Specify nucleus (e.g., 'O16', 'Ca40', etc.)
        nucleon : str
            Specify 'proton' or 'neutron'.
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
            
        Returns
        -------
        output : tuple
            Tuple of functions that depend only on momentum q [fm^-1] where
            each function corresponds to contributions to n_\lambda(q): total,
            1, \delta U, and \delta U^2.

        """
        
        # Directory for distributions data
        data_directory = 'Data/snmd/kvnn_%d' % self.kvnn
        
        # Get file name
        file_name = '%s_%s_channels' % (nucleus, nucleon)
        # Add each channel to file name
        for channel in self.channels:
            file_name += '_%s' % channel
        file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
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
                                fill_value='extrapolate')
        n_1_func = interp1d(q_array, n_1_array, bounds_error=False,
                            fill_value='extrapolate')
        n_delU_func = interp1d(q_array, n_delU_array, bounds_error=False,
                               fill_value='extrapolate')
        n_delU2_func = interp1d(q_array, n_delU2_array, bounds_error=False,
                                fill_value='extrapolate')
        
        # Return all contributions with total first
        # Note, these are functions of q
        return n_total_func, n_1_func, n_delU_func, n_delU2_func


if __name__ == '__main__':
    

    # Generate all data for single-nucleon momentum distributions
    # Currently this takes about ~60 hours to run
    import time
    
    # Potentials
    kvnns_list = [6, 222, 224]
    
    # Channels to include in calculation (S-waves only or higher partial waves)
    channels_list = [ ('1S0', '3S1'), ('1S0', '3S1', '3P0', '1P1', '3P1') ]
    
    # SRG \lambda values
    lambdas_list = [6.0, 3.0, 2.0, 1.5, 1.35]
    
    # Momentum mesh details
    kmax, kmid, ntot = 15.0, 3.0, 120 # Default

    # Nuclei to calculate
    nuclei_list = [ ('He4', 2, 2), ('C12', 6, 6), ('O16', 8, 8),
                    ('Ca40', 20, 20), ('Ca48', 20, 28), ('Fe56', 26, 30),
                    ('Pb208', 82, 126) ]
    
    # Loop over every case generating data files
    for kvnn in kvnns_list:
        
        t0_k = time.time()
        
        for ic, channels in enumerate(channels_list):
            
            t0_c = time.time()
            
            for lamb in lambdas_list:
                
                t0_l = time.time()
                
                for nuclei in nuclei_list:
                    
                    t0_N = time.time()
                    
                    # Initialize class
                    snmd = single_nucleon_momentum_distributions(kvnn,
                           channels, lamb, kmax, kmid, ntot, interp=False)

                    nucleus = nuclei[0]
                    if nucleus == 'He4':
                        edf = 'AV18'
                    else:
                        edf = 'SLY4'
                    Z = nuclei[1]
                    N = nuclei[2]
                    
                    # Write proton and neutron file for given nucleus
                    snmd.write_files(nucleus, 'proton', Z, N, edf)
                    snmd.write_files(nucleus, 'neutron', Z, N, edf)
                    
                    # Time for each nucleus to run
                    t1_N = time.time()
                    mins_N = (t1_N-t0_N)/60
                    print('\t\t\tDone with %s after %.5f minutes.' % (nucleus,
                                                                mins_N) )
            
                # Time for each \lambda to run   
                t1_l = time.time()
                mins_l = (t1_l-t0_l)/60
                print( '\n\t\tDone with \lambda=%s after %.5f minutes.\n' % (
                       ff.convert_number_to_string(lamb), mins_l) )
               
            # Time for each channels to run
            t1_c = time.time()
            hours_c = (t1_c-t0_c)/3600
            if ic == 0:
                print('\tDone with S-waves after %.3f hours.\n' % hours_c)
            else:
                print('\tDone with all channels after %.3f hours.\n' % hours_c)
        
        # Time for each potential to run
        t1_k = time.time()
        hours_k = (t1_k-t0_k)/3600
        print( 'Done with kvnn=%d after %.5f hours.\n' % (kvnn, hours_k) )
