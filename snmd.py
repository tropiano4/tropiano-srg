#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: snmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 25, 2021
# 
# Calculates SRG-evolved single-nucleon momentum distributions for nuclei
# assuming the evolved wave function is given by a free Fermi gas. (Note,
# 'snmd' stands for single-nucleon momentum distribution.) Combine these
# functions with LDA.py for nuclear-averaged momentum distributions. This
# script is based off single_particle_momentum_dist.py in Old_codes.
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
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline, UnivariateSpline
# Scripts made by A.T.
from Misc.integration import gaussian_quadrature_mesh
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class single_nucleon_momentum_distributions(object):
    
    
    def __init__(self, kvnn, channels, lamb, kmax=0.0, kmid=0.0, ntot=0):
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
            
        """
        
        # --- Set up --- #
        
        # Save highest allowed L based on input channels
        highest_L = 0
        for channel in channels:
            next_L = vnn.channel_L_value(channel)
            if next_L > highest_L:
                highest_L = next_L
        
        # Load and save momentum and angle arrays for integration
        
        # Relative momentum k [fm^-1]
        k_array, k_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
        if ntot == 0:
            ntot = len(k_array) # Make sure ntot is the number of k-points
        
        # cos(\theta) angles for averaging
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

        # Interpolate pp and pn < k | \delta U | k > matrix elements
        self.deltaU_pp_func = UnivariateSpline( k_array, np.diag(deltaU_pp) )
        self.deltaU_pn_func = UnivariateSpline( k_array, np.diag(deltaU_pn) )
        
        # Interpolate pp and pn < k | \delta U \delta U^{\dagger} | k' > 
        # matrix elements
        self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pp)
        self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pn)
        
        
        # # --- Save the r values and nucleonic densities --- #
        # self.r_array = r_array
        # self.rho_p_array = rho_p_array
        # self.rho_n_array = rho_n_array
        
        
    def select_number_integration_points(self, k_max):
        """
        Select the number of integration points over momenta k given the upper
        limit of integration. Assumes Gaussian quadrature integration.
        
        Parameters
        ----------
        k_max : float
            Upper limit of integration over momentum [fm^-1].
            
        Returns
        -------
        ntot_k : int
            Number of integration points in Gaussian quadrature mesh.
        
        """
        
        # Basing these numbers off expected kF values
        if k_max >= 1.2:
            ntot_k = 60
        elif 1.2 > k_max >= 1.0:
            ntot_k = 50
        elif 1.0 > k_max >= 0.8:
            ntot_k = 40
        elif 0.8 > k_max >= 0.6:
            ntot_k = 30
        elif 0.6 > k_max >= 0.4:
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
        output : 1-D ndarray
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
            if q < kF and 2*k < kF-q:
                theta_deltaU[i] = 1
            
            # Case 2: q < kF and kF-q < 2k < kF+q
            # -> h(q,k) = 1/2 * ( kF^2 - (q-2k)^2 ) / (4*q*k)
            elif q < kF and kF-q < 2*k < kF+q:
                theta_deltaU[i] = 0.5 * ( kF**2 - ( q - 2*k )**2 ) / ( 4*k*q )
                
            # Case 3: q-kF < 2k < kF+q
            # -> h(q,k) = 1/2 * ( kF^2 - (q-2k)^2 ) / (4*q*k)
            elif q-kF < 2*k < kF+q:
                theta_deltaU[i] = 0.5 * ( kF**2 - ( q - 2*k )**2 ) / ( 4*k*q )
                
            # Otherwise, h(q,k) = 0
                
        return theta_deltaU


    def theta_deltaU2(self, kF, K, k_array, ntot_k):
        """
        Evaluates \theta( k_F - \abs(K_vec/2 + k_vec) ) \times
        \theta( k_F - \abs(K_vec/2 - k_vec) ) for every momentum k.
        This function appears in the \delta U \delta U^\dagger term.

        Parameters
        ----------
        kF : float
            Fermi momentum [fm^-1].
        K : float
            C.o.M. momentum [fm^-1].
        k_array : 1-D ndarray
            Momentum values [fm^-1] of \int dk k^2 integral in \delta U term.
        ntot_k : int
            Number of momentum points in momentum mesh k_array.

        Returns
        -------
        output : 1-D ndarray
            \theta function evaluated at each momenta k. This is a unitless
            array of length ntot where ntot = len(k_array).
            
        Notes
        -----
        This function assumes that the k_array does not exceed kF - K/2
        which is the maximum limit given by the \theta function, q < kF_1
        (with the kF in this function corresponding to kF_2), and that K < 2kF.

        """
    
        # Make \theta( k_F - \abs(K_vec/2 +(-) k_vec) ) the same length as
        # k_array
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
    
    
    def n_1(self, q, kF_1):
        """
        Evaluates first term in U n(q) U^\dagger ~ I n(q) I which gives
        2 \theta( kF(r) - q ).
        
        Parameters
        ----------
        q : float
            Single-nucleon momentum [fm^-1].
        kF_1 : float
            Fermi momentum [fm^-1].
            
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

            # Create integration mesh k_array up to kF_2 + q where kF_2 + q
            # corresponds to the upper limit of \theta(kF_2(r)-|q_vec-2k_vec|)
            kmax_delU = kF_2 + q

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

            # Overall factor in front of integral where the first factor of 2
            # is for combining the \delta U and \delta U^\dagger terms and the
            # factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
            deltaU_factor = 2 * 8 * 2/np.pi * 2
            
            # Integrate over \int dk k^2
            return deltaU_factor * np.sum(integrand_k)
        
        else:
            
            return 0


    def n_deltaU2_K(self, q, K, kF):
        """
        Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger
        up to integration over \int dK K^2.
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        K : float
            C.o.M. momentum [fm^-1].
        kF : float
            Fermi momentum [fm^-1].

        Returns
        -------
        output : float
            Momentum distribution from \delta U \delta U^\dagger term before
            integration over \int dr r^2 \int dK K^2.
            
        """
        
        # Create integration mesh k_array up to \sqrt( kF^2 - K^2/4 ) which
        # corresponds to the upper limit of \theta(kF(r)-|K_vec/2 +(-) k_vec|)
        kmax_delU2 = np.sqrt(kF**2-K**2/4)

        # Select number of integration points based on kmax_delU2
        ntot_delU2 = self.select_number_integration_points(kmax_delU2)

        # Get Gaussian quadrature mesh
        k_array_delU2, k_weights_delU2 = gaussian_quadrature_mesh(kmax_delU2,
                                                                  ntot_delU2)
        
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
        
        # Average \theta(kF(r)-|K_vec/2 +(-) k_vec|) functions assuming
        # kF_1=kF_2
        # This is a (ntot_delU2, 1) array
        theta_array = self.theta_deltaU2(kF, K, k_array_delU2, ntot_delU2)
        
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
        
        # For now use minimum kF value and assume symmetric nucleus (N=Z)
        kF = min(kF_1, kF_2)
        
        # Create K_mesh from 0-2*kF
        Kmax = 2*kF # Max momentum
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
        for iK, K in enumerate(K_array):
            
            integrand_K[iK] = self.n_deltaU2_K(q, K, kF)
            
        # Attach integration measure
        integrand_K *= K_array**2 * K_weights
        
        # Overall factor in front of integral
        deltaU2_factor = 1/2 * (2/np.pi)**2 * 2**2
        
        # Integrate over \int dK K^2
        return deltaU2_factor * np.sum(integrand_K)


    def n_lambda_temp(self, q, kF_1, kF_2, contributions='total'):
        
        n_1 = self.n_1(q, kF_1)
        n_delU = self.n_deltaU(q, kF_1, kF_2)
        n_delU2 = self.n_deltaU2(q, kF_1, kF_2)
        total = n_1 + n_delU + n_delU2
        
        if contributions == 'q_contributions':
            return total, n_1, n_delU, n_delU2
        else:
            return total
        
        
    # def n_lambda(self, q_array, r_array, rho_1_array, rho_2_array):
        
    #     # Set shape of return array with 'axes'
    #     axes = 4
            
    #     # Number of r_array points
    #     mtot = len(r_array)
    #     r2_grid, _ = np.meshgrid(r_array**2, np.zeros(axes), indexing='ij')
    #     dr = 0.1 # Spacing between r-points
            
    #     # Length of q_array
    #     ntot = len(q_array)
        
    #     # Evaluate f(q, kF) at each point in q_array and kF
    #     expectation_values = np.zeros( (ntot, axes) )

    #     # Evaluate k_F at each point in r_array
    #     kFp_array = (3*np.pi**2 * rho_p_array)**(1/3)
    #     kFn_array = (3*np.pi**2 * rho_n_array)**(1/3)
            
    #     # Loop over q
    #     for i, q in enumerate(q_array):
    
    #         function_array = np.zeros( (mtot, axes) )
            
    #         # Loop over r for k_F values
    #         for j, (kFp, kFn) in enumerate( zip(kFp_array, kFn_array) ):

    #             function_array[j, :] = func_q(q, kFn, kFp,
    #                                            contributions=contributions)

    #             # Integrate over r for each contribution (summing over axis=0)
    #             expectation_values[i, :] = 4*np.pi*dr * \
    #                                        np.sum(r2_grid * function_array,
    #                                               axis=0)

        
if __name__ == '__main__':
    
    kvnn = 6
    channels = ('1S0', '3S1')
    lamb = 1.35
    kmax, kmid, ntot = 15.0, 3.0, 120
    
    snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax,
                                                  kmid, ntot)
    
    
    # # --- Test \theta function in \delta U term --- #
    # kF = 1.35
    # q = 0.1
    # # q = 0.5
    # # q = 1.4
    # if kF >= 1.2:
    #     ntot_k = 60
    # elif 1.2 > kF >= 1.0:
    #     ntot_k = 50
    # elif 1.0 > kF >= 0.8:
    #     ntot_k = 40
    # elif 0.8 > kF >= 0.6:
    #     ntot_k = 30
    # elif 0.6 > kF >= 0.4:
    #     ntot_k = 20
    # else:
    #     ntot_k = 10
    # k_array, k_weights = gaussian_quadrature_mesh(kF, ntot_k)
    # k_integration_measure = k_weights * k_array**2
    
    # theta_delU = snmd.theta_deltaU(kF, q, k_array, ntot_k)
    # print(theta_delU)
    # print(len(theta_delU))


    # # --- Compare \delta U term to old \delta U term --- #
    # from snmd_v1 import single_nucleon_momentum_distributions as snmd_v1
    # import matplotlib.pyplot as plt
    
    # q_array, _ = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    
    # snmd_old = snmd_v1(kvnn, channels, lamb, kmax, kmid, ntot, interp=False)
    
    # n_new_array = np.zeros(ntot)
    # n_old_array = np.zeros(ntot)
    
    # for iq, q in enumerate(q_array):
        
    #     _, _, delU_new, _ = snmd.n_lambda(q, 1.0, 1.0)
    #     n_new_array[iq] = delU_new
        
    #     _, _, delU_old, _ = snmd_old.n_lambda(q, 1.0, 1.0, 'q_contributions')
    #     n_old_array[iq] = delU_old
        
    # for q, i, j in zip(q_array, n_new_array, n_old_array):
    #     print(q, i, j)
    
    # plt.plot(q_array, n_new_array, label='New')
    # plt.plot(q_array, n_old_array, label='Old')
    # plt.xlabel('q')
    # plt.ylabel('\delta U term')
    # plt.xlim((0.0, 2.0))
    # plt.ylim((-0.3, 0.1))
    # plt.legend(loc='upper left')
    # plt.show()
    
    
    # # --- Test \theta function in \delta U^2 term --- #
    # # kF = 1.35
    # kF = 0.5
    # # K = 0.1*kF
    # # K = 1.0*kF
    # K = 1.9*kF
    # kmax_k = np.sqrt(kF**2-K**2/4)
    # if kmax_k >= 1.2:
    #     ntot_k = 60
    # elif 1.2 > kmax_k >= 1.0:
    #     ntot_k = 50
    # elif 1.0 > kmax_k >= 0.8:
    #     ntot_k = 40
    # elif 0.8 > kmax_k >= 0.6:
    #     ntot_k = 30
    # elif 0.6 > kmax_k >= 0.4:
    #     ntot_k = 20
    # else:
    #     ntot_k = 10
    
    # k_array, k_weights = gaussian_quadrature_mesh(kmax_k, ntot_k)
    # k_integration_measure = k_weights * k_array**2
    
    # theta_delU2 = snmd.theta_deltaU2(kF, K, k_array, ntot_k)
    # print(theta_delU2)
    # print(len(theta_delU2))
    
    
    # --- Compare \delta U^2 term to old \delta U^2 term --- #
    from snmd_v1 import single_nucleon_momentum_distributions as snmd_v1
    import matplotlib.pyplot as plt
    import time
    
    q_array, _ = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    
    snmd_old = snmd_v1(kvnn, channels, lamb, kmax, kmid, ntot, interp=False)
    
    n_new_array = np.zeros(ntot)
    n_old_array = np.zeros(ntot)
    
    # kF = 0.5
    kF = 1.0
    # kF = 1.3

    t0 = time.time()
    for iq, q in enumerate(q_array):
        
        delU2_new = snmd.n_deltaU2(q, kF, kF)
        n_new_array[iq] = delU2_new
        
    t1 = time.time()
    mins = (t1-t0)/60
    print('%.2f mins'%mins)
        
    t0 = time.time()
    for iq, q in enumerate(q_array):
        
        _, _, _, delU2_old = snmd_old.n_lambda(q, kF, kF, 'q_contributions')
        n_old_array[iq] = delU2_old
        
    t1 = time.time()
    mins = (t1-t0)/60
    print('%.2f mins'%mins)
        
    for q, i, j in zip(q_array, n_new_array, n_old_array):
        r = i/j
        print(q, i, j, r)
    
    plt.semilogy(q_array, n_new_array, label='New')
    plt.semilogy(q_array, n_old_array, label='Old')
    plt.xlabel('q')
    plt.ylabel('\delta U^2 term')
    plt.xlim((0.0, 6.0))
    #plt.ylim((-0.3, 0.1))
    plt.legend(loc='upper right')
    plt.show()