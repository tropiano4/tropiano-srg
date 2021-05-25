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
            kmax = round( max(k_array), 0 )
        self.k_array, self.k_weights, self.ntot = k_array, k_weights, ntot
        self.k_integration_measure = k_weights * k_array**2
        
        # TESTING: remove this part since K_mesh will be created for each
        # q, kF_1, kF_2
        # Total momentum K [fm^-1] where we put more points toward K=0 fm^-1
        Kmax = 3.0 # Max momentum
        Kmid = 1.0 # Mid-point
        Ntot = 20 # Total number of points
        Nmod = 10 # Total number of points in the low-K region
        K_array, K_weights = gaussian_quadrature_mesh(Kmax, Ntot, xmid=Kmid,
                                                      nmod=Nmod)
        self.K_array, self.K_weights, self.Ntot = K_array, K_weights, Ntot
        # K integration measure should be a 2-D grid of dimension (ntot, Ntot)
        # since we'll integrate over k last
        _, self.K_integration_measure = np.meshgrid(k_array,
                                        K_weights * K_array**2, indexing='ij')
        
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
        
        
        # --- Set up grids for integration over K and x --- #
        
        # TESTING: remove possibly all of these since you create meshes within
        # n_lambda functions
        
        # Make these grids of k, K, and x for evaluating \theta functions
        self.k_grid_2d, self.x_grid_2d = np.meshgrid(k_array, x_array,
                                                     indexing='ij')
        self.k_grid_3d, self.K_grid_3d, self.x_grid_3d = np.meshgrid(k_array,
                                               K_array, x_array, indexing='ij')
        
        # The next few are grids for attaching integration weights
        _, self.x_weights_2d = np.meshgrid(k_array, x_weights, indexing='ij')
        _, _, self.x_weights_3d = np.meshgrid(k_array, K_array, x_weights,
                                              indexing='ij')
        _, self.K_weights_2d = np.meshgrid(k_array, K_weights, indexing='ij')


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


    def theta_deltaU2(self, kF, sign=1):

        return None
    
    
    def n_1(self, q, kF_1):
        """
        """
        
        if q < kF_1:
            
            return 2
        
        else:
            
            return 0
    
        """
        Evaluates \theta( k_F - \abs(K_vec/2 +(-) k_vec) ) for every momentum k
        and angle x where
            \abs(K_vec/2 +(-) k_vec) = \sqrt(K^2/4 + k^2 +(-) K*k*x).
        These functions appear in the \delta U \delta U^\dagger term.

        Parameters
        ----------
        kF : float
            Fermi momentum [fm^-1].
        sign : int, optional
            Sign of k_vec in the magnitude. Only other option is -1.

        Returns
        -------
        output : 3-D ndarray
            \theta function evaluated at each point in k, K, and x. This is a
            (ntot, Ntot, xtot), unitless array.

        """
    
        K_k_magnitude = np.sqrt( self.K_grid_3d**2/4 + self.k_grid_3d**2 + \
                                 sign*self.k_grid_3d*self.K_grid_3d* \
                                 self.x_grid_3d )
        
        # This returns a (ntot, Ntot, xtot) array of boolean values at every
        # point converted to 1's or 0's by multiplying 1
        theta_K_k = ( K_k_magnitude < kF ) * 1

        # Return the 3-D array with x integration weights
        return theta_K_k * self.x_weights_3d


    def n_lambda(self, q, kF_1, kF_2):
        """
        Single-nucleon momentum distribution where kF_1 specifies proton or
        neutron and kF_2 accounts for correlations to the opposite nucleon.
        (Note, q is not a relative momentum.)
        
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
            
        Return
        ------
        output : tuple
            Single-nucleon momentum distribution [unitless] evaluated at
            momentum q [fm^-1]. (Note, this function will return a tuple of
            floats corresponding to each contribution 1, \delta U, and 
            \delta U^2.)
            
        """
        
        # Split into low- and high-q terms
        
        # Low-q terms
        # Term 1: 1 * n(q) * 1
        # Term 2: \deltaU * n(q) * 1
        # Term 3: 1 * n(q) * \deltaU^\dagger = Term 2
        # \delta U term : Term 2 + Term 3
        if q < kF_1:
            
            term_1 = 2 # 2*\theta(kF_1 - q)
            
            # Create integration mesh k_array up to kF_2+q
            kmax_delU = kF_2 + q

            # Select number of integration points based on kmax_delU
            if kmax_delU >= 1.2:
                ntot_delU = 60
            elif 1.2 > kmax_delU >= 1.0:
                ntot_delU = 50
            elif 1.0 > kmax_delU >= 0.8:
                ntot_delU = 40
            elif 0.8 > kmax_delU >= 0.6:
                ntot_delU = 30
            elif 0.6 > kmax_delU >= 0.4:
                ntot_delU = 20
            else:
                ntot_delU = 10

            # Get Gaussian quadrature mesh
            k_array_delU, k_weights_delU = gaussian_quadrature_mesh(kmax_delU,
                                                                    ntot_delU)
            
            # Average \theta( kF_2 - |q_vec-2k_vec| ) for kF_1 and kF_2
            # (meaning \tau=\tau' and \tau=-\tau')
            theta_pp_array = self.theta_deltaU(kF_1, q, k_array_delU,
                                               ntot_delU)
            theta_pn_array = self.theta_deltaU(kF_2, q, k_array_delU,
                                               ntot_delU)
            
            # Evaluate pp and pn < k | \delta U | k > matrix elements
            deltaU_pp_k = self.deltaU_pp_func(k_array_delU)
            deltaU_pn_k = self.deltaU_pn_func(k_array_delU)
            
            # Build integrand for k integration where we split terms according
            # to pp and pn (or nn and np if kF_1 corresponds to a neutron)
            integrand_k = k_array_delU**2 * k_weights_delU * ( \
                              deltaU_pp_k * theta_pp_array + \
                              deltaU_pn_k * theta_pn_array )

            # Integrate over k where the first factor of 2 is for combining the
            # \delta U and \delta U^\dagger terms and factor of 8 from
            # evaluating \int d^3K \delta(K/2 - ...)
            deltaU_factor = 2 * 8 * 2/np.pi * 2
            term_deltaU = deltaU_factor * np.sum(integrand_k)
            
        # q > kF_1
        else:
            
            term_1 = 0
            term_deltaU = 0
        
        # # High-q term: \deltaU * n(q) * \deltaU^\dagger

        # # Create a grid for evaluation of \delta U( k, abs(q_vec - K_vec/2) )
        # # * \delta U^\dagger( abs(q_vec - K_vec/2), k )
        # # Dimensions (ntot, Ntot, xtot)
        # q_K_grid = np.sqrt( q**2 + self.K_grid_3d**2/4 - 
        #                     q*self.K_grid_3d*self.x_grid_3d )
        
        # # Evaluate \delta U^2( k, \abs(q_vec + K_vec/2) ) and integrate over
        # # angles
        # deltaU2_pp_x = self.deltaU2_pp_func.ev(self.k_grid_3d, q_K_grid) * \
        #                self.x_weights_3d
        # deltaU2_pn_x = self.deltaU2_pn_func.ev(self.k_grid_3d, q_K_grid) * \
        #                self.x_weights_3d
        
        # # Integrate over angles
        # # These are now (ntot, Ntot)
        # deltaU2_pp_array = np.sum(deltaU2_pp_x, axis=-1)/2
        # deltaU2_pn_array = np.sum(deltaU2_pn_x, axis=-1)/2
        
        # # Build x-dependent part and integrate with resized K
        
        # # These are (ntot, Ktot, xtot) arrays of
        # # \theta( k_F - \abs( 1/2*K_vec +/- k_vec ) ) for every total momentum
        # # K, relative momentum k, and angle x where sign specifies the sign of
        # # k_vec
        # theta_kF1_K_plus_k_x = self.theta_deltaU2(kF_1, sign=1)
        # theta_kF1_K_minus_k_x = self.theta_deltaU2(kF_1, sign=-1)
        # # theta_kF2_K_plus_k_x = self.theta_deltaU2(kF_2, sign=1)
        # theta_kF2_K_minus_k_x = self.theta_deltaU2(kF_2, sign=-1)
            
        # # Integrate over x first
        # # This sum collapses theta_K_k_x to a matrix dependent only on K and 
        # # k: (ntot, Ktot, xtot) -> (ntot, Ktot)
            
        # # \int dx/2 \theta_pp (or \theta_nn)
        # theta_kF1_kF1_K_k = np.sum(theta_kF1_K_plus_k_x*theta_kF1_K_minus_k_x,
        #                            axis=-1) / 2
        # # Take out the extra 1/2 in the isospin CG's and remove second term in
        # # theta_kF1_kF2_K_k for better clarity
        # theta_kF1_kF2_K_k = np.sum( 
        #                     theta_kF1_K_plus_k_x * theta_kF2_K_minus_k_x,
        #                     axis=-1 ) / 2
        
        # # Overall factor in front of \delta U^2 term
        # # deltaU2_factor = 1/2 * (2/np.pi)**2 * 2**4
        # # Using expression from Dickhoff and Overleaf appendix
        # deltaU2_factor = 1/2 * (2/np.pi)**2 * 2**2

        # # Build K integrand (ntot, K_cutoff_index) array
        # integrand_k_K = self.K_integration_measure * (\
        #                     deltaU2_pp_array * theta_kF1_kF1_K_k + \
        #                     deltaU2_pn_array * theta_kF1_kF2_K_k )

        # # Integrate over K and build k integrand
        # integrand_k = np.sum( integrand_k_K, axis=-1 ) * \
        #               self.k_integration_measure
                      
        # # Integrate over k
        # term_deltaU2 = deltaU2_factor * np.sum(integrand_k)
        
        term_deltaU2 = 0.0
        
        # Add up each contribution for total
        total = term_1 + term_deltaU + term_deltaU2

        # Return contributions and total or just total
        return total, term_1, term_deltaU, term_deltaU2

        
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

    # --- Compare \delta U term to old \delta U term --- #
    from snmd_v1 import single_nucleon_momentum_distributions as snmd_v1
    import matplotlib.pyplot as plt
    
    q_array, _ = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    
    snmd_old = snmd_v1(kvnn, channels, lamb, kmax, kmid, ntot, interp=False)
    
    n_new_array = np.zeros(ntot)
    n_old_array = np.zeros(ntot)
    
    for iq, q in enumerate(q_array):
        
        _, _, delU_new, _ = snmd.n_lambda(q, 1.0, 1.0)
        n_new_array[iq] = delU_new
        
        _, _, delU_old, _ = snmd_old.n_lambda(q, 1.0, 1.0, 'q_contributions')
        n_old_array[iq] = delU_old
        
    for q, i, j in zip(q_array, n_new_array, n_old_array):
        print(q, i, j)
    
    plt.plot(q_array, n_new_array, label='New')
    plt.plot(q_array, n_old_array, label='Old')
    plt.xlabel('q')
    plt.ylabel('\delta U term')
    plt.xlim((0.0, 2.0))
    plt.ylim((-0.3, 0.1))
    plt.legend(loc='upper left')
    plt.show()