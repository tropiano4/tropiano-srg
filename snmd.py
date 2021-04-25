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
# Revision history:
#   03/29/21 --- Finalized and produces same results (up to factors of 2*\pi)
#                as single_particle_momentum_dist.py. Also, fixed a bug
#                involving the integration over total momentum K.
#   03/31/21 --- Moved channel_L_value function to vnn.py. See vnn.py for
#                details of the function.
#   04/02/21 --- Added option to return pp and pn (or nn and np) contributions
#                to \delta U \delta U^\dagger term along with total.
#   04/22/21 --- Testing changes to overall factors in front of \delta U
#                and \delta U^2 terms. Update this update when it's correct!
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline
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
            Maximum value in the momentum mesh [fm^-1].
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
            ntot = len(k_array)
        self.k_array, self.k_weights, self.ntot = k_array, k_weights, ntot
        self.k_integration_measure = k_weights * k_array**2
        
        # Total momentum K [fm^-1] where we put more points toward K=0 fm^-1
        Kmax = 3.0 # Max momentum
        Kmid = 1.5 # Mid-point
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
                deltaU2_pn += (2*J+1)/4 * ( delta_U_matrix[:ntot, :ntot]**2 \
                                          + delta_U_matrix[:ntot, ntot:]**2 )

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

        # Evaluate pp and pn < k | \delta U | k > terms
        # These are (ntot, 1) arrays
        self.deltaU_pp_k = np.diag(deltaU_pp)
        self.deltaU_pn_k = np.diag(deltaU_pn)
        
        # Interpolate pp and pn < k | \delta U \delta U^{\dagger} | k' > 
        # matrix elements to evaluate at |q_vec-K_vec/2|
        self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pp)
        self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pn)
        
        
        # --- Set up grids for integration over K and x --- #
        
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


    def theta_q_2k(self, kF, q):
        """
        Evaluates \theta( k_F - \abs(q_vec - 2k_vec) ) for every momentum k
        and angle x where \abs(q_vec - 2k_vec) = \sqrt(q^2 + 4k^2 + 4kqx).

        Parameters
        ----------
        kF : float
            Fermi momentum [fm^-1].
        q : float
            Single-nucleon momentum [fm^-1].

        Returns
        -------
        output : 2-D ndarray
            \theta function evaluated at each point in k and x. This is a
            (ntot, xtot), unitless array.

        """
    
        q_2k_magnitude = np.sqrt( q**2 + 4 * self.k_grid_2d**2 - \
                                  4 * q * self.k_grid_2d * self.x_grid_2d )
        
        # This returns a (ntot, xtot) array of boolean values at every point
        # converted to 1's or 0's by multiplying 1
        theta_q_2k = ( q_2k_magnitude < kF ) * 1

        # Return the weighted matrix for integration purposes
        return theta_q_2k * self.x_weights_2d
    
    
    def theta_K_k(self, kF, sign=1):
        """
        Evaluates \theta( k_F - \abs(K_vec/2 + sign*k_vec) ) for every
        momentum k and angle x where \abs(K_vec/2 + sign * k_vec) = 
        \sqrt(K^2/4 + k^2 + sign*Kkx).

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
                                 sign * self.k_grid_3d * self.K_grid_3d * \
                                 self.x_grid_3d )
        
        # This returns a (ntot, Ntot, xtot) array of boolean values at every
        # point converted to 1's or 0's by multiplying 1
        theta_K_k = ( K_k_magnitude < kF ) * 1

        # Return the weighted matrix for integration purposes
        return theta_K_k * self.x_weights_3d


    def n_lambda(self, q, kF_1, kF_2, contributions='total'):
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
        contributions : str, optional
            Option to return different contributions to the momentum
            distribution.
            1. Default is 'total' which only returns the total momentum
               distribution.
            2. Specify 'NN_contributions' for total, pp, and pn (or nn, np)
               where the nucleon-nucleon contributions are isolated in the
               high-q term.
            3. Specify 'q_contributions' for total, 1, \delta U, and 
               \delta U \delta U^\dagger.
            
        Return
        ------
        output : float or tuple
            Single-nucleon momentum distribution [unitless] evaluated at
            momentum q [fm^-1]. (Note, this function will return a tuple of
            floats if contributions is not 'total'.)
            
        """
        
        # Split into low- and high-q terms
        
        # Low-q terms
        # Term 1: 1 * n(q) * 1
        # Term 2: \deltaU * n(q) * 1
        # Term 3: 1 * n(q) * \deltaU^\dagger = Term 2
        # \delta U term : Term 2 + Term 3
        if q < kF_1:
            
            term_1 = 2 # 2*\theta(kF_1 - q)
            
            # These are (ntot, xtot) matrices (including x_weights) of 
            # \theta( kF1 - \abs( q_vec - 2k_vec ) ) for every momentum k and 
            # angle x
            theta_kF1_k_x_matrix = self.theta_q_2k(kF_1, q)
            theta_kF2_k_x_matrix = self.theta_q_2k(kF_2, q)
            
            # Integrate over x first (angle-averaging) where the integration
            # weights of x are already built-in
            # This sum collapses theta_k_x_matrix to a vector dependent only
            # on k: (ntot, xtot) -> (ntot, 1)
            theta_kF1_k_vector = np.sum( theta_kF1_k_x_matrix, axis=-1 ) / 2
            theta_kF2_k_vector = np.sum( theta_kF2_k_x_matrix, axis=-1 ) / 2
            
            # Build integrand for k integration where we split terms according
            # to pp and pn (or nn and np if kF_1 corresponds to a neutron)
            integrand_k = self.k_integration_measure * ( \
                              self.deltaU_pp_k * theta_kF1_k_vector + \
                              self.deltaU_pn_k * theta_kF2_k_vector )

            # Integrate over k where the factor of 2 is for combining the
            # second and third terms
            deltaU_factor = 2 * 2/np.pi * 2**2
            term_deltaU = deltaU_factor * np.sum(integrand_k)
            
        # q > kF_1
        else:
            
            term_1 = 0
            term_deltaU = 0
        
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
        
        # Overall factor in front of \delta U^2 term
        deltaU2_factor = 1/2 * (2/np.pi)**2 * 2**4
        
        # Split pp and np up to isolate contributions
        if contributions == 'NN_contributions':
            
            # Build K integrand (ntot, K_cutoff_index) array isolating pp and
            # pn terms
            integrand_k_K_pp = self.K_integration_measure[:, :K_cutoff_index] * \
                               deltaU2_pp_array * theta_kF1_kF1_K_k
            integrand_k_K_pn = self.K_integration_measure[:, :K_cutoff_index] * \
                               deltaU2_pn_array * theta_kF1_kF2_K_k
                               
            # Integrate over K and build k integrand
            integrand_k_pp = np.sum( integrand_k_K_pp, axis=-1 ) * \
                             self.k_integration_measure
            integrand_k_pn = np.sum( integrand_k_K_pn, axis=-1 ) * \
                             self.k_integration_measure
                             
            # Integrate over k
            term_deltaU2_pp = deltaU2_factor * np.sum(integrand_k_pp)
            term_deltaU2_pn = deltaU2_factor * np.sum(integrand_k_pn)
            
            # Add for total \delta U^2 term then calculate total
            term_deltaU2 = term_deltaU2_pp + term_deltaU2_pn
            total = total = term_1 + term_deltaU + term_deltaU2
            
            # Return \delta U^2 pp and pn (or nn and np) contributions along
            # with total
            return total, term_deltaU2_pp, term_deltaU2_pn
            
        # Otherwise, add pp and pn
        else:

            # Build K integrand (ntot, K_cutoff_index) array
            integrand_k_K = self.K_integration_measure[:, :K_cutoff_index] * (\
                            deltaU2_pp_array * theta_kF1_kF1_K_k + \
                            deltaU2_pn_array * theta_kF1_kF2_K_k )

            # Integrate over K and build k integrand
            integrand_k = np.sum( integrand_k_K, axis=-1 ) * \
                          self.k_integration_measure
                      
            # Integrate over k
            term_deltaU2 = deltaU2_factor * np.sum(integrand_k)
        
            # Add up each contribution for total
            total = term_1 + term_deltaU + term_deltaU2

            # Return contributions and total or just total
            if contributions == 'q_contributions':
                return total, term_1, term_deltaU, term_deltaU2
            else: # Default
                return total
    
    
if __name__ == '__main__':
    
    
    # --- Set up --- #
    
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText
    import time
    from lda import load_density, LDA
    
    # AV18 with \lambda=1.35 fm^-1
    kvnn = 6
    lamb = 1.35
    channels = ('1S0', '3S1')
    # channels = ('1S0', '3S1', '3P0', '1P1', '3P1', '3P2')
    # channels = ('1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3')
    kmax, kmid, ntot = 10.0, 2.0, 120
    # kmax, kmid, ntot = 30.0, 4.0, 120
    q_array, q_weights = vnn.load_momentum(kvnn, '1S0', kmax, kmid, ntot)
    factor_array = 2/np.pi * q_weights * q_array**2
    
    # Load n_\lambda_pp(q, k_F) for AV18
    snmd = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax,
                                                 kmid, ntot)
    
    # Details of example nuclei (format is [nuclei, Z, N])
    # nuclei_list = ['C12', 6, 6]
    # nuclei_list = ['O16', 8, 8]
    nuclei_list = ['Ca40', 20, 20]
    # nuclei_list = ['Ca48', 20, 28]
    # nuclei_list = ['Fe56', 26, 30]
    # nuclei_list = ['Pb208', 82, 126]
    nucleus = nuclei_list[0]
    Z = nuclei_list[1]
    N = nuclei_list[2]
    A = N + Z
    r_array, rho_p_array = load_density(nucleus, 'proton', Z, N)
    r_array, rho_n_array = load_density(nucleus, 'neutron', Z, N)
    
    # Call LDA class
    lda = LDA(r_array, rho_p_array, rho_n_array)
    
    
    # --- Calculate momentum distributions --- #
    
    # n(q) p single-particle
    t0 = time.time() # Start time
    n_p_array = lda.local_density_approximation( q_array, snmd.n_lambda, 'p' )
    t1 = time.time() # End time
    mins = round( (t1 - t0) / 60.0, 2)
    print('n_lambda^p(q) done after %.2f minutes' % mins)
    
    t0 = time.time()
    n_n_array = lda.local_density_approximation( q_array, snmd.n_lambda, 'n' )
    t1 = time.time()
    print('n_lambda^n(q) done after %.2f minutes' % mins)
    mins = round( (t1 - t0) / 60.0, 2)
    
    # This is 4 \pi for d3p and 1/(2*\pi)^3 for converting
    # from sums to integrals
    overall_factor = 4*np.pi * 1/(2*np.pi)**3
    
    p_a_p = q_array**2 * n_p_array / A * overall_factor
    p_a_n = q_array**2 * n_n_array / A * overall_factor

    p_a = p_a_p + p_a_n

    print( np.sum(p_a*q_weights) )
    
    
    # --- Plot --- #
    
    plt.clf()
    f, ax = plt.subplots( 1, 1, figsize=(4, 4) )
    
    # n(q) p single-particle
    plt.semilogy(q_array, p_a_p, color='red', linestyle='dashdot', 
                  label=nucleus+' (p)')
    plt.semilogy(q_array, p_a_n, color='blue', linestyle='dashed', 
                  label=nucleus+' (n)')
    plt.semilogy(q_array, p_a, color='black', label=nucleus+' (total)')
    plt.ylabel(r'$P^A(q)$')
    ax.set_xlim( [min(q_array), 4] )
    ax.set_ylim( [1e-9, 2e0])
    
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('q [fm' + r'$^{-1}$' + ']')
    
    lambda_label = 'AV18\n' + r'$\lambda=%.2f$' % lamb + ' fm' + r'$^{-1}$'
    lambda_label_location = 'center right'
    anchored_text = AnchoredText(lambda_label, loc=lambda_label_location, 
                                  frameon=False)
    ax.add_artist(anchored_text)
    
    plt.show()
    
    # --- Further tests and bug checks --- #
    # 1. How much do things change with different K and x meshes?
    # 2. Check that the q_K_cutoff is working correctly. (This is still giving
    #    me 0 at the end of the final nuclear-averaged array.)