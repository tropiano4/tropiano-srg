#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: dmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 30, 2021
# 
# Calculates SRG-evolved deuteron momentum distributions. (Note, 'dmd' stands
# for deuteron momentum distribution.) Combine these functions with lda.py for
# nuclear-averaged momentum distributions which can then be used for
# calculation of A/d ratios with snmd.py.
#
# Revision history:
#   03/31/21 --- Moved hankel_transformation function to fourier_transform.py.
#   04/06/21 --- Added pair momentum distribution function called
#                'n_lambda_pair' and renamed original single-nucleon momentum
#                distribution from 'n_lambda' to 'n_lambda_single'. Also
#                merged pmd_deuteron_test.py with this script.
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
from Misc.fourier_transform import hankel_transformation_k2r
from Misc.integration import gaussian_quadrature_mesh
import observables as ob
from operators import find_q_index, momentum_projection_operator
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class deuteron_momentum_distributions(object):
    
    
    def __init__( self, kvnn, lamb, kmax=0.0, kmid=0.0, ntot=0,
                  r_array=np.linspace(0.1, 20.0, 200) ):
        """
        Saves momentum arrays, grids, and matrix elements of \delta U and
        \delta U^{\dagger} given the input potential and SRG \lambda.
        Evaluates and saves matrix elements. Computes nucleonic density in
        deuteron given wave function.
        
        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        lamb : float
            SRG evolution parameter lambda [fm^-1].
        kmax : float, optional
            Maximum value in the momentum mesh [fm^-1].
        kmid : float, optional
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int, optional
            Number of momentum points in mesh.
        r_array : 1-D ndarray, optional
            Relative position values [fm]. Default is the same as r_array from
            Densities codes. (Note, we are assuming a linearly-spaced array.)
            
        """
        
        # --- Set up --- #
        
        # Channel is 3S1-3D1 for deuteron
        channel = '3S1'
        
        # Load and save momentum and angle arrays for integration
        
        # Relative momentum k [fm^-1]
        k_array, k_weights = vnn.load_momentum(kvnn, channel)
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

        # Load SRG transformation
        H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
        H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
                                         method='srg', generator='Wegner', 
                                         lamb=lamb)
        # Load U(k, k') [unitless]
        U_matrix_unitless = SRG_unitary_transformation(H_initial, H_evolved)

        # Isolate 2-body term, save, and convert to fm^3
        I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
        row, col = np.meshgrid(factor_array_cc, factor_array_cc)
        delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
        self.delta_U_matrix = delta_U_matrix_unitless
        delta_U_matrix = delta_U_matrix_unitless / row / col
            
        # 2J+1 factor
        J = 1
            
        # Evaluate matrix elements
        deltaU = (2*J+1)/2 * ( delta_U_matrix[:ntot, :ntot] + \
                               delta_U_matrix[ntot:, ntot:] )
        deltaU2 = (2*J+1)/4 * ( delta_U_matrix[:ntot, :ntot]**2 + \
                                delta_U_matrix[:ntot, ntot:]**2 + \
                                delta_U_matrix[ntot:, :ntot]**2 + \
                                delta_U_matrix[ntot:, ntot:]**2 )
            
        # Save \delta U matrix elements (these are (ntot, ntot) arrays)
        # < k | \delta U | k' >
        self.deltaU = deltaU
        # < k | \delta U \delta U^\dagger | k' >
        self.deltaU2 = deltaU2

        # Evaluate < k | \delta U | k > term
        # This a (ntot, 1) arrays
        self.deltaU_k = np.diag(deltaU)

        # Interpolate < k | \delta U \delta U^{\dagger} | k' > matrix elements
        # to evaluate at |q_vec-K_vec/2|
        self.deltaU2_func = RectBivariateSpline(k_array, k_array, deltaU2)
        
        
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
        
        
        # --- Set up deuteron density for averaging --- #
        
        # Unitless wave function in momentum space
        psi_k_unitless = ob.wave_function(H_initial, U=U_matrix_unitless)
        self.psi = psi_k_unitless
    
        # Divide out momenta/weights
        psi_k = psi_k_unitless / factor_array_cc # [fm^3/2]
    
        # Transform to coordinate space
        hank_trans_3S1 = hankel_transformation_k2r('3S1', k_array, k_weights,
                                                   r_array)
        hank_trans_3D1 = hankel_transformation_k2r('3D1', k_array, k_weights,
                                                   r_array)
    
        # Get 3S1 and 3D1 waves [fm^-3/2]
        psi_r_3S1 = hank_trans_3S1 @ psi_k[:ntot]
        psi_r_3D1 = hank_trans_3D1 @ psi_k[ntot:]
    
        # Save deuteron nucleonic density [fm^-3], r_array [fm], and dr [fm]
        self.rho_array = psi_r_3S1**2 + psi_r_3D1**2
        self.r_array = r_array
        self.dr = r_array[1] - r_array[0]
        
        
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
    
    
    def n_lambda_single(self, q, kF, contributions='total'):
        """
        Single-nucleon momentum distribution of deuteron given a Fermi
        momentum for both the proton and neutron. kF(r) is set by the deuteron
        relative position wave function \psi(r). (Note, q is not a relative
        momentum.)
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        kF : float
            Fermi momentum [fm^-1] for the proton and neutron.
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
        if q < kF:
            
            term_1 = 2 # 2*\theta(kF_1 - q)
            
            # This is a (ntot, xtot) matrix (including x_weights) of 
            # \theta( kF - \abs( q_vec - 2k_vec ) ) for every momentum k and 
            # angle x
            theta_kF_k_x_matrix = self.theta_q_2k(kF, q)
            
            # Integrate over x first (angle-averaging) where the integration
            # weights of x are already built-in
            # This sum collapses theta_k_x_matrix to a vector dependent only
            # on k: (ntot, xtot) -> (ntot, 1)
            theta_kF_k_vector = np.sum( theta_kF_k_x_matrix, axis=-1 ) / 2
            
            # Build integrand for k integration
            integrand_k = self.k_integration_measure * self.deltaU_k * \
                          theta_kF_k_vector

            # Integrate over k where the factor of 2 is for combining the
            # second and third terms
            deltaU_factor = 2 / (2*np.pi)**3 * 2/np.pi
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
        deltaU2_array = self.deltaU2_func.ev(k_grid, q_K_grid)
        
        # Build x-dependent part and integrate with resized K
        
        # These are (ntot, Ktot, xtot) arrays of
        # \theta( k_F - \abs( 1/2*K_vec +/- k_vec ) ) for every total momentum
        # K, relative momentum k, and angle x where sign specifies the sign of
        # k_vec
        theta_kF_K_plus_k_x = self.theta_K_k(kF, sign=1)
        theta_kF_K_minus_k_x = self.theta_K_k(kF, sign=-1)
            
        # Integrate over x first
        # This sum collapses theta_K_k_x to a matrix dependent only on K and 
        # k: (ntot, Ktot, xtot) -> (ntot, K_cutoff_index)
            
        # \int dx/2
        theta_kF_K_k = ( np.sum( 
                              theta_kF_K_plus_k_x * theta_kF_K_minus_k_x,
                              axis=-1 ) )[:, :K_cutoff_index] / 2

        # Build K integrand (ntot, K_cutoff_index) array spliting pp and np
        # contributions (or nn and np)
        integrand_k_K = self.K_integration_measure[:, :K_cutoff_index] * \
                        deltaU2_array * theta_kF_K_k

        # Integrate over K and build k integrand
        integrand_k = np.sum( integrand_k_K, axis=-1 ) * \
                      self.k_integration_measure
                      
        # Integrate over k (no 1/2 factor as in snmd.py because the two pn
        # \theta's are the same in this case and give a factor of 2)
        deltaU2_factor = 1/(2*np.pi)**6 * (2/np.pi)**2
        term_deltaU2 = deltaU2_factor * np.sum(integrand_k)
        
        # Add up each contribution for total
        total = term_1 + term_deltaU + term_deltaU2

        # Return contributions and total or just total
        if contributions == 'q_contributions':
            return total, term_1, term_deltaU, term_deltaU2
        else: # Default
            return total
    

    def n_lambda_pair(self, q, kF, contributions='total'):
        """
        Pair momentum distribution of deuteron given a Fermi momentum for both
        the proton and neutron. kF(r) is set by the deuteron relative position
        wave function \psi(r).
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        kF : float
            Fermi momentum [fm^-1] for the proton and neutron.
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
            
            deltaU_factor = 1/(4*np.pi) * 2/np.pi
            # delta U evaluated at q
            term_deltaU = deltaU_factor * self.deltaU[q_index, q_index]
            
        # q > kF_1
        else:
            
            term_1 = 0
            term_deltaU = 0

        # High-q term: \deltaU * n(q) * \deltaU^\dagger
        
        # Evaluate \delta(k, q) * \delta U^\dagger(q, k) and multiply by
        # dk k^2
        integrand_k = self.deltaU2[:, q_index] * self.k_integration_measure

        # Index of where to stop the integral in the momentum mesh k_array
        # (This is from \theta(kFp-k))
        kF_cutoff = find_q_index(kF, self.k_array)
                      
        # Integrate over k
        deltaU2_factor = 1/4 * 1/(2*np.pi)**3 * (2/np.pi)**2 * 1/(4*np.pi) * 2
        # Last factor of 2 is for \theta^p \theta^n + \theta^p \theta^n
        term_deltaU2 = deltaU2_factor * np.sum( integrand_k[:kF_cutoff] )
        
        # Add up each term for total
        total = term_1 + term_deltaU + term_deltaU2

        # Return contributions and total or just total
        if contributions == 'q_contributions':
            return total, term_1, term_deltaU, term_deltaU2
        else: # Default
            return total
    

    def local_density_approximation(self, q_array, distribution):
        """
        Evaluates the nuclear-averaged single-nucleon momentum distribution in
        deuteron assuming LDA.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1]. These are not relative momenta.
        distribution : str
            Type of distribution: 'single-nucleon' or 'pair'.

        Returns
        -------
        expectation_values : 1-D ndarray
            Array of expectation values of the distribution evaluated at each
            momentum q.

        """
        
        # Number of q points
        M = len(q_array)
        
        # Load r_array and rho_array
        r_array = self.r_array
        rho_array = self.rho_array
        # Number of r points
        N = len(r_array)
    
        # kF values for pn pair
        kF_array = ( 3*np.pi**2 * rho_array )**(1/3)
        
        # Loop over q_array and average with respect to r
        expectation_values = np.zeros(M)
        for i, q in enumerate(q_array):
    
            # Now evaluate n(q, kF)
            n_lambda_array = np.zeros(N)
            for j, k_F in enumerate(kF_array):

                if distribution == 'single-nucleon':
                    n_lambda_array[j] = self.n_lambda_single(q, k_F)
                elif distribution == 'pair':
                    n_lambda_array[j] = self.n_lambda_pair(q, k_F)
                    
            # Deuteron is normalization is \int dr r^2 \rho_d(r) = 1 without
            # any factor of 4*\pi
            expectation_values[i] = self.dr * np.sum( r_array**2 *
                                                      n_lambda_array )
  
        return expectation_values
    
    
    def n_lambda_pair_exact(self, q, contributions='total'):
        """
        Computes the exact pair momentum distribution in deuteron using the
        wave function from the input Hamiltonian.

        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        contributions : str, optional
            Option to return different contributions to the pair momentum
            distribution.
            1. Default is 'total' which only returns the total momentum
               distribution.
            2. Specify 'q_contributions' for 1, \delta U, and 
               \delta U \delta U^\dagger along with the total.
            3. Specify 'partial_wave_ratio' for the ratio of 3S1-3S1 to full
               3S1-3D1 along with the absolute total.

        Returns
        -------
        output : float or tuple
            Pair momentum distribution [unitless] evaluated at momentum q
            [fm^-1]. (Note, this function will return a tuple of floats if
            contributions is not 'total'.)

        """

        # Load deuteron wave function and delta_U_matrix (both unitless)
        psi_vector = self.psi
        delta_U_matrix = self.delta_U_matrix
        
        # Load bare momentum projection operator [fm^3]
        bare_operator = momentum_projection_operator(q, self.k_array,
                                                      self.k_weights, '3S1',
                                                      smeared=False)
        
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
            
            # Find index of q
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
    
    
if __name__ == '__main__':
    
    
    # --- Compare pair momentum distribution to exact result --- #
    
    import matplotlib.pyplot as plt
    
    channel = '3S1'
    
    kvnn = 6
    lamb = 1.35
    kmax, kmid, ntot = 10.0, 2.0, 120
    
    # Load momentum and weights
    q_array, q_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
    factor_array = q_array**2 * q_weights
    
    # Load hamiltonian
    H_matrix = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    
    # Load exact wave function [unitless]
    psi_exact_unitless = ob.wave_function(H_matrix)
    psi_squared_exact = ( psi_exact_unitless[:ntot]**2 + \
                          psi_exact_unitless[ntot:]**2 ) / factor_array
        
    # Calculate using LDA
    dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot)
    n_d_array = dmd.local_density_approximation(q_array, 'pair')
    
    # Plot pair momentum distributions
    plt.semilogy(q_array, psi_squared_exact, label='Exact')
    plt.semilogy(q_array, n_d_array, label='LDA')
    plt.semilogy(q_array, n_d_array * (2*np.pi)**3,
                  label=r'$(2\pi)^3$' + 'LDA')
    plt.xlim( (0.0, 5.0) )
    plt.ylim( (1e-5, 1e3) )
    plt.xlabel(r'$q$' + ' [fm' + r'$^{-1}$' + ']')
    plt.ylabel(r'$n_d(q)$' + ' [fm' + r'$^3$' + ']')
    plt.legend(loc=0)