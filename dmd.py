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
# Notes on normalizations:
#   1. The deuteron wave function describing relative position or momentum is
#      normalized according to
#        \int dr r^2 ( |\psi_3S1(r)|^2 + |\psi_3D1(r)|^2 ) = 1,
#        2/\pi * \int dk k^2 ( |\psi_{3S1}(k)|^2 + |\psi_{3D1}(k)|^2 ) = 1.
#   2. Under LDA, we adopt the normalization
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
    
    
    def __init__(self, kvnn, lamb, kmax=0.0, kmid=0.0, ntot=0, interp=True):
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
        interp : bool, optional
            Option to use interpolated n_\lambda(q) functions.
            
        """
        
        # Part of the data directory name
        self.kvnn = kvnn
        # Part of the file name
        self.lamb = lamb
        self.kmax = kmax
        
        # Calculate directly and/or generate data files
        if interp == False:
        
            # Channel is 3S1-3D1 for deuteron
            channel = '3S1'
        
            # Relative momentum k [fm^-1]
            k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid,
                                                   ntot)
            self.k_array, self.k_weights = k_array, k_weights
            if ntot == 0:
                ntot = len(k_array) # Make sure ntot is the number of k-points
            self.ntot = ntot
        
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
            U_matrix_unitless = SRG_unitary_transformation(H_initial,
                                                           H_evolved)

            # Isolate 2-body term, save, and convert to fm^3
            I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
            row, col = np.meshgrid(factor_array_cc, factor_array_cc)
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            # Save for exact unitary calculation
            self.delta_U_matrix = delta_U_matrix_unitless
            delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # No 2*J+1 factor since we're fixing M_J for deuteron
            
            # Evaluate matrix elements
            deltaU = 1/2 * ( delta_U_matrix[:ntot, :ntot] + \
                             delta_U_matrix[ntot:, ntot:] )
            deltaU2 = 1/2 * ( delta_U_matrix[:ntot, :ntot]**2 + \
                              delta_U_matrix[:ntot, ntot:]**2 + \
                              delta_U_matrix[ntot:, :ntot]**2 + \
                              delta_U_matrix[ntot:, ntot:]**2 )

            # Interpolate < k | \delta U | k > matrix elements
            # self.deltaU_func = UnivariateSpline( k_array, np.diag(deltaU) )
            self.deltaU_func = interp1d( k_array, np.diag(deltaU),
                                         bounds_error=False,
                                         fill_value='extrapolate' )

            # Interpolate < k | \delta U \delta U^{\dagger} | k' > matrix
            # elements to evaluate at |q_vec-K_vec/2|
            self.deltaU2_func = RectBivariateSpline(k_array, k_array, deltaU2)
        
        
            # --- Set up deuteron density for averaging --- #
        
            # Unitless wave function in momentum space
            psi_k_unitless = ob.wave_function(H_initial, U=U_matrix_unitless)
            self.psi = psi_k_unitless
    
            # Divide out momenta/weights and save
            self.psi_k = psi_k_unitless / factor_array_cc # [fm^3/2]

        
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
            if q < kF and 2*k <= kF-q:
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
    
    
    def n_1(self, q, kF):
        """
        Evaluates first term in U n(q) U^\dagger ~ I n(q) I which gives
        2 \theta( kF(r) - q ).
        
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
        
        # This is simply a \theta function where the factor of two comes from
        # summing over particle spin projection \sigma
        if q < kF:
            
            return 2
        
        else:
            
            return 0
        
        
    def n_deltaU(self, q, kF):
        """
        Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
        we are combining \delta U and \delta U^\dagger, hence the factor of 2.
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        kF : float
            Fermi momentum [fm^-1].
            
        Returns
        -------
        output : float
            Momentum distribution from \delta U term before integration over
            \int dr r^2.
            
        """

        # Initial \theta( kF(r) - q )
        if q < kF:

            # Create integration mesh k_array up to (kF + q)/2 where
            # (kF + q)/2 corresponds to the upper limit of 
            # \theta(kF(r)-|q_vec-2k_vec|)
            kmax_delU = (kF + q)/2

            # Select number of integration points based on kmax_delU
            ntot_delU = self.select_number_integration_points(kmax_delU)

            # Get Gaussian quadrature mesh
            k_array_delU, k_weights_delU = gaussian_quadrature_mesh(kmax_delU,
                                                                    ntot_delU)
            
            # Average \theta( kF - |q_vec-2k_vec| )
            theta_array = self.theta_deltaU(kF, q, k_array_delU, ntot_delU)
            
            # Evaluate < k | \delta U | k > matrix elements
            deltaU_k = self.deltaU_func(k_array_delU)
            
            # Build integrand for k integration
            integrand_k = k_array_delU**2 * k_weights_delU * deltaU_k * \
                          theta_array

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
        
        # Create integration mesh k_array from minimum and maximum k values
        # corresponding to the limits of \theta(kF(r)-|K_vec/2 +(-) k_vec|)
        kmin_delU2 = max(K/2 - kF, 0)
        kmax_delU2 = min( np.sqrt( kF**2 - K**2/4 ), kF + K/2 )

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
        # This is a (ntot_delU2, xtot) array
        deltaU2_k_x = self.deltaU2_func.ev(k_grid, q_K_grid) * \
                      x_weights_grid
        
        # Integrate over \int dx/2 -> (ntot_k, 1)
        deltaU2_k = np.sum(deltaU2_k_x, axis=-1)/2
        
        # Average \theta(kF(r)-|K_vec/2 +(-) k_vec|) functions
        # This is a (ntot_delU2, 1) array
        theta_array = self.theta_deltaU2(kF, K, k_array_delU2, ntot_delU2)
        
        # Build integrand for k integration
        integrand_k = k_array_delU2**2 * k_weights_delU2 * deltaU2_k * \
                      theta_array

        # Integrate over \int dk k^2 leaving K-dependent part and keeping the
        # overall factor in the n_deltaU2 function below
        return np.sum(integrand_k)
    
    
    def n_deltaU2(self, q, kF):
        """
        Evaluates integration over \int dK K^2 in fourth term in
        U n(q) U^\dagger ~ \delta U \delta U^\dagger.
        
        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        kF : float
            Fermi momentum [fm^-1].
            
        Returns
        -------
        output : float
            Momentum distribution from \delta U \delta U^\dagger term before
            integration over \int dr r^2.
            
        """
        
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

        # Loop over K
        for iK, K in enumerate(K_array):
            
            integrand_K[iK] = self.n_deltaU2_K(q, K, kF)
            
        # Attach integration measure
        integrand_K *= K_array**2 * K_weights
        
        # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2
        
        # Integrate over \int dK K^2
        return deltaU2_factor * np.sum(integrand_K)


    def n_lambda(self, q_array, r_array):
        """
        Deuteron single-nucleon momentum distribution under LDA.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1]. These are not relative momenta.
        r_array : 1-D ndarray
            Coordinates array [fm].

        Returns
        -------
        n_lambda : 2-D ndarray
            Array of contributions to the deuteron momentum distribution
            ordered according to total, 1, \delta U, and \delta U^2 [fm^3].

        """
        
        # Number of columns for output array (total, 1, \delta U, \delta U^2)
        axes = 4
        
        # First compute deuteron density in coordinate space
        
        # Get momentum mesh from before
        k_array, k_weights, ntot_k = self.k_array, self.k_weights, self.ntot

        # Transform to coordinate space
        hank_trans_3S1 = hankel_transformation_k2r('3S1', k_array, k_weights,
                                                   r_array)
        hank_trans_3D1 = hankel_transformation_k2r('3D1', k_array, k_weights,
                                                   r_array)
    
        # Get 3S1 and 3D1 waves [fm^-3/2] in coordinate-space
        psi_r_3S1 = hank_trans_3S1 @ self.psi_k[:ntot_k]
        psi_r_3D1 = hank_trans_3D1 @ self.psi_k[ntot_k:]
    
        # Calculate the deuteron nucleonic density [fm^-3]
        rho_array = psi_r_3S1**2 + psi_r_3D1**2
        
        # Number of r_array points
        mtot = len(r_array)
        r2_grid, _ = np.meshgrid(r_array**2, np.zeros(axes), indexing='ij')
        dr = 0.1 # Spacing between r-points
            
        # Length of q_array
        ntot_q = len(q_array)
        
        # Evaluate n_\lambda(q, kF) contributions at each point in q_array
        # and rho_array
        n_lambda = np.zeros( (ntot_q, axes) )
    
        # kF values for each point in r_array
        kF_array = ( 3*np.pi**2 * rho_array )**(1/3)
        
        # Loop over q
        for i, q in enumerate(q_array):
    
            # n_lambda contributions before integrating over r
            n_lambda_r = np.zeros( (mtot, axes) )

            # Loop over kF values
            for j, kF in enumerate(kF_array):

                term_1 = self.n_1(q, kF)
                term_deltaU = self.n_deltaU(q, kF)
                term_deltaU2 = self.n_deltaU2(q, kF)
                total = term_1 + term_deltaU + term_deltaU2
                
                n_lambda_r[j, :] = total, term_1, term_deltaU, term_deltaU2
                    
            # Deuteron normalization is \int dr r^2 \rho_d(r) = 1 without
            # any factor of 4*\pi
            # For nucleus A, normalization would be
            #     4*\pi \int dr r^2 \rho_A(r) = Z or N depending on nucleon
            
            # Integrate over \int dr r^2 for each contribution
            n_lambda[i, :] = dr * np.sum(r2_grid * n_lambda_r, axis=0)
            
        # Return (ntot, 4) array of contributions to n_\lambda^d(q)
        return n_lambda
    
    
    def write_files(self):
        """
        Write deuteron momentum distribution files for interpolation purposes.
        Split things into total, 1, \delta U, and \delta U^2 contributions.

        """
        
        # Directory for distributions data
        data_directory = 'Data/dmd/kvnn_%d' % self.kvnn
        
        # Create file name
        file_name = 'deuteron_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = ff.replace_periods(file_name) + '.dat'
        
        # Use r_array in accordance with densities code from densities.py
        r_array = np.linspace(0.1, 20.0, 200)

        # Calculate n_\lambda^d(k) for each k in k_array
        n_array = self.n_lambda(self.k_array, r_array)
    
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
            
            
    def n_lambda_interp(self):
        """
        Interpolate the deuteron momentum distribution for the specified file.

        Returns
        -------
        output : tuple
            Tuple of functions that depend only on momentum q [fm^-1] where
            each function corresponds to contributions to n_\lambda(q): total,
            1, \delta U, and \delta U^2.

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
    
    
    def n_lambda_exact(self, q, contributions='total'):
        """
        Computes the exact deuteron momentum distribution using the wave
        function from the input Hamiltonian.

        Parameters
        ----------
        q : float
            Momentum value [fm^-1].
        contributions : str, optional
            Option to return different contributions to the momentum
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
    

if __name__ == '__main__':
    
    
    # # --- Compare LDA momentum distribution to wave function result --- #
    
    # import matplotlib.pyplot as plt
    # import time
    
    # channel = '3S1'
    
    # kvnn = 6
    # lamb = 1.35
    # kmax, kmid, ntot = 15.0, 3.0, 120
    
    # # Load momentum and weights
    # q_array, q_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
    # factor_array = q_array**2 * q_weights
    
    # # Load hamiltonian
    # H_matrix = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    
    # # Load exact wave function [unitless]
    # psi_exact_unitless = ob.wave_function(H_matrix)
    # psi_squared_exact = ( psi_exact_unitless[:ntot]**2 + \
    #                       psi_exact_unitless[ntot:]**2 ) / factor_array
        
    # # Calculate using LDA
    # t0 = time.time()
    # dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot, False)
    # r_array = np.linspace(0.1, 20.0, 200)
    # n_d_array_full = dmd.n_lambda(q_array, r_array)
    # t1 = time.time()
    
    # # Print minutes elapsed
    # mins = (t1-t0)/60
    # print('%.5f minutes elapsed.' % mins)
    
    # # Get just the total
    # n_d_array = n_d_array_full[:, 0]
    
    # # Normalization of wave function
    # norm = np.sum(factor_array * psi_squared_exact)
    # print('Normalization of exact: \dq q^2 n_d(q) = %.5f' % norm)
    
    # # Normalization of LDA deuteron momentum distribution
    # lda_factor = 4*np.pi * 1/(2*np.pi)**3
    # norm_lda = lda_factor * np.sum(factor_array * n_d_array)
    # print('Normalization of LDA: 4\pi/(2\pi)^3 \dq q^2 <n_d^N(q)> = %.5f'
    #       % norm_lda) 
    
    # # Plot pair momentum distributions
    # plt.semilogy(q_array, psi_squared_exact, label='AV18 wave function')
    # plt.semilogy(q_array, n_d_array * lda_factor, label='LDA')
    # plt.xlim( (0.0, 5.0) )
    # plt.ylim( (1e-5, 1e4) )
    # plt.xlabel(r'$q$' + ' [fm' + r'$^{-1}$' + ']')
    # plt.ylabel(r'$n_d(q)$' + ' [fm' + r'$^3$' + ']')
    # plt.legend(loc=0)
    
    
    # --- Generate all data for deuteron momentum distributions --- #
    # Currently this takes about ... hours to run
    
    import time
    
    # Potentials
    kvnns_list = [6, 222, 224]
    
    # SRG \lambda values
    lambdas_list = [6.0, 3.0, 2.0, 1.5, 1.35]
    
    # Momentum mesh details
    kmax, kmid, ntot = 15.0, 3.0, 120 # Default
    
    # Loop over every case generating data files
    for kvnn in kvnns_list:
        
        t0_k = time.time()

        for lamb in lambdas_list:
                
            t0_l = time.time()

            # Initialize class
            dmd = deuteron_momentum_distributions(kvnn, lamb, kmax, kmid, ntot,
                                                  interp=False)

            # Write deuteron file
            dmd.write_files()

            # Time for each \lambda to run   
            t1_l = time.time()
            mins_l = (t1_l-t0_l)/60
            print( '\tDone with \lambda=%s after %.5f minutes.' % (
                    ff.convert_number_to_string(lamb), mins_l) )
        
        # Time for each potential to run
        t1_k = time.time()
        mins_k = (t1_k-t0_k)/60
        print( 'Done with kvnn=%d after %.5f mins.\n' % (kvnn, mins_k) )