#!/usr/bin/env python3

"""
File: snmd.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 17, 2022

The SingleNucleon class calculates single-nucleon momentum distributions (snmd)
SRG-evolving the operator assuming the evolved wave function is given by HF 
treated in LDA. This class is a sub-class of the MomentumDistribution class
from momentum_distributions.py.

Last update: April 20, 2022

"""

# To-do: Make sure R_array parameter is described correctly.

# Python imports
import numpy as np
from numpy.polynomial.legendre import leggauss

# Imports from A.T. codes
from densities import load_density
from modules.integration import gaussian_quadrature_mesh
from modules.labels import replace_periods
from momentum_distributions import MomentumDistribution


class SingleNucleon(MomentumDistribution):
    
    def __init__(
            self, kvnn, kmax, kmid, ntot, channels, generator, lamb,
            lambda_initial=None, kvnn_inv=None, delta_lambda=None):
        """
        Sets \delta U(k,k') and \delta U^2(k,k') functions as attributes.

        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        kmax : float
            Maximum value in the momentum mesh [fm^-1].
        kmid : float
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int
            Number of momentum points in mesh.
        channels : tuple
            Partial wave channels to include in the calculation.
        generator : str
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float
            SRG evolution parameter \lambda [fm^-1].
        lambda_initial : float, optional
            SRG evolution parameter \lambda for initial Hamiltonian [fm^-1].
            This allows one to use an SRG-evolved potential as the starting
            point.
        kvnn_inv : int, optional
            This number specifies a potential for which inverse-SRG
            transformations will be applied to the initial Hamiltonian
            H_initial = U_{kvnn_inv}^{\dagger} H_kvnn U_{kvnn_inv},
            where the transformations are evaluated at \delta \lambda.
        delta_lambda : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and delta_lambda must be specified
            for this to run.

        """

        super().__init__(kvnn, kmax, kmid, ntot)
        
        # Set-up calculation
        self.save_deltaU_funcs(channels, generator, lamb, lambda_initial,
                               kvnn_inv, delta_lambda)
        
        # Set instance attributes for saving files
        self.channels = channels
        self.generator = generator
        self.lamb = lamb
        self.lambda_initial = lambda_initial
        self.kvnn_inv = kvnn_inv
        self.delta_lambda = delta_lambda

    def get_theta_I(self, q_grid, kF1_grid):
        """
        Evaluates the Heaviside step function \theta(kF1(R) - q) which appears
        in the I term.

        Parameters
        ----------
        q_grid : 2-D ndarray
            Meshgrid of momentum values [fm^-1].
        kF1_grid : 2-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau with respect to R.

        Returns
        -------
        theta_grid : 2-D ndarray
            Heaviside step function [unitless] evaluated on meshgrids of q
            and kF1(R).

        """
        
        theta_grid = np.zeros_like(q_grid)
        
        # Gives 1 if q < kF1(R)
        theta_grid[q_grid < kF1_grid] = 1
        
        # This is a (ntot_q, ntot_R) shape array
        return theta_grid
    
    def angle_avg_deltaU(self, q_grid, kF1_grid, kF2_grid, k_grid):
        """
        Evaluate angle-average involving Heaviside step function in the
        \delta U term:
            \int dx/2 \theta(kF2(R) - |q-2k|),
        where x is the angle between q_vector and k_vector.

        Parameters
        ----------
        q_grid : 3-D ndarray
            Meshgrid of momentum values [fm^-1].
        kF1_grid : 3-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau with respect to R.
        kF2_grid : 3-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau' with respect to R.
        k_grid : 3-D ndarray
            Meshgrid of relative momentum values [fm^-1].

        Returns
        -------
        angle_avg_grid : 3-D ndarray
            Angle-averaging integral [unitless] that weights the \delta U term
            evaluated on meshgrids of q, kF1(R), kF2(R), and k.
            
        Notes
        -----
        Not sure why the cases had to be computed in reverse order. Does that
        mean there is an overlap of truth values in case 3 and case 1 when
        there shouldn't be?

        """
        
        angle_avg_grid = np.zeros_like(q_grid)
        
        # Evaluate each boolean case and use these to fill in the meshgrid
        
        # This condition applies to each case: q < kF1
        common_constraint = q_grid < kF1_grid
        
        # Case 3: q-kF < 2k < kF+q
        mask_3 = (common_constraint * (q_grid-kF2_grid < 2*k_grid)
                  * (2*k_grid < kF2_grid+q_grid))
        angle_avg_grid[mask_3] = ((kF2_grid**2 - (q_grid-2*k_grid)**2)
                                  / (8*k_grid*q_grid))[mask_3]
            
        # Case 2: q < kF and kF-q <= 2k < kF+q
        mask_2 = (common_constraint * (q_grid < kF2_grid) *
                  (kF2_grid-q_grid <= 2*k_grid) * (2*k_grid < kF2_grid+q_grid))
        angle_avg_grid[mask_2] = ((kF2_grid**2 - (q_grid-2*k_grid)**2)
                                  / (8*k_grid*q_grid))[mask_2]
        
        # Case 1: q < kF and 2k < kF-q
        mask_1 = (common_constraint * (q_grid < kF2_grid)
                  * (2*k_grid < kF2_grid-q_grid))
        angle_avg_grid[mask_1] = 1

        # This is a (ntot_q, ntot_R, ntot_k) shape array
        return angle_avg_grid

    def angle_avg_deltaU2(self, kF1_grid, kF2_grid, K_grid, k_grid):
        """
        Evaluate angle-average involving Heaviside step functions in the
        \delta U \delta U^\dagger term:
            \int dz/2 \theta(kF1(R) - |K/2+k|) \theta(kF2(R) - |K/2-k|),
        where z is the angle between K_vector and k_vector.

        Parameters
        ----------
        kF1_grid : 3-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau with respect to R.
        kF2_grid : 3-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau' with respect to R.
        K_grid : 4-D ndarray
            Meshgrid of C.o.M. momentum values [fm^-1].
        k_grid : 4-D ndarray
            Meshgrid of relative momentum values [fm^-1].

        Returns
        -------
        angle_avg_grid : 4-D ndarray
            Angle-averaging integral [unitless] that weights the \delta U^2
            term evaluated on meshgrids of kF1(R), kF2(R), K, and k. Note, the
            first dimension of this meshgrid corresponds to q even though it
            doesn't depend on q. This is to match the shape of the \delta U
            matrix elements.

        """
        
        angle_avg_grid = np.zeros_like(kF1_grid)
        
        # Evaluate each boolean case and use these to fill in the meshgrid
        
        # Case 4: 2k+K < 2kF1 and -4 < (4k^2-4kF2^2+K^2)/(kK) < 4
        mask_4 = ((2*k_grid+K_grid <= 2*kF1_grid)
                  * (-4 < (4*k_grid**2-4*kF2_grid**2+K_grid**2)
                     / (k_grid*K_grid))
                  * ((4*k_grid**2-4*kF2_grid**2+K_grid**2)
                     /(k_grid*K_grid) <= 4))
        angle_avg_grid[mask_4] = ((4*kF2_grid**2-(K_grid-2*k_grid)**2 )
                                  / (8*k_grid*K_grid))[mask_4]
            
        # Case 3: 2k+K < 2kF2 and -4 < (4k^2 - 4kF1^2 + K^2)/(kK) < 4
        mask_3 = ((2*k_grid+K_grid <= 2*kF2_grid)
                  * (-4 < (4*k_grid**2-4*kF1_grid**2+K_grid**2)
                     / (k_grid*K_grid))
                  * ((4*k_grid**2-4*kF1_grid**2+K_grid**2)
                     / (k_grid*K_grid) <= 4 ))
        angle_avg_grid[mask_3] = ((4*kF1_grid**2-(K_grid-2*k_grid)**2)
                                  / (8*k_grid*K_grid))[mask_3]
            
        # Case 2: 2k+K > 2kF1 and 2k+K > 2kF2 and 4k^2+K^2 < 2(kF1^2+kF2^2)
        mask_2 = ((2*k_grid+K_grid > 2*kF1_grid)
                  * (2*k_grid+K_grid > 2*kF2_grid)
                  * (4*k_grid**2+K_grid**2 <= 2*(kF1_grid**2+kF2_grid**2))) 
        angle_avg_grid[mask_2] = ((2*(kF1_grid**2+kF2_grid**2) 
                                   - 4*k_grid**2-K_grid**2)
                                  / (4*k_grid*K_grid))[mask_2]

        # Case 1: 2k+K < 2kF1 and 2k+K < 2kF2
        mask_1 = ((2*k_grid+K_grid <= 2*kF1_grid)
                  * (2*k_grid+K_grid <= 2*kF2_grid))
        angle_avg_grid[mask_1] = 1
        
        # This is a (ntot_q, ntot_R, ntot_K, ntot_k) shape array
        return angle_avg_grid
        
    def get_I_term(self, q_array, R_array, dR, kF1_array):
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
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from I term as a function of q.

        """
        
        # Initialize 2-D meshgrids (q, R)
        q_grid, R_grid = np.meshgrid(q_array, R_array, indexing='ij')
        
        # Get 2-D kF1(R) meshgrid
        _, kF1_grid = np.meshgrid(q_array, kF1_array, indexing='ij')
        
        # Evaluate the Heaviside step function in I term
        theta_grid = self.get_theta_I(q_grid, kF1_grid)
        
        # Calculate R integrand (ntot_q, ntot_R)
        integrand_R = theta_grid * R_grid**2 * dR

        # Integrate over R leaving a (ntot_q, 1) shape array
        # Factor of 2 is overall factor for summation over spin projections
        return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def get_deltaU_term(self, q_array, R_array, dR, kF1_array, kF2_array):
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
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U term as a function of q.

        """
        
        # Set number of k points for integration over k
        ntot_k = 40  # Typical integration goes up to kF ~ 1.3 fm^-1
        
        # Initialize 3-D meshgrids (q, R, k) with k values equal to 0
        q_grid, R_grid, k_grid = np.meshgrid(
            q_array, R_array, np.zeros(ntot_k), indexing='ij')
        
        # Get 3-D kF1(R) and kF2(R) meshgrids and initialize k weights mesh
        _, kF1_grid, dk_grid = np.meshgrid(
            q_array, kF1_array, np.zeros(ntot_k), indexing='ij')
        _, kF2_grid, _ = np.meshgrid(
            q_array, kF2_array, np.zeros(ntot_k), indexing='ij')

        # Loop over q and kF2 to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iq, q in enumerate(q_array):
            for ikF2, kF2 in enumerate(kF2_array):

                # Create integration meshgrid k_array up to (kF2 + q)/2 which
                # corresponds to the upper limit of \theta(kF2(R) - |q-2k|)
                k_max = (kF2 + q)/2
 
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max, ntot_k)
                
                # Fill in k_grid and dk_grid given the specific k_array
                k_grid[iq, ikF2, :] = k_array
                dk_grid[iq, ikF2, :] = k_weights
        
        # Evaluate angle-average over Heaviside step functions in \delta U 
        # term for \tau and \tau'
        angle_avg_pp_grid = self.angle_avg_deltaU(q_grid, kF1_grid, kF1_grid,
                                                  k_grid)
        angle_avg_pn_grid = self.angle_avg_deltaU(q_grid, kF1_grid, kF2_grid,
                                                  k_grid)
        
        # Evaluate \delta U(k,k) for \tau and \tau'
        deltaU_pp_grid = self.deltaU_pp_func(k_grid)
        deltaU_pn_grid = self.deltaU_pn_func(k_grid)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # Factor of 2 from \delta U + \delta U^\dagger
        # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        deltaU_factor = 2 * 8 * 2/np.pi * 2
        
        # Calculate the k integrand where we split terms according to pp and
        # pn (ntot_q, ntot_R, ntot_k)
        integrand_k = (deltaU_factor * k_grid**2 * dk_grid * R_grid**2 * dR
                       * (deltaU_pp_grid*angle_avg_pp_grid
                          + deltaU_pn_grid*angle_avg_pn_grid))
        
        # Integrate over k leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R leaving a (ntot_q, 1) shape array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def get_deltaU2_term(self, q_array, R_array, dR, kF1_array, kF2_array):
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
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U \delta U^\dagger term as
            a function of q.

        """
        
        ntot_q = len(q_array)
        # Set number of K and k points for integration over K and k
        ntot_K = 40
        ntot_k = 40
        
        # y = cos(\theta) angles for averaging integration over angle between
        # K/2 and k
        ntot_y = 7
        y_array, y_weights = leggauss(ntot_y)
        
        # Initialize 4-D meshgrids (q, R, K, k) with k and K equal to 0
        q_grid, R_grid, K_grid, k_grid = np.meshgrid(
            q_array, R_array, np.zeros(ntot_K), np.zeros(ntot_k),
            indexing='ij')
        
        # Get 4-D kF1(R) and kF2(R) meshgrids and initialize K and k weights
        # meshgrids
        _, kF1_grid, dK_grid, dk_grid = np.meshgrid(
            q_array, kF1_array, np.zeros(ntot_K), np.zeros(ntot_k),
            indexing='ij')
        _, kF2_grid, _, _ = np.meshgrid(q_array, kF2_array, np.zeros(ntot_K),
                                        np.zeros(ntot_k), indexing='ij')

        # Loop over q, R, and k to find limits of K integration and then
        # create K integration mesh using Gaussian quadrature
        for iR, R in enumerate(R_array):
                
            kF1, kF2 = kF1_array[iR], kF2_array[iR]
    
            # K integration goes from 0 to kF1+kF2
            K_max = kF1 + kF2
                
            # Get Gaussian quadrature mesh for K integration
            K_array, K_weights = gaussian_quadrature_mesh(K_max, ntot_K)
              
            # Loop over remaining variables and fill in 4-D array
            for iq in range(ntot_q):
                for ik in range(ntot_k):
                    
                    # Fill in K_grid and dK_grid given the specific K_array
                    K_grid[iq, iR, :, ik] = K_array
                    dK_grid[iq, iR, :, ik] = K_weights

        # Loop over q, R, and K to find limits of k integration and then
        # create k_array using Gaussian quadrature
        for iR, R in enumerate(R_array):
        
            kF1, kF2 = kF1_array[iR], kF2_array[iR]
        
            # Get minimum kF value
            kF_min = min(kF1, kF2)
            
            # K_array only depends on R, so loop over K_grid[0, iR, :, 0]
            for iK, K in enumerate(K_grid[0, iR, :, 0]):
                
                # Lower limit of k integration
                k_min = max(K/2 - kF_min, 0)
                
                # Upper limit of k integration
                if K**2/4 < (kF1**2 + kF2**2)/2:
                    k_max = min(np.sqrt((kF1**2+kF2**2)/2-K**2/4), kF_min+K/2)
                else:
                    k_max = kF_min + K/2

                # Get Gaussian quadrature mesh for k integration
                k_array, k_weights = gaussian_quadrature_mesh(
                    k_max, ntot_k, xmin=k_min)
                
                # Loop over remaining variables and fill in 4-D meshgrid
                for iq in range(ntot_q):
                    
                    # Fill in k_grid and dk_grid given the specific k_array
                    k_grid[iq, iR, iK, :] = k_array
                    dk_grid[iq, iR, iK, :] = k_weights

        # Evaluate angle-average over Heaviside step functions in \delta U^2
        # term for \tau and \tau'
        angle_avg_pp_grid = self.angle_avg_deltaU2(kF1_grid, kF1_grid, K_grid,
                                                   k_grid)
        angle_avg_pn_grid = self.angle_avg_deltaU2(kF1_grid, kF2_grid, K_grid,
                                                   k_grid)

        # Set-up 4-D \delta U \delta U^\dagger(k, |q-K/2|) for \tau and \tau'
        deltaU2_pp_grid = np.zeros_like(angle_avg_pp_grid)
        deltaU2_pn_grid = np.zeros_like(angle_avg_pn_grid)
        
        # Integrate over y (angle between q_vector and K_vector)
        for y, dy in zip(y_array, y_weights):
            
            # Evaluate |q-K/2| meshgrid
            q_K_grid = np.sqrt(q_grid**2 + K_grid**2/4 - q_grid*K_grid*y)
            
            deltaU2_pp_grid += self.deltaU2_pp_func.ev(k_grid, q_K_grid) * dy/2
            deltaU2_pn_grid += self.deltaU2_pn_func.ev(k_grid, q_K_grid) * dy/2
        
        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2

        # Calculate the k integrand where we split terms according to pp and
        # pn leaving a (ntot_q, ntot_R, ntot_K, ntot_k) shape array
        integrand_k = (deltaU2_factor * k_grid**2 * dk_grid * K_grid**2 
                       * dK_grid * R_grid**2 * dR
                       * (deltaU2_pp_grid*angle_avg_pp_grid
                          +deltaU2_pn_grid*angle_avg_pn_grid))
                      
        # Integrate over k leaving K integrand (ntot_q, ntot_R, ntot_K)
        integrand_K = np.sum(integrand_k, axis=-1)
        
        # Integrate over K leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_K, axis=-1)
        
        # Integrate over R leaving a (ntot_q, 1) shape array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def compute_momentum_distribution(
            self, q_array, nucleon, nucleus_name, Z, N, density='Gogny',
            save=False):
        """
        Single-nucleon momentum distribution for a specified nucleus.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        nucleon : str
            Specify 'proton' or 'neutron'.
        nucleus_name : str
            Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
        density : str, optional
            Name of nucleonic density (e.g., 'SLy4', 'Gogny').
        save : bool, optional
            Option to save data in data/momentum_distributions directory.

        Returns
        -------
        n_total : 1-D ndarray
            Single-nucleon momentum distribution n^\tau(q) [fm^3].

        """
        
        # Load R values and nucleonic densities (the R_array's are the same)
        R_array, rho_p_array = load_density('proton', nucleus_name, Z, N,
                                            density)
        R_array, rho_n_array = load_density('neutron', nucleus_name, Z, N,
                                            density)
        dR = R_array[2] - R_array[1]  # Assuming linear spacing
        
        # Evaluate kF values at each point in R_array
        if nucleon == 'proton':
            kF1_array = (3*np.pi**2 * rho_p_array)**(1/3)
            kF2_array = (3*np.pi**2 * rho_n_array)**(1/3)
        else:
            kF1_array = (3*np.pi**2 * rho_n_array)**(1/3)
            kF2_array = (3*np.pi**2 * rho_p_array)**(1/3)
        
        # Get each contribution with respect to q (ntot_q, 1)
        n_I = self.get_I_term(q_array, R_array, dR, kF1_array)
        n_deltaU = self.get_deltaU_term(q_array, R_array, dR, kF1_array,
                                        kF2_array)
        n_deltaU2 = self.get_deltaU2_term(q_array, R_array, dR, kF1_array,
                                          kF2_array)
        
        n_total = n_I + n_deltaU + n_deltaU2
        
        if save:
            self.save_momentum_distribution(
                q_array, n_total, n_I, n_deltaU, n_deltaU2, nucleon, 
                nucleus_name, density)

        return n_total
    
    def save_momentum_distribution(
            self, q_array, n_total, n_I, n_deltaU, n_deltaU2, nucleon,
            nucleus_name, density):
        """
        Saves momentum distribution in data/momentum_distributions. Here the
        specifications of the potential, SRG evolution, nucleon, and density
        dictate the file name.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        n_total : 1-D ndarray
            Total momentum distribution [fm^3].
        n_I : 1-D ndarray
            Contribution to the momentum distribution from the I term [fm^3].
        n_deltaU : 1-D ndarray
            Contribution to the momentum distribution from the \delta U term
            [fm^3].
        n_deltaU2 : 1-D ndarray
            Contribution to the momentum distribution from the \delta U^2 term
            [fm^3].
        nucleon : str
            Specify 'proton' or 'neutron'.
        nucleus_name : str
            Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
        density : str
            Name of nucleonic density (e.g., 'SLy4', 'Gogny').

        """
        
        # Directory for distributions data
        data_directory = f'../data/momentum_distributions/{nucleus_name}/'
        
        # Create file name
        file_name = (f'n_{nucleon}_{density}_kvnn_{self.kvnn}_kmax_{self.kmax}'
                     f'_kmid_{self.kmid}_ntot_{self.ntot}')
        
        for channel in self.channels:
            file_name += f'_{channel}'
        
        file_name += f'_{self.generator}'
        if self.generator == 'Block-diag':
            file_name += f'_LambdaBD_{self.lamb}'
        else:
            file_name += f'_lambda_{self.lamb}'
        
        if self.lambda_initial != None:
            file_name += f'_lambda_initial_{self.lambda_initial}'
            
        if self.kvnn_inv != None:
            file_name += (f'_kvnn_inv_{self.kvnn_inv}'
                          f'_delta_lamb_{self.delta_lambda}')
        
        file_name = replace_periods(file_name) + '.dat'
        
        # Open file and write header where we allocate 18 centered spaces for
        # each label
        f = open(data_directory + file_name, 'w')
        header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}'.format(
            'q', 'total', '1', '\delta U', '\delta U^2')
        f.write(header + '\n')
    
        # Loop over momenta q and write each contribution
        for iq, q in enumerate(q_array):

            # Write to data file following the format from the header
            line = (f'{q:^18.6f}{n_total[iq]:^18.6e}{n_I[iq]:^18.6e}'
                    f'{n_deltaU[iq]:^18.6e}{n_deltaU2[iq]:^18.6e}')
            f.write('\n' + line)

        f.close()