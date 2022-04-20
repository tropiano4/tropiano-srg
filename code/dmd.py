#!/usr/bin/env python3

"""
File: dmd.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 17, 2022

The Deuteron class calculates deuteron momentum distributions (dmd) SRG-
evolving the operator assuming the evolved wave function is given by HF
treated in LDA. This class is a sub-class of the MomentumDistribution class 
from momentum_distributions.py.

Notes on normalizations:
  1. The deuteron wave function describing relative position or momentum is
      normalized according to
        \int dr r^2 (|\psi_{3S1}(r)|^2 + |\psi_{3D1}|^2) = 1,
        2/\pi * \int dk k^2 (|\psi_{3S1}(k)|^2 + |\psi_{3D1}(k)|^2) = 1.
  2. Under HF+LDA, we adopt the normalization
        4\pi / (2\pi)^3 \int dk k^2 < n_d(k) > = 1,
     where angled-brackets indicate nuclear-averaging (integration over R).

Last update: April 20, 2022

"""

# To-do: Make sure R_array parameter is described correctly.
# To-do: There are several things that overlap (if not slightly) with the
# single-nucleon momentum distribution code. This suggests making Deuteron
# inherit SingleNucleon.
# To-do: Think about integration over K. Does this even make sense?

# Python imports
import numpy as np
from numpy.polynomial.legendre import leggauss

# Imports from A.T. codes
from densities import load_density
from modules.fourier_transform import hankel_transformation_k2r
from modules.integration import gaussian_quadrature_mesh
from modules.labels import replace_periods
from momentum_distributions import MomentumDistribution


class Deuteron(MomentumDistribution):

    def __init__(
            self, kvnn, kmax, kmid, ntot, generator, lamb,
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
        self.save_deuteron_deltaU_funcs(generator, lamb, lambda_initial,
                                        kvnn_inv, delta_lambda)
        
        # Set instance attributes for saving files
        self.generator = generator
        self.lamb = lamb
        self.lambda_initial = lambda_initial
        self.kvnn_inv = kvnn_inv
        self.delta_lambda = delta_lambda
        
    def get_theta_I(self, q_grid, kF_grid):
        """
        Evaluates the Heaviside step function \theta(kF(R) - q) which appears
        in the I term.

        Parameters
        ----------
        q_grid : 2-D ndarray
            Meshgrid of momentum values [fm^-1].
        kF_grid : 2-D ndarray
            Meshgrid of deuteron Fermi momenta [fm^-1] with respect to R.

        Returns
        -------
        theta_grid : 2-D ndarray
            Heaviside step function [unitless] evaluated on meshgrids of q
            and kF(R).

        """
        
        theta_grid = np.zeros_like(q_grid)
        
        # Gives 1 if q < kF1(R)
        theta_grid[q_grid < kF_grid] = 1
        
        # This is a (ntot_q, ntot_R) shape array
        return theta_grid
    
    def angle_avg_deltaU(self, q_grid, kF_grid, k_grid):
        """
        Evaluate angle-average involving Heaviside step function in the
        \delta U term:
            \int dx/2 \theta(kF(R) - |q-2k|),
        where x is the angle between q_vector and k_vector.

        Parameters
        ----------
        q_grid : 3-D ndarray
            Meshgrid of momentum values [fm^-1].
        kF_grid : 2-D ndarray
            Meshgrid of deuteron Fermi momenta [fm^-1] with respect to R.
        k_grid : 3-D ndarray
            Meshgrid of relative momentum values [fm^-1].

        Returns
        -------
        angle_avg_grid : 3-D ndarray
            Angle-averaging integral [unitless] that weights the \delta U term
            evaluated on meshgrids of q, kF(R), and k.
            
        Notes
        -----
        Not sure why the cases had to be computed in reverse order. Does that
        mean there is an overlap of truth values in case 3 and case 1 when
        there shouldn't be?

        """
        
        angle_avg_grid = np.zeros_like(q_grid)
        
        # Evaluate each boolean case and use these to fill in the meshgrid
        
        # This condition applies to each case: q < kF
        common_constraint = q_grid < kF_grid
        
        # Case 3: q-kF < 2k < kF+q
        mask_3 = (common_constraint * (q_grid-kF_grid < 2*k_grid)
                  * (2*k_grid < kF_grid+q_grid))
        angle_avg_grid[mask_3] = ((kF_grid**2 - (q_grid-2*k_grid)**2)
                                  / (8*k_grid*q_grid))[mask_3]
            
        # Case 2: q < kF and kF-q <= 2k < kF+q
        mask_2 = (common_constraint * (q_grid < kF_grid) *
                  (kF_grid-q_grid <= 2*k_grid) * (2*k_grid < kF_grid+q_grid))
        angle_avg_grid[mask_2] = ((kF_grid**2 - (q_grid-2*k_grid)**2)
                                  / (8*k_grid*q_grid))[mask_2]
        
        # Case 1: q < kF and 2k < kF-q
        mask_1 = (common_constraint * (q_grid < kF_grid)
                  * (2*k_grid < kF_grid-q_grid))
        angle_avg_grid[mask_1] = 1

        # This is a (ntot_q, ntot_R, ntot_k) shape array
        return angle_avg_grid

    def angle_avg_deltaU2(self, kF_grid, K_grid, k_grid):
        """
        Evaluate angle-average involving Heaviside step functions in the
        \delta U \delta U^\dagger term:
            \int dz/2 \theta(kF(R) - |K/2+k|) \theta(kF(R) - |K/2-k|),
        where z is the angle between K_vector and k_vector.

        Parameters
        ----------
        kF_grid : 2-D ndarray
            Meshgrid of deuteron Fermi momenta [fm^-1] with respect to R.
        K_grid : 4-D ndarray
            Meshgrid of C.o.M. momentum values [fm^-1].
        k_grid : 4-D ndarray
            Meshgrid of relative momentum values [fm^-1].

        Returns
        -------
        angle_avg_grid : 4-D ndarray
            Angle-averaging integral [unitless] that weights the \delta U^2
            term evaluated on meshgrids of kF(R), K, and k. Note, the first
            dimension of this meshgrid corresponds to q even though it doesn't
            depend on q. This is to match the shape of the \delta U matrix
            elements.

        """
        
        angle_avg_grid = np.zeros_like(kF_grid)
        
        # Evaluate each boolean case and use these to fill in the meshgrid
        
        # Case 2: 2k+K > 2kF and 4k^2+K^2 < 4kF^2
        mask_2 = ((2*k_grid+K_grid > 2*kF_grid)
                  * (4*k_grid**2+K_grid**2 <= 4*kF_grid**2))
        angle_avg_grid[mask_2] = ((4*kF_grid**2-4*k_grid**2-K_grid**2)
                                  / (4*k_grid*K_grid))[mask_2]

        # Case 1: 2k+K < 2kF
        mask_1 = (2*k_grid+K_grid <= 2*kF_grid)
        angle_avg_grid[mask_1] = 1
        
        # This is a (ntot_q, ntot_R, ntot_K, ntot_k) shape array
        return angle_avg_grid

    def get_I_term(self, q_array, R_array, dR, kF_array):
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
        q_grid, R_grid = np.meshgrid(q_array, R_array, indexing='ij')
        
        # Get 2-D kF(R) meshgrid
        _, kF_grid = np.meshgrid(q_array, kF_array, indexing='ij')
        
        # Evaluate the Heaviside step function in I term
        theta_grid = self.get_theta_I(q_grid, kF_grid)
        
        # Calculate R integrand (ntot_q, ntot_R)
        integrand_R = theta_grid * R_grid**2 * dR

        # Integrate over R leaving a (ntot_q, 1) shape array
        # Factor of 2 is overall factor for summation over spin projections
        return 2 * np.sum(integrand_R, axis=-1)
    
    def get_deltaU_term(self, q_array, R_array, dR, kF_array):
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
        ntot_k = 40  # Typical integration goes up to kF ~ 1.3 fm^-1
        
        # Initialize 3-D meshgrids (q, R, k) with k values equal to 0
        q_grid, R_grid, k_grid = np.meshgrid(q_array, R_array,
                                             np.zeros(ntot_k), indexing='ij')
        
        # Get 3-D kF(R) meshgrid and initialize k weights mesh
        _, kF_grid, dk_grid = np.meshgrid(q_array, kF_array, np.zeros(ntot_k),
                                          indexing='ij')

        # Loop over q and kF to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iq, q in enumerate(q_array):
            for ikF, kF in enumerate(kF_array):

                # Create integration meshgrid k_array up to (kF + q)/2 which
                # corresponds to the upper limit of \theta(kF(R) - |q-2k|)
                k_max = (kF + q)/2
 
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max, ntot_k)
                
                # Fill in k_grid and dk_grid given the specific k_array
                k_grid[iq, ikF, :] = k_array
                dk_grid[iq, ikF, :] = k_weights

        # Evaluate angle-average over Heaviside step functions in \delta U term
        angle_avg_grid = self.angle_avg_deltaU(q_grid, kF_grid, k_grid)
        
        # Evaluate \delta U(k,k)
        deltaU_grid = self.deltaU_func(k_grid)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # Factor of 2 from \delta U + \delta U^\dagger
        # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        deltaU_factor = 2 * 8 * 2/np.pi * 2
        
        # Calculate the k integrand leaving (ntot_q, ntot_R, ntot_k)
        integrand_k = (deltaU_factor * k_grid**2 * dk_grid * R_grid**2 * dR
                       * deltaU_grid * angle_avg_grid)
        
        # Integrate over k leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R leaving a (ntot_q, 1) shape array
        return np.sum(integrand_R, axis=-1)
    
    def get_deltaU2_term(self, q_array, R_array, dR, kF_array):
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
            Momentum distribution [fm^3] from \delta U \delta U^\dagger term 
            as a function of q.

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
        
        # Get 4-D kF(R) meshgrid and initialize K and k weights meshgrids
        _, kF_grid, dK_grid, dk_grid = np.meshgrid(
            q_array, kF_array, np.zeros(ntot_K), np.zeros(ntot_k),
            indexing='ij')

        # Loop over q, R, and k to find limits of K integration and then
        # create K integration mesh using Gaussian quadrature
        for iR, R in enumerate(R_array):
                
            kF = kF_array[iR]
    
            # K integration goes from 0 to 2*kF
            K_max = 2*kF
                
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
        
            kF = kF_array[iR]
            
            # K_array only depends on R, so loop over K_grid[0, iR, :, 0]
            for iK, K in enumerate(K_grid[0, iR, :, 0]):
                
                # Lower limit of k integration
                k_min = max(K/2 - kF, 0)
                
                # Upper limit of k integration
                if K**2/4 < kF**2:
                    k_max = min(np.sqrt(kF**2-K**2/4), kF+K/2)
                else:
                    k_max = kF + K/2

                # Get Gaussian quadrature mesh for k integration
                k_array, k_weights = gaussian_quadrature_mesh(k_max, ntot_k,
                                                              xmin=k_min)
                
                # Loop over remaining variables and fill in 4-D meshgrid
                for iq in range(ntot_q):
                    
                    # Fill in k_grid and dk_grid given the specific k_array
                    k_grid[iq, iR, iK, :] = k_array
                    dk_grid[iq, iR, iK, :] = k_weights
        
        # Evaluate angle-average over Heaviside step functions in \delta U^2
        # term
        angle_avg_grid = self.angle_avg_deltaU2(kF_grid, K_grid, k_grid)

        # Set-up 4-D \delta U \delta U^\dagger(k, |q-K/2|)
        deltaU2_grid = np.zeros_like(angle_avg_grid)
        
        # Integrate over y
        for y, dy in zip(y_array, y_weights):
            
            # Evaluate |q-K/2| meshgrid
            q_K_grid = np.sqrt(q_grid**2 + K_grid**2/4 - q_grid*K_grid*y)
            
            deltaU2_grid += self.deltaU2_func.ev(k_grid, q_K_grid) * dy/2
        
        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2

        # Calculate the k integrand (ntot_q, ntot_R, ntot_K, ntot_k)
        integrand_k = (deltaU2_factor * k_grid**2 * dk_grid * K_grid**2 
                       * dK_grid * R_grid**2 * dR * deltaU2_grid
                       * angle_avg_grid)
                      
        # Integrate over k leaving K integrand (ntot_q, ntot_R, ntot_K)
        integrand_K = np.sum(integrand_k, axis=-1)
        
        # Integrate over K leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_K, axis=-1)
        
        # Integrate over R leaving a (ntot_q, 1) shape array
        return np.sum(integrand_R, axis=-1)
    
    def compute_momentum_distribution(self, q_array, save=False):
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
            Name of nucleonic density (e.g., 'SLY4', 'Gogny').
        save : bool, optional
            Option to save data in data/momentum_distributions directory.

        Returns
        -------
        n_total : 1-D ndarray
            Deuteron momentum distribution n(q) [fm^3].

        """
        
        # Get momentum mesh of input Hamiltonian and number of points
        k_array, k_weights, ntot_k = self.k_array, self.k_weights, self.ntot
        
        # Use R values in accordance with densities code from densities.py
        # Nuclei arguments don't matter here - just getting R values
        R_array, _ = load_density('proton', 'C12', 6, 6)
        dR = R_array[2] - R_array[1]  # Assuming linear spacing
        
        # Transform wave function to coordinate space
        hank_trans_3S1 = hankel_transformation_k2r(0, k_array, k_weights,
                                                   R_array)
        hank_trans_3D1 = hankel_transformation_k2r(2, k_array, k_weights,
                                                   R_array)
    
        # Get 3S1 and 3D1 waves [fm^-3/2] in coordinate space
        psi_R_3S1 = hank_trans_3S1 @ self.psi_k[:ntot_k]
        psi_R_3D1 = hank_trans_3D1 @ self.psi_k[ntot_k:]
    
        # Calculate the deuteron density [fm^-3]
        rho_array = psi_R_3S1**2 + psi_R_3D1**2

        # Evaluate kF values at each point in R_array
        kF_array = (3*np.pi**2 * rho_array)**(1/3)
        
        # Get each contribution with respect to q (ntot_q, 1)
        n_I = self.get_I_term(q_array, R_array, dR, kF_array)
        n_deltaU = self.get_deltaU_term(q_array, R_array, dR, kF_array)
        n_deltaU2 = self.get_deltaU2_term(q_array, R_array, dR, kF_array)

        n_total = n_I + n_deltaU + n_deltaU2
        
        if save:
            self.save_momentum_distribution(q_array, n_total, n_I, n_deltaU,
                                            n_deltaU2)

        return n_total
    
    def save_momentum_distribution(
            self, q_array, n_total, n_I, n_deltaU, n_deltaU2):
        """
        Saves momentum distribution in data/momentum_distributions/H2. Here 
        the specifications of the potential and SRG evolution dictate the 
        file name.

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

        """
        
        # Directory for distributions data
        data_directory = '../data/momentum_distributions/H2/'
        
        # Create file name
        file_name = (f'n_kvnn_{self.kvnn}_kmax_{self.kmax}_kmid_{self.kmid}'
                     f'_ntot_{self.ntot}_{self.generator}')
        
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
        

# Run this script to compute and save momentum distributions
if __name__ == '__main__':
    
    # Specific imports
    import time
    from potentials import Potential
    
    # Default inputs
    kvnn = 6
    kmax, kmid, ntot = 15.0, 3.0, 120
    generator = 'Wegner'
    lamb = 1.35
    
    dmd = Deuteron(kvnn, kmax, kmid, ntot, generator, lamb)
    
    # Get momentum values (channel argument doesn't matter here)
    potential = Potential(kvnn, '1S0', kmax, kmid, ntot)
    q_array, _ = potential.load_mesh() 
    
    t0 = time.time()
        
    n_array = dmd.compute_momentum_distribution(q_array, save=True)
            
    t1 = time.time()
    mins = (t1-t0)/60
        
    print(f'Done after {mins:.2f} minutes.')


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
#                integrations. Saved old version as dmd_v1.py in Old_codes.
#   07/15/21 --- Switched interpolation functions interp1d and 
#                RectBivariateSpline from cubic to linear since \delta U^2(k,q)
#                was returning negative values.
#   02/15/22 --- Added optional lambda_init argument to class. This allows one
#                to use an SRG-evolved potential as the initial interaction.
#   02/24/22 --- Added optional kvnn_hard argument to class. This allows one
#                to SRG-evolve the initial potential back to a harder one
#                using transformations of the kvnn_hard potential.
#
#------------------------------------------------------------------------------


# import numpy as np
# from numpy.polynomial.legendre import leggauss
# from scipy.interpolate import interp1d, RectBivariateSpline
# # Scripts made by A.T.
# # from densities import load_density
# # from figures import figures_functions as ff
# # from misc.fourier_transform import hankel_transformation_k2r
# # from misc.integration import gaussian_quadrature_mesh
# # import observables as ob
# # from operators import momentum_projection_operator
# # from potentials.vsrg_macos import vnn
# # from srg.srg_unitary_transformation import SRG_unitary_transformation
# from .figure_labels import replace_periods
# from .fourier_transform import hankel_transformation_k2r
# from .integration import gaussian_quadrature_mesh
# from .momentum_projection_operator import momentum_projection_operator
# from ..vnn import Potential
# from .srg_transformation import get_transformation
# from .wave_function import wave_function


# def load_density(nucleus_name, nucleon, Z, N, edf='SLY4'):
#     """
#     Loads a nucleonic density for the given nucleus. Densities are normalized
#     according to
#         4*\pi \int_0^\infty dR R^2 \rho_A(R) = Z or N.
    
#     Parameters
#     ----------
#     nucleus_name : str
#         Specify the nucleus (e.g., 'O16', 'Ca40', etc.)
#     nucleon : str
#         Specify 'proton' or 'neutron'.
#     Z : int
#         Proton number of the nucleus.
#     N : int
#         Neutron number of the nucleus.
#     edf : str, optional
#         Name of EDF (e.g., 'SLY4').
        
#     Returns
#     -------
#     R_array : 1-D ndarray
#         C.o.M. coordinates [fm].
#     rho_array : 1-D ndarray
#         Nucleonic density as a function of R [# of nucleons / vol].
                                              
#     Notes
#     -----
#     Momentum distributions code compute intermediate integration arrays in 
#     relative k and C.o.M. K which rely on kF(R) values. These values can be
#     zero if the density \rho(R) = 0. We must replace zeros in \rho(R) with an
#     extremely small number so the codes run correctly to avoid zero division
#     errors. (This only happens for edf = 'AV18' densities.)
    
#     """


#     # Go to directory corresponding to specified nucleus and EDF
#     if edf == 'SLY4':
        
#         densities_directory = f'../../densities/HFBRAD_{edf}/{nucleus_name}/'
#         file_name = f'{nucleon}_{N:d}_{Z:d}.dens'
#         column_number = 1
        
#     elif edf == 'Gogny':
        
#         densities_directory = f'../../densities/{edf}/{nucleus_name}/'
#         file_name = 'DensityQP.dat'
#         if nucleon == 'proton':
#             column_number = 1
#         elif nucleon == 'neutron':
#             column_number = 2
    
#     # Technically it doesn't make sense to have a case edf == AV18 since
#     # AV18 does not use an EDF. It would also make more sense to call it VMC.
#     elif edf == 'AV18':
        
#         densities_directory = '../../densities/{edf}/'
#         file_name = '{nucleus_name}_densities_{N:d}_{Z:d}.txt'
        
#         # AV18 files either have single \rho column for N=Z nuclei or
#         # two columns for proton (1) and neutron (3)
#         if N == Z:
#             column_number = 1
#         else:
#             if nucleon == 'proton':
#                 column_number = 1
#             elif nucleon == 'neutron':
#                 column_number = 3 
        
#     # Load file
#     table = np.loadtxt(densities_directory + file_name)
    
#     R_array = table[:, 0]
#     rho_array = table[:, column_number]
    
#     # Avoiding zero division errors
#     zero_case = rho_array == 0
#     rho_array[zero_case] = 1e-30
    
#     return R_array, rho_array


# class deuteron_momentum_distributions(object):
    
    
#     def __init__(self, kvnn, lamb, kmax, kmid, ntot, interp=False,
#                  lambda_init=np.inf, kvnn_hard=0):
#         """
#         Evaluates and saves 3S1-3D1 matrix elements of \delta U and
#         \delta U^{\dagger} given the input potential and SRG \lambda.
        
#         Parameters
#         ----------
#         kvnn : int
#             This number specifies the potential.
#         lamb : float
#             SRG evolution parameter \lambda [fm^-1].
#         kmax : float
#             Maximum value in the momentum mesh [fm^-1].
#         kmid : float
#             Mid-point value in the momentum mesh [fm^-1].
#         ntot : int
#             Number of momentum points in mesh.
#         interp : bool, optional
#             Option to use interpolated n_\lambda(q) functions.
#         lambda_init : float, optional
#             SRG evolution parameter \lambda for initial Hamiltonian [fm^-1].
#             This allows one to use an SRG-evolved potential as the starting
#             point.
#         kvnn_hard : int, optional
#             Inputing a nonzero argument here will evolve the initial potential
#             (corresponding to kvnn) back to a harder scale using 
#             transformations from a harder potential (kvnn_hard).
            
#         """
        
#         # Get relevant info for file and directory names
#         # Part of data directory name
#         self.kvnn = kvnn 
#         # Part of file name
#         self.lamb = lamb
#         self.kmax = kmax
        
#         # Set-up for calculation
#         if interp == False:
            
#             # Channel is 3S1-3D1 for deuteron
#             channel = '3S1'
            
#             potential = Potential(kvnn, channel, kmax, kmid, ntot)

#             # Load and save momentum arrays
#             k_array, k_weights = potential.load_momentum(kvnn, channel, kmax,
#                                                          kmid, ntot)
#             # Save k_array, k_weights, ntot
#             self.k_array, self.k_weights, self.ntot = k_array, k_weights, ntot
            
#             # For dividing out momenta/weights
#             factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
#             # For coupled-channel matrices
#             factor_array_cc = np.concatenate( (factor_array, factor_array) )

#             # Load SRG transformation
                
#             # If lambda_init = np.inf, then take the initial Hamiltonian of
#             # kvnn as the starting point
#             if lambda_init == np.inf:
                
#                 H_initial = potential.load_hamiltonian()
                    
#             # Otherwise, this splits into two cases
#             else:
                    
#                 # If kvnn_hard is zero, take an SRG-evolved Hamiltonian
#                 # corresponding to kvnn and lambda_init as the initial
#                 # Hamiltonian
#                 if kvnn_hard == 0:
                    
#                     H_initial = potential.load_hamiltonian('srg', 'Wegner',
#                                                            lambda_init)
                    
#                 # If kvnn_hard is nonzero, SRG-evolve the initial Hamiltonian 
#                 # of kvnn back using transformations from kvnn_hard at
#                 # lambda_initial
#                 else:
                    
#                     potential_hard = Potential(kvnn_hard, channel, kmax, kmid,
#                                                ntot)
#                     H_hard_initial = potential_hard.load_hamiltonian()
#                     H_hard_evolved = potential_hard.load_hamiltonian(
#                         'srg', 'Wegner', lambda_init
#                     )
                            
#                     # Get SRG transformation from hard potential
#                     U_hard = get_transformation(H_hard_initial, H_hard_evolved)
                        
#                     # Get initial Hamiltonian for kvnn
#                     H_matrix = potential.load_hamiltonian()
                        
#                     # Do inverse transformation on softer Hamiltonian
#                     H_initial = U_hard.T @ H_matrix @ U_hard
                
#             H_evolved = potential.load_hamiltonian('srg', 'Wegner', lamb)
#             # Load U(k, k') [unitless]
#             U_matrix_unitless = get_transformation(H_initial, H_evolved)
            
#             # Unitless wave function in momentum space
#             psi_k_unitless = wave_function(H_initial, U=U_matrix_unitless)
#             self.psi = psi_k_unitless
#             # Divide out momenta/weights and save
#             self.psi_k = psi_k_unitless / factor_array_cc # [fm^3/2]

#             # Isolate 2-body term and convert to fm^3
#             I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
#             row, col = np.meshgrid(factor_array_cc, factor_array_cc)
#             delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
#             # Save for exact unitary calculation
#             self.delta_U_matrix = delta_U_matrix_unitless
#             delta_U_matrix = delta_U_matrix_unitless / row / col # fm^3
            
#             # No 2*J+1 factor since we're fixing M_J for deuteron
            
#             # Evaluate matrix elements
#             deltaU = 1/2 * ( delta_U_matrix[:ntot, :ntot] + \
#                              delta_U_matrix[ntot:, ntot:] )
#             deltaU2 = 1/2 * ( delta_U_matrix[:ntot, :ntot]**2 + \
#                               delta_U_matrix[:ntot, ntot:]**2 + \
#                               delta_U_matrix[ntot:, :ntot]**2 + \
#                               delta_U_matrix[ntot:, ntot:]**2 )
                
#             # Interpolate < k | \delta U | k >
#             self.deltaU_func = interp1d( k_array, np.diag(deltaU),
#                                          kind='linear', bounds_error=False,
#                                          fill_value='extrapolate' )

#             # Interpolate < k | \delta U \delta U^{\dagger} | k' >
#             self.deltaU2_func = RectBivariateSpline(k_array, k_array, deltaU2,
#                                                     kx=1, ky=1)


#     def theta_I(self, q_mesh, kF_mesh):
#         """
#         Evaluates \theta( kF(R) - q ). This function appears in the I term.

#         Parameters
#         ----------
#         q_mesh : 2-D ndarray
#             Momentum values [fm^-1].
#         kF_mesh : 2-D ndarray
#             Deuteron Fermi momentum [fm^-1] with respect to R.

#         Returns
#         -------
#         theta_mesh : 2-D ndarray
#             \theta function [unitless] evaluated for each q and kF(R).

#         """
        
#         # Initialize 2-D array
#         theta_mesh = np.zeros( (self.ntot_q, self.ntot_R) )
        
#         # Gives 1 if q < kF(R)
#         theta_mesh[ q_mesh < kF_mesh ] = 1
        
#         # This is a (ntot_q, ntot_R) size array
#         return theta_mesh
        

    # def theta_deltaU(self, q_mesh, kF_mesh, k_mesh):
    #     """
    #     Evaluates angle-average of \theta( kF(R) - \abs(q - 2k) ). This
    #     function appears in the \delta U term.

    #     Parameters
    #     ----------
    #     q_mesh : 3-D ndarray
    #         Momentum values [fm^-1].
    #     kF_mesh : 3-D ndarray
    #         Deuteron Fermi momentum [fm^-1] with respect to R.
    #     k_mesh : 3-D ndarray
    #         Relative momentum values [fm^-1].

    #     Returns
    #     -------
    #     theta_mesh : 3-D ndarray
    #         \theta function [unitless] evaluated for each q, kF(R), and k.
            
    #     Notes
    #     -----
    #     Not sure why the cases had to be computed in reverse order. Does that
    #     mean there is an overlap of truth values in case 3 and case 1 when
    #     there shouldn't be?

    #     """
        
    #     # Initialize 3-D array
    #     theta_mesh = np.zeros( (self.ntot_q, self.ntot_R, self.ntot_k) )
        
    #     # Evaluate each boolean case and use these to fill in the theta_mesh
        
    #     # This applies to each case: q < kF
    #     case_all = q_mesh < kF_mesh
        
    #     # Case 3: q-kF < 2k < kF+q
    #     case_3 = case_all * ( q_mesh - kF_mesh < 2*k_mesh ) * \
    #              ( 2*k_mesh < kF_mesh + q_mesh )
    #     theta_mesh[case_3] = ( ( kF_mesh**2 - ( q_mesh - 2*k_mesh )**2 ) / \
    #                            ( 8*k_mesh*q_mesh ) )[case_3]
            
    #     # Case 2: kF-q <= 2k < kF+q
    #     case_2 = case_all * ( kF_mesh - q_mesh <= 2*k_mesh ) * \
    #              ( 2*k_mesh < kF_mesh + q_mesh )
    #     theta_mesh[case_2] = ( ( kF_mesh**2 - ( q_mesh - 2*k_mesh )**2 ) / \
    #                            ( 8*k_mesh*q_mesh ) )[case_2]
        
    #     # Case 1: 2k < kF-q
    #     case_1 = case_all * ( 2*k_mesh < kF_mesh - q_mesh )
    #     theta_mesh[case_1] = 1

    #     # This is a (ntot_q, ntot_R, ntot_k) size array
    #     return theta_mesh


    # def theta_deltaU2(self, kF_mesh, K_mesh, k_mesh):
    #     """
    #     Evaluates angle-average of \theta( kF(R) - \abs(K/2 + k) ) x
    #     \theta( kF(R) - \abs(K/2 - k) ). This function appears in the
    #     \delta U \delta U^\dagger term.

    #     Parameters
    #     ----------
    #     kF_mesh : 4-D ndarray
    #         Deuteron Fermi momentum [fm^-1] with respect to R.
    #     K_mesh : 4-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     k_mesh : 4-D ndarray
    #         Relative momentum values [fm^-1].

    #     Returns
    #     -------
    #     theta_mesh : 4-D ndarray
    #         \theta function [unitless] evaluated for each q, kF(R), K, and k.
    #         Note, this function does not depend on q but we need the array to
    #         match the size of \delta U \delta U^\dagger matrix elements.

    #     """
        
    #     # Initialize 4-D array
    #     theta_mesh = np.zeros( (self.ntot_q, self.ntot_R, self.ntot_K,
    #                             self.ntot_k) )

    #     # Evaluate each boolean case and use these to fill in the theta_mesh
        
    #     # Case 2: 2k+K > 2kF and 4k^2+K^2 < 4kF^2
    #     case_2 = ( 2*k_mesh + K_mesh > 2*kF_mesh ) * \
    #              ( 4*k_mesh**2 + K_mesh**2 <= 4*kF_mesh**2 )    
    #     theta_mesh[case_2] = ( ( 4*kF_mesh**2 - 4*k_mesh**2 - K_mesh**2 ) / \
    #                            (4*k_mesh*K_mesh) )[case_2]

    #     # Case 1: 2k+K < 2kF
    #     case_1 = ( 2*k_mesh + K_mesh <= 2*kF_mesh )
    #     theta_mesh[case_1] = 1
        
    #     # This is a (ntot_q, ntot_R, ntot_K, ntot_k) size array
    #     return theta_mesh
    
    
    # def n_I(self, q_array, R_array, dR, kF_array):
    #     """
    #     Evaluates the I term in U n(q) U^\dagger ~ I n(q) I.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     kF_array : 1-D ndarray
    #         Deuteron Fermi momentum [fm^-1] with respect to R.

    #     Returns
    #     -------
    #     output : 1-D ndarray
    #         Momentum distribution [fm^3] from I term as a function of q.

    #     """
        
    #     # Initialize 2-D meshgrids (q, R)
    #     q_mesh, R_mesh = np.meshgrid(q_array, R_array, indexing='ij')
        
    #     # Get 2-D kF(R) mesh
    #     _, kF_mesh = np.meshgrid(q_array, kF_array, indexing='ij')
        
    #     # Evaluate the \theta-function in I term
    #     theta_mesh = self.theta_I(q_mesh, kF_mesh)
        
    #     # Calculate R integrand (ntot_q, ntot_R)
    #     integrand_R = theta_mesh * R_mesh**2 * dR

    #     # Integrate over R
    #     # This is a (ntot_q, 1) size array
    #     # Factor of 2 is overall factor
    #     return 2 * np.sum(integrand_R, axis=-1)
    
    
    # def n_deltaU(self, q_array, R_array, dR, kF_array):
    #     """
    #     Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
    #     we are combining \delta U and \delta U^\dagger, hence the factor of 2.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     kF_array : 1-D ndarray
    #         Deuteron Fermi momentum [fm^-1] with respect to R.

    #     Returns
    #     -------
    #     output : 1-D ndarray
    #         Momentum distribution [fm^3] from \delta U term as a function of q.

    #     """
        
    #     # Set number of k points for integration over k
    #     self.ntot_k = 40
        
    #     # Initialize 3-D meshgrids (q, R, k) with k values equal to 0
    #     q_mesh, R_mesh, k_mesh = np.meshgrid(q_array, R_array,
    #                                          np.zeros(self.ntot_k),
    #                                          indexing='ij')
        
    #     # Get 3-D kF(R) mesh and initialize k weights mesh
    #     _, kF_mesh, dk_mesh = np.meshgrid(q_array, kF_array,
    #                                        np.zeros(self.ntot_k),
    #                                        indexing='ij')

    #     # Loop over q and R to find limits of k integration and then create
    #     # k_array using Gaussian quadrature
    #     for iq, q in enumerate(q_array):
    #         for iR, R in enumerate(R_array):
                
    #             kF = kF_array[iR]

    #             # Create integration mesh k_array up to (kF + q)/2 which
    #             # corresponds to the upper limit of \theta( kF2(R) - |q-2k| )
    #             k_max = (kF + q)/2
 
    #             # Get Gaussian quadrature mesh
    #             k_array, k_weights = gaussian_quadrature_mesh(k_max,
    #                                                           self.ntot_k)
                
    #             # Fill in k_mesh and dk_mesh given the specific k_array
    #             k_mesh[iq, iR, :] = k_array
    #             dk_mesh[iq, iR, :] = k_weights
        
    #     # Evaluate angle-average of \theta-functions in \delta U term
    #     theta_mesh = self.theta_deltaU(q_mesh, kF_mesh, k_mesh)
        
    #     # Evaluate < k | \delta U | k >
    #     deltaU_mesh = self.deltaU_func(k_mesh)

    #     # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
    #     # combine to give 2
    #     # Factor of 2 from \delta U + \delta U^\dagger
    #     # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
    #     # 2/\pi for two | k_vec > -> | k J L S ... > changes
    #     deltaU_factor = 2 * 8 * 2/np.pi * 2
        
    #     # Calculate the k integrand where we split terms according to pp and
    #     # pn (or nn and np if kF_1 corresponds to a neutron)
    #     # (ntot_q, ntot_R, ntot_k)
    #     integrand_k = deltaU_factor * k_mesh**2 * dk_mesh * R_mesh**2 * dR * \
    #                   deltaU_mesh * theta_mesh
        
    #     # Integrate over k leaving R integrand (ntot_q, ntot_R)
    #     integrand_R = np.sum(integrand_k, axis=-1)

    #     # Integrate over R
    #     # This is a (ntot_q, 1) size array
    #     return np.sum(integrand_R, axis=-1)


    # def n_deltaU2(self, q_array, R_array, dR, kF_array):
    #     """
    #     Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     kF_array : 1-D ndarray
    #         Deuteron Fermi momentum [fm^-1] with respect to R.

    #     Returns
    #     -------
    #     output : 1-D ndarray
    #         Momentum distribution [fm^3] from \delta U \delta U^\dagger term as
    #         a function of q.

    #     """
        
    #     # Set number of K and k points for integration over K and k
    #     self.ntot_K = 40
    #     self.ntot_k = 40
        
    #     # y = cos(\theta) angles for averaging integration
    #     ntot_y = 7
    #     y_array, y_weights = leggauss(ntot_y)
        
    #     # Initialize 4-D meshgrids (q, R, K, k) with k and K equal to 0
    #     q_mesh, R_mesh, K_mesh, k_mesh = np.meshgrid(q_array, R_array,
    #                                                  np.zeros(self.ntot_K),
    #                                                  np.zeros(self.ntot_k),
    #                                                  indexing='ij')
        
    #     # Get 4-D kF(R) mesh and initialize K and k weights mesh
    #     _, kF_mesh, dK_mesh, dk_mesh = np.meshgrid(q_array, kF_array,
    #                                                 np.zeros(self.ntot_K),
    #                                                 np.zeros(self.ntot_k),
    #                                                 indexing='ij')

    #     # Loop over q, R, and k to find limits of K integration and then create
    #     # K_array using Gaussian quadrature
    #     for iR, R in enumerate(R_array):
                
    #         kF = kF_array[iR]
    
    #         # K integration goes from 0 to 2*kF
    #         K_max = 2*kF
                
    #         # Get Gaussian quadrature mesh for K integration
    #         K_array, K_weights = gaussian_quadrature_mesh(K_max, self.ntot_K)
              
    #         # Loop over remaining variables and fill in 4-D array
    #         for iq in range(self.ntot_q):
    #             for ik in range(self.ntot_k):
                    
    #                 # Fill in K_mesh and dK_mesh given the specific K_array
    #                 K_mesh[iq, iR, :, ik] = K_array
    #                 dK_mesh[iq, iR, :, ik] = K_weights

    #     # Loop over q, R, and K to find limits of k integration and then create
    #     # k_array using Gaussian quadrature
    #     for iR, R in enumerate(R_array):
        
    #         kF = kF_array[iR]
            
    #         # K_array only depends on R, so loop over K_mesh[0, iR, :, 0]
    #         for iK, K in enumerate( K_mesh[0, iR, :, 0] ):
                
    #             # Lower limit of k integration
    #             k_min = max(K/2 - kF, 0)
                
    #             # Upper limit of k integration
    #             if K**2/4 < kF**2:
    #                 k_max = min( np.sqrt(kF**2 - K**2/4), kF + K/2 )
    #             else:
    #                 k_max = kF + K/2

    #             # Get Gaussian quadrature mesh for k integration
    #             k_array, k_weights = gaussian_quadrature_mesh(k_max,
    #                                                           self.ntot_k,
    #                                                           xmin=k_min)
                
    #             # Loop over remaining variables and fill in 4-D array    
    #             for iq in range(self.ntot_q):
                    
    #                 # Fill in k_mesh and dk_mesh given the specific k_array
    #                 k_mesh[iq, iR, iK, :] = k_array
    #                 dk_mesh[iq, iR, iK, :] = k_weights
                        
    #     # Evaluate angle-average of \theta-functions in \delta U^2 term
    #     theta_mesh = self.theta_deltaU2(kF_mesh, K_mesh, k_mesh)

    #     # Evaluate 4-D < k | \delta U | |q_vec-K_vec/2| >^2 while also
    #     # averaging over angle y
    #     deltaU2_mesh = np.zeros_like(theta_mesh) # (q, R, K, k)
        
    #     # Integrate over y
    #     for y, dy in zip(y_array, y_weights):
            
    #         # Evaluate |q-K/2| mesh
    #         q_K_mesh = np.sqrt( q_mesh**2 + K_mesh**2/4 - q_mesh*K_mesh*y )
            
    #         deltaU2_mesh += self.deltaU2_func.ev(k_mesh, q_K_mesh) * dy/2
        
    #     # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
    #     # combine to give 2
    #     # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
    #     deltaU2_factor = 2 * (2/np.pi)**2

    #     # Calculate the k integrand
    #     # (ntot_q, ntot_R, ntot_K, ntot_k)
    #     integrand_k = deltaU2_factor * k_mesh**2 * dk_mesh * K_mesh**2 * \
    #                   dK_mesh * R_mesh**2 * dR * deltaU2_mesh * theta_mesh

    #     # Integrate over k leaving K integrand (ntot_q, ntot_R, ntot_K)
    #     integrand_K = np.sum(integrand_k, axis=-1)
        
    #     # Integrate over K leaving R integrand (ntot_q, ntot_R)
    #     integrand_R = np.sum(integrand_K, axis=-1)
        
    #     # Integrate over R
    #     # This is a (ntot_q, 1) size array
    #     return np.sum(integrand_R, axis=-1)
            
    
    # def n_total(self, q_array, R_array, dR):
    #     """
    #     Deuteron momentum distribution under HF+LDA.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].

    #     Returns
    #     -------
    #     n_total : 1-D ndarray
    #         Deuteron momentum distribution [fm^3] for each q.
            
    #     Notes
    #     -----
    #     Should R_array be called r_array for relative coordinate?

    #     """
        
    #     # Save lengths of q_array and R_array
    #     self.ntot_q = len(q_array)
    #     self.ntot_R = len(R_array)
        
    #     # Get deuteron density
    #     k_array, k_weights, ntot_k = self.k_array, self.k_weights, self.ntot
        
    #     # Transform wave function to coordinate space
    #     hank_trans_3S1 = hankel_transformation_k2r('3S1', k_array, k_weights,
    #                                                R_array)
    #     hank_trans_3D1 = hankel_transformation_k2r('3D1', k_array, k_weights,
    #                                                R_array)
    
    #     # Get 3S1 and 3D1 waves [fm^-3/2] in coordinate-space
    #     psi_R_3S1 = hank_trans_3S1 @ self.psi_k[:ntot_k]
    #     psi_R_3D1 = hank_trans_3D1 @ self.psi_k[ntot_k:]
    
    #     # Calculate the deuteron nucleonic density [fm^-3]
    #     rho_array = psi_R_3S1**2 + psi_R_3D1**2

    #     # Evaluate kF values at each point in R_array
    #     kF_array = (3*np.pi**2 * rho_array)**(1/3)
            
    #     # Get each contribution with respect to q
    #     n_I = self.n_I(q_array, R_array, dR, kF_array)
    #     n_deltaU = self.n_deltaU(q_array, R_array, dR, kF_array)
    #     n_deltaU2 = self.n_deltaU2(q_array, R_array, dR, kF_array)

    #     # Return total (ntot_q, 1)
    #     return n_I + n_deltaU + n_deltaU2
    
    
    # def n_contributions(self, q_array, R_array, dR):
    #     """
    #     Contributions to the deuteron momentum distribution. This function
    #     isolates the I, \delta U, and \delta U \delta U^\dagger terms.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].

    #     Returns
    #     -------
    #     n_contributions : tuple
    #         Tuple of 1-D ndarrays corresponding to the I, \delta U, and
    #         \delta U^\dagger terms of the pair momentum distribution [fm^3] for
    #         each q.

    #     """
        
    #     # Save lengths of q_array and R_array
    #     self.ntot_q = len(q_array)
    #     self.ntot_R = len(R_array)
        
    #     # Get deuteron density
    #     k_array, k_weights, ntot_k = self.k_array, self.k_weights, self.ntot
        
    #     # Transform wave function to coordinate space
    #     hank_trans_3S1 = hankel_transformation_k2r('3S1', k_array, k_weights,
    #                                                R_array)
    #     hank_trans_3D1 = hankel_transformation_k2r('3D1', k_array, k_weights,
    #                                                R_array)
    
    #     # Get 3S1 and 3D1 waves [fm^-3/2] in coordinate-space
    #     psi_R_3S1 = hank_trans_3S1 @ self.psi_k[:ntot_k]
    #     psi_R_3D1 = hank_trans_3D1 @ self.psi_k[ntot_k:]
    
    #     # Calculate the deuteron nucleonic density [fm^-3]
    #     rho_array = psi_R_3S1**2 + psi_R_3D1**2

    #     # Evaluate kF values at each point in R_array
    #     kF_array = (3*np.pi**2 * rho_array)**(1/3)
            
    #     # Get each contribution with respect to q
    #     n_I = self.n_I(q_array, R_array, dR, kF_array)
    #     n_deltaU = self.n_deltaU(q_array, R_array, dR, kF_array)
    #     n_deltaU2 = self.n_deltaU2(q_array, R_array, dR, kF_array)
        
    #     # Return tuple of contributions ( (ntot_q, ntot_Q), ... )
    #     return n_I, n_deltaU, n_deltaU2
    
    
    # def write_file(self):
    #     """
    #     Write deuteron momentum distribution file for interpolation purposes.
    #     Split things into total, I, \delta U, and \delta U^2 contributions.

    #     """
        
    #     # Get momentum values
    #     q_array = self.k_array
        
    #     # Directory for distributions data
    #     data_directory = 'data/dmd/kvnn_%d' % self.kvnn
        
    #     # Create file name
    #     file_name = 'deuteron_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
    #     file_name = replace_periods(file_name) + '.dat'
        
        
    #     # Use R values in accordance with densities code from densities.py
    #     # Nuclei arguments don't matter here - just getting R values
    #     R_array, _ = load_density('C12', 'proton', 6, 6)
    #     dR = R_array[2] - R_array[1] # Assuming linear spacing

    #     # Calculate n_\lambda^d(q) for each q in q_array
    #     n_I_array, n_delU_array, n_delU2_array = self.n_contributions(q_array,
    #                                                                   R_array,
    #                                                                   dR)
            
    #     # Total momentum distribution
    #     n_total_array = n_I_array + n_delU_array + n_delU2_array
    
    #     # Open file and write header where we allocate roughly 18 centered
    #     # spaces for each label
    #     f = open(data_directory + '/' + file_name, 'w')
    #     header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}'.format('q',
    #               'total', '1', '\delta U', '\delta U^2')
    #     f.write(header + '\n')
    
    #     # Loop over momenta q
    #     for iq, q in enumerate(q_array):

    #         # Write to data file following the format from the header
    #         line = '{:^18.6f}{:^18.6e}{:^18.6e}{:^18.6e}{:^18.6e}'.format( q,
    #                         n_total_array[iq], n_I_array[iq], n_delU_array[iq],
    #                         n_delU2_array[iq] )
    #         f.write('\n' + line)

    #     # Close file
    #     f.close()
        
        
    # def n_lambda_interp(self):
    #     """
    #     Interpolate the deuteron momentum distribution for the specified file.
    
    #     Returns
    #     -------
    #     output : tuple
    #         Tuple of functions that depend only on momentum q [fm^-1] where
    #         each function corresponds to contributions to n_\lambda(q): total,
    #         I, \delta U, and \delta U^2.

    #     """
        
    #     # Directory for distributions data
    #     data_directory = 'data/dmd/kvnn_%d' % self.kvnn
        
    #     # Get file name
    #     file_name = 'deuteron_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
    #     file_name = replace_periods(file_name) + '.dat'
        
    #     # Load data which includes all contributions to n_\lambda(q)
    #     data = np.loadtxt(data_directory + '/' + file_name)
        
    #     # Split data into 1-D arrays for each column
    #     q_array = data[:, 0] # Momentum in fm^-1
    #     n_total_array = data[:, 1] # Total distribution
    #     n_I_array = data[:, 2] # 1 term
    #     n_delU_array = data[:, 3] # \delta U term
    #     n_delU2_array = data[:, 4] # \delta U^2 term
        
    #     # Interpolate each array (UnivariateSpline is for smoothing whereas
    #     # interp1d gives closer value to the actual calculation)
    #     n_total_func = interp1d(q_array, n_total_array, bounds_error=False,
    #                             kind='linear', fill_value='extrapolate')
    #     n_I_func = interp1d(q_array, n_I_array, bounds_error=False,
    #                         kind='linear', fill_value='extrapolate')
    #     n_delU_func = interp1d(q_array, n_delU_array, bounds_error=False,
    #                            kind='linear', fill_value='extrapolate')
    #     n_delU2_func = interp1d(q_array, n_delU2_array, bounds_error=False,
    #                             kind='linear', fill_value='extrapolate')
        
    #     # Return all contributions with total first
    #     # Note, these are functions of q
    #     return n_total_func, n_I_func, n_delU_func, n_delU2_func
    
    
    # def n_lambda_exact(self, q, contributions='total'):
    #     """
    #     Computes the exact deuteron momentum distribution using the wave
    #     function from the input Hamiltonian. Note, set interp = False to use
    #     this function.

    #     Parameters
    #     ----------
    #     q : float
    #         Momentum value [fm^-1].
    #     contributions : str, optional
    #         Option to return different contributions to the momentum
    #         distribution.
    #         1. Default is 'total' which only returns the total momentum
    #            distribution.
    #         2. Specify 'q_contributions' for I, \delta U, and 
    #            \delta U \delta U^\dagger along with the total.
    #         3. Specify 'partial_wave_ratio' for the ratio of 3S1-3S1 to full
    #            3S1-3D1 along with the absolute total.

    #     Returns
    #     -------
    #     output : float or tuple
    #         Pair momentum distribution [fm^3] evaluated at momentum q
    #         [fm^-1]. (Note, this function will return a tuple of floats if
    #         contributions is not 'total'.)

    #     """

    #     # Load deuteron wave function and delta_U_matrix (both unitless)
    #     psi_vector = self.psi
    #     delta_U_matrix = self.delta_U_matrix
        
    #     # Load bare momentum projection operator [fm^3]
    #     bare_operator = momentum_projection_operator(q, self.k_array,
    #                     self.k_weights, '3S1', smeared=False)
        
    #     # Term 1: 1 * n(q) * 1
    #     term_1 = psi_vector.T @ bare_operator @ psi_vector
        
    #     # Term 2: \deltaU * n(q) * 1
    #     term_deltaU = psi_vector.T @ delta_U_matrix @ bare_operator @ \
    #                   psi_vector
        
    #     # Term 3: 1 * n(q) * \deltaU^\dagger = Term 2
    #     term_deltaU *= 2
        
    #     # High-q term: \deltaU * n(q) * \deltaU^\dagger
    #     term_deltaU2 = psi_vector.T @ delta_U_matrix @ bare_operator @ \
    #                     delta_U_matrix.T @ psi_vector
                       
    #     if contributions == 'partial_wave_ratio':
            
    #         # Length of q_array
    #         ntot = self.ntot
            
    #         # 3S1-3S1 only
    #         numerator = abs( psi_vector.T[:ntot] @ ( \
    #                              delta_U_matrix[:ntot, :ntot] @ \
    #                              bare_operator[:ntot, :ntot] @ \
    #                              delta_U_matrix.T[:ntot, :ntot] +
    #                              delta_U_matrix[:ntot, ntot:] @ \
    #                              bare_operator[ntot:, ntot:] @ \
    #                              delta_U_matrix.T[ntot:, :ntot] ) @ \
    #                          psi_vector[:ntot] )
                        
    #         # Full 3S1-3D1 taking absolute values 
    #         denominator = numerator + \
    #                       abs( psi_vector.T[:ntot] @ ( \
    #                                delta_U_matrix[:ntot, :ntot] @ \
    #                                bare_operator[:ntot, :ntot] @ \
    #                                delta_U_matrix.T[:ntot, ntot:] +
    #                                delta_U_matrix[:ntot, ntot:] @ \
    #                                bare_operator[ntot:, ntot:] @ \
    #                                delta_U_matrix.T[ntot:, ntot:] ) @ \
    #                            psi_vector[ntot:] ) + \
    #                       abs( psi_vector.T[ntot:] @ ( \
    #                                delta_U_matrix[ntot:, :ntot] @ \
    #                                bare_operator[:ntot, :ntot] @ \
    #                                delta_U_matrix.T[:ntot, :ntot] +
    #                                delta_U_matrix[ntot:, ntot:] @ \
    #                                bare_operator[ntot:, ntot:] @ \
    #                                delta_U_matrix.T[ntot:, :ntot] ) @ \
    #                            psi_vector[:ntot] ) + \
    #                       abs( psi_vector.T[ntot:] @ ( \
    #                                delta_U_matrix[ntot:, :ntot] @ \
    #                                bare_operator[:ntot, :ntot] @ \
    #                                delta_U_matrix.T[:ntot, ntot:] +
    #                                delta_U_matrix[ntot:, ntot:] @ \
    #                                bare_operator[ntot:, ntot:] @ \
    #                                delta_U_matrix.T[ntot:, ntot:] ) @ \
    #                            psi_vector[ntot:] )
                                   
    #     # Add up each term for total
    #     total = term_1 + term_deltaU + term_deltaU2

    #     # Return contributions and total or just total
    #     if contributions == 'q_contributions':
    #         return total, term_1, term_deltaU, term_deltaU2
    #     elif contributions == 'partial_wave_ratio':
    #         return total, numerator/denominator
    #     else: # Default
    #         return total