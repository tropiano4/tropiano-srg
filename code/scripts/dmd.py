#!/usr/bin/env python3

"""
File: dmd.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 17, 2022

The Deuteron class calculates deuteron momentum distributions (dmd) SRG-
evolving the operator assuming the evolved wave function is given by HF
treated in LDA. This class is a subclass of the MomentumDistribution class
from momentum_distributions.py.

Notes on normalizations:
  1. The deuteron wave function describing relative position or momentum is
      normalized according to
        \int dr r^2 (|\psi_{3S1}(r)|^2 + |\psi_{3D1}|^2) = 1,
        2/\pi * \int dk k^2 (|\psi_{3S1}(k)|^2 + |\psi_{3D1}(k)|^2) = 1.
  2. Under HF+LDA, we adopt the normalization
        4\pi / (2\pi)^3 \int dk k^2 < n_d(k) > = 1,
     where angled-brackets indicate nuclear-averaging (integration over R).

Last update: July 28, 2022

"""

# Todo: There are several things that overlap (if not slightly) with the
#  single-nucleon momentum distribution code. This suggests making Deuteron
#  inherit SingleNucleon.
# Todo: Think about integration over K. Does this even make sense?

# Python imports
import numpy as np
from numpy.polynomial.legendre import leggauss

# Imports from A.T. codes
from .densities import load_density
from .fourier_transform import hankel_transformation_k2r
from .integration import gaussian_quadrature_mesh
from .momentum_distributions import MomentumDistribution
from .tools import replace_periods


class Deuteron(MomentumDistribution):
    """
    Computes deuteron momentum distributions SRG-evolving the operator
    assuming the evolved wave function is given by HF treated in LDA.

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
        SRG evolution parameter \lambda for initial Hamiltonian [fm^-1]. This
        allows one to use an SRG-evolved potential as the starting point.
    kvnn_inv : int, optional
        This number specifies a potential for which inverse-SRG transformations
        will be applied to the initial Hamiltonian
            H_initial = U_{kvnn_inv}^{\dagger} H_kvnn U_{kvnn_inv},
        where the transformations are evaluated at \delta \lambda.
    lambda_m : float, optional
        SRG evolution parameter \lambda for inverse-SRG transformations
        [fm^-1]. Note, both kvnn_inv and lambda_m must be specified for
        this to run.

    """

    def __init__(
            self, kvnn, kmax, kmid, ntot, generator, lamb, lambda_initial=None,
            kvnn_inv=None, lambda_m=None):
        """Sets \delta U(k,k') and \delta U^2(k,k') functions as attributes."""

        super().__init__(kvnn, kmax, kmid, ntot)

        # Set-up calculation
        self.save_deuteron_deltaU_funcs(generator, lamb, lambda_initial,
                                        kvnn_inv, lambda_m)

        # Set instance attributes for saving files
        self.generator = generator
        self.lamb = lamb
        self.lambda_initial = lambda_initial
        self.kvnn_inv = kvnn_inv
        self.lambda_m = lambda_m

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
        mask_3 = (common_constraint * (q_grid - kF_grid < 2 * k_grid)
                  * (2 * k_grid < kF_grid + q_grid))
        angle_avg_grid[mask_3] = ((kF_grid ** 2 - (q_grid - 2 * k_grid) ** 2)
                                  / (8 * k_grid * q_grid))[mask_3]

        # Case 2: q < kF and kF-q <= 2k < kF+q
        mask_2 = (common_constraint * (q_grid < kF_grid) *
                  (kF_grid - q_grid <= 2 * k_grid) * (
                          2 * k_grid < kF_grid + q_grid))
        angle_avg_grid[mask_2] = ((kF_grid ** 2 - (q_grid - 2 * k_grid) ** 2)
                                  / (8 * k_grid * q_grid))[mask_2]

        # Case 1: q < kF and 2k < kF-q
        mask_1 = (common_constraint * (q_grid < kF_grid)
                  * (2 * k_grid < kF_grid - q_grid))
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
        mask_2 = ((2 * k_grid + K_grid > 2 * kF_grid)
                  * (4 * k_grid ** 2 + K_grid ** 2 <= 4 * kF_grid ** 2))
        angle_avg_grid[mask_2] = \
            ((4 * kF_grid ** 2 - 4 * k_grid ** 2 - K_grid ** 2)
             / (4 * k_grid * K_grid))[mask_2]

        # Case 1: 2k+K < 2kF
        mask_1 = (2 * k_grid + K_grid <= 2 * kF_grid)
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
            C.o.M. position [fm].
        dR : float
            C.o.M. position step-size (assuming linearly-spaced array) [fm].
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
        integrand_R = theta_grid * R_grid ** 2 * dR

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
            C.o.M. position [fm].
        dR : float
            C.o.M. position step-size (assuming linearly-spaced array) [fm].
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
                k_max = (kF + q) / 2

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
        deltaU_factor = 2 * 8 * 2 / np.pi * 2

        # Calculate the k integrand leaving (ntot_q, ntot_R, ntot_k)
        integrand_k = (deltaU_factor * k_grid ** 2 * dk_grid * R_grid ** 2 * dR
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
            C.o.M. position [fm].
        dR : float
            C.o.M. position step-size (assuming linearly-spaced array) [fm].
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
            K_max = 2 * kF

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
                k_min = max(K / 2 - kF, 0)

                # Upper limit of k integration
                if K ** 2 / 4 < kF ** 2:
                    k_max = min(np.sqrt(kF ** 2 - K ** 2 / 4), kF + K / 2)
                else:
                    k_max = kF + K / 2

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
            q_K_grid = np.sqrt(
                q_grid ** 2 + K_grid ** 2 / 4 - q_grid * K_grid * y)

            deltaU2_grid += self.deltaU2_func.ev(k_grid, q_K_grid) * dy / 2

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2 / np.pi) ** 2

        # Calculate the k integrand (ntot_q, ntot_R, ntot_K, ntot_k)
        integrand_k = (deltaU2_factor * k_grid ** 2 * dk_grid * K_grid ** 2
                       * dK_grid * R_grid ** 2 * dR * deltaU2_grid
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
        R_array, _ = load_density('proton', 'C12', 6, 6, 'SLy4')
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
        rho_array = psi_R_3S1 ** 2 + psi_R_3D1 ** 2

        # Evaluate kF values at each point in R_array
        kF_array = (3 * np.pi ** 2 * rho_array) ** (1 / 3)

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

        if self.lambda_initial is not None:
            file_name += f'_lambda_initial_{self.lambda_initial}'

        if self.kvnn_inv is not None:
            file_name += f'_kvnn_inv_{self.kvnn_inv}_lamb_m_{self.lambda_m}'

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
