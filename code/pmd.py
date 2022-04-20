#!/usr/bin/env python3

"""
File: pmd.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 17, 2022

The Pair class calculates pair momentum distributions (pmd) SRG-evolving the 
operator assuming the evolved wave function is given by HF treated in LDA.
This class is a sub-class of the MomentumDistribution class from
momentum_distributions.py.

Last update: April 20, 2022

"""

# To-do: Make sure R_array parameter is described correctly.

# Python imports
import numpy as np

# Imports from A.T. codes
from densities import load_density
from modules.integration import gaussian_quadrature_mesh
from modules.labels import replace_periods
from momentum_distributions import MomentumDistribution


class Pair(MomentumDistribution):

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
        
    def angle_avg(self, q_grid, Q_grid, kF1_grid, kF2_grid):
        """
        Evaluate angle-average involving Heaviside step functions:
            \int dz/2 \theta(kF1 - |Q/2+q|) \theta(kF2 - |Q/2-q|),
        where z is the angle between Q_vector and q_vector. Note, this 
        function will work for all three terms in the momentum distribution.
        
        Parameters
        ----------
        q_grid : 3-D or 4-D ndarray
            Meshgrid of relative momentum values [fm^-1].
        Q_grid : 3-D or 4-D ndarray
            Meshgrid of C.o.M. momentum values [fm^-1].
        kF1_grid : 3-D or 4-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau with respect to R.
        kF2_grid : 3-D or 4-D ndarray
            Meshgrid of Fermi momenta [fm^-1] for the nucleon corresponding to
            \tau' with respect to R or R'.

        Returns
        -------
        angle_avg_grid : 3-D 4-D ndarray
            Angle-averaging integral [unitless] that weights the corresponding
            term in the momentum distribution evaluated on meshgrids of kF1(R),
            kF2(R) (or kF2(R')), Q, and q (or k).

        """
        
        angle_avg_grid = np.zeros_like(q_grid)
        
        # Evaluate each boolean case and use these to fill in the meshgrid
        
        # Case 1: 2q+Q < 2kF1 and 2q+Q < 2kF2
        mask_1 = ((2*q_grid+Q_grid <= 2*kF1_grid)
                  * (2*q_grid+Q_grid <= 2*kF2_grid))
        angle_avg_grid[mask_1] = 1
        
        # If Q = 0, skip the following cases
        if self.ntot_Q == 1:
            return angle_avg_grid
        
        # Case 2: 2q+Q > 2kF1 and 2q+Q > 2kF2 and 4q^2+Q^2 < 2(kF1^2+kF2^2)
        mask_2 = ((2*q_grid+Q_grid > 2*kF1_grid)
                  * (2*q_grid+Q_grid > 2*kF2_grid)
                  * (4*q_grid**2+Q_grid**2 <= 2*(kF1_grid**2+kF2_grid**2))) 
        angle_avg_grid[mask_2] = ((2*(kF1_grid**2+kF2_grid**2) 
                                   - 4*q_grid**2-Q_grid**2)
                                  / (4*q_grid*Q_grid))[mask_2]
        
        # Case 3: 2q+Q < 2kF2 and -4 < (4q^2 - 4kF1^2 + Q^2)/(qQ) < 4
        mask_3 = ((2*q_grid+Q_grid <= 2*kF2_grid)
                  * (-4 < (4*q_grid**2-4*kF1_grid**2+Q_grid**2)
                     / (q_grid*Q_grid))
                  * ((4*q_grid**2-4*kF1_grid**2+Q_grid**2)
                     / (q_grid*Q_grid) <= 4 ))
        angle_avg_grid[mask_3] = ((4*kF1_grid**2-(Q_grid-2*q_grid)**2)
                                  / (8*q_grid*Q_grid))[mask_3]
        
        # Case 4: 2q+Q < 2kF1 and -4 < (4q^2-4kF2^2+Q^2)/(qQ) < 4
        mask_4 = ((2*q_grid+Q_grid <= 2*kF1_grid)
                  * (-4 < (4*q_grid**2-4*kF2_grid**2+Q_grid**2)
                     / (q_grid*Q_grid))
                  * ((4*q_grid**2-4*kF2_grid**2+Q_grid**2)
                     /(q_grid*Q_grid) <= 4))
        angle_avg_grid[mask_4] = ((4*kF2_grid**2-(Q_grid-2*q_grid)**2 )
                                  / (8*q_grid*Q_grid))[mask_4]

        # Shape of this meshgrid depends on the term we're evaluating for
        # I term -> (ntot_q, ntot_Q, ntot_R, ntot_Rp)
        # \delta U term -> (ntot_q, ntot_Q, ntot_R)
        # \delta U^2 term -> (ntot_q, ntot_Q, ntot_R, ntot_k)
        return angle_avg_grid
        
    def get_I_term(
            self, q_array, Q_array, R_array, dR, kF1_array, kF2_array=None):
        """
        Evaluates the I term in U n(q,Q) U^\dagger ~ I n(q,Q) I.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_array : 1-D ndarray, optional
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R'. If nothing is input, the function will
            evaluate for pp or nn assuming Fermi momentum from kF1_array.

        Returns
        -------
        output : 2-D ndarray
            Momentum distribution [fm^6] from I term as a function of q and Q.

        """
        
        if kF2_array is None:
            kF2_array = kF1_array  # pp or nn pair
        
        # Initialize 4-D meshgrids (q, Q, R, R')
        q_grid, Q_grid, R_grid, Rp_grid = np.meshgrid(
            q_array, Q_array, R_array, R_array, indexing='ij')
        # Get 4-D kF1(R) and kF2(R') meshgrids
        _, _, kF1_grid, kF2_grid = np.meshgrid(
            q_array, Q_array, kF1_array, kF2_array, indexing='ij')
        
        # Evaluate angle-average over Heaviside step functions in I term
        angle_avg_grid = self.angle_avg(q_grid, Q_grid, kF1_grid, kF2_grid)
        
        # Calculate R' integrand (ntot_q, ntot_Q, ntot_R, ntot_R)
        integrand_Rp = angle_avg_grid * Rp_grid**2 * dR * R_grid**2 * dR
        
        # Integrate over R' leaving R integrand (ntot_q, ntot_Q, ntot_R)
        integrand_R = 4*np.pi * np.sum(integrand_Rp, axis=-1)
        
        # Integrate over R leaving a (ntot_q, ntot_Q) shape array
        # Factor of 2 is overall factor for summation over spin projections
        return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def get_deltaU_term(
            self, q_array, Q_array, R_array, dR, kF1_array, kF2_array=None):
        """
        Evaluates second and third terms in U n(q,Q) U^\dagger ~ \delta U.
        Here we are combining \delta U and \delta U^\dagger, hence the factor
        of 2.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_array : 1-D ndarray, optional
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R. If nothing is input, the function will
            evaluate for pp or nn assuming Fermi momentum from kF1_array.

        Returns
        -------
        output : 2-D ndarray
            Momentum distribution [fm^6] from \delta U term as a function of q
            and Q.

        """
        
        # Set \delta U function based on pn vs pp (or nn)
        if kF2_array is None:
            kF2_array = kF1_array  # pp or nn pair
            deltaU_func = self.deltaU_pp_func
        else:
            deltaU_func = self.deltaU_pn_func
        
        # Initialize 3-D meshgrids (q, Q, R)
        q_grid, Q_grid, R_grid = np.meshgrid(q_array, Q_array, R_array,
                                             indexing='ij')
        
        # Get 3-D kF1(R) and kF2(R) meshes
        _, _, kF1_grid = np.meshgrid(q_array, Q_array, kF1_array,
                                     indexing='ij')
        _, _, kF2_grid = np.meshgrid(q_array, Q_array, kF2_array,
                                     indexing='ij')
        
        # Evaluate angle-average with Heaviside step functions in \delta U term
        angle_avg_grid = self.angle_avg(q_grid, Q_grid, kF1_grid, kF2_grid)
        
        # Evaluate \delta U(q,q) for \tau and \tau'
        deltaU_grid = deltaU_func(q_grid)
            
        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 1
        # Factor of 2 from \delta U + \delta U^\dagger
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        # 1/(4\pi) for averaging over \int d\Omega_q
        deltaU_factor =  2 * 2/np.pi * (2*np.pi)**3/(4*np.pi)
        
        # Calculate R integrand (ntot_q, ntot_Q, ntot_R)
        integrand_R = (deltaU_factor * deltaU_grid * angle_avg_grid * R_grid**2
                       * dR)
        
        # Integrate over R leaving a (ntot_q, ntot_Q) shape array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def get_deltaU2_term(
            self, q_array, Q_array, R_array, dR, kF1_array, kF2_array=None):
        """
        Evaluates fourth term in U n(q,Q) U^\dagger ~ \delta U^2.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
            with respect to R.
        kF2_array : 1-D ndarray, optional
            Fermi momentum [fm^-1] for the second nucleon corresponding to
            \tau' with respect to R. If nothing is input, the function will
            evaluate for pp or nn assuming Fermi momentum from kF1_array.

        Returns
        -------
        output : 2-D ndarray
            Momentum distribution [fm^6] from \delta U \delta U^\dagger term
            as a function of q and Q.

        """
        
        # Set \delta U^2 function based on pn vs pp (or nn)
        if kF2_array is None:
            kF2_array = kF1_array  # pp or nn pair
            deltaU2_func = self.deltaU2_pp_func
        else:
            deltaU2_func = self.deltaU2_pn_func
        
        ntot_q = len(q_array)
        # Set number of k points for integration over k
        ntot_k = 40
        
        # Initialize 4-D meshgrids (q, Q, R, k) with k values equal to 0
        q_grid, Q_grid, R_grid, k_grid = np.meshgrid(
            q_array, Q_array, R_array, np.zeros(ntot_k), indexing='ij')
        
        # Get 4-D kF1(R) and kF2(R) meshgrids and initialize k weights
        _, _, kF1_grid, dk_grid = np.meshgrid(
            q_array, Q_array, kF1_array, np.zeros(ntot_k), indexing='ij')
        _, _, kF2_grid, _ = np.meshgrid(q_array, Q_array, kF2_array,
                                        np.zeros(ntot_k), indexing='ij')

        # Loop over q, Q, and R to find limits of k integration and then create
        # k integration mesh using Gaussian quadrature
        for iQ, Q in enumerate(Q_array):
            for iR, R in enumerate(R_array):
                
                kF1, kF2 = kF1_array[iR], kF2_array[iR]
                
                # Minimum kF value
                kF_min = min(kF1, kF2)
                
                # Lower limit of integration
                k_min = max(Q/2 - kF_min, 0)
                
                # Upper limit of integration
                if Q**2/4 < (kF1**2 + kF2**2)/2:
                    k_max = min(np.sqrt((kF1**2+kF2**2)/2-Q**2/4), kF_min+Q/2)
                else:
                    k_max = kF_min + Q/2
                    
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max, ntot_k,
                                                              xmin=k_min)
                
                # Fill in k_grid and dk_grid given the specific k_array
                for iq in range(ntot_q):
                    k_grid[iq, iQ, iR, :] = k_array
                    dk_grid[iq, iQ, iR, :] = k_weights
                     
        # Evaluate angle-average with Heaviside step functions in \delta U^2
        # term (inputing k here NOT q)
        angle_avg_grid = self.angle_avg(k_grid, Q_grid, kF1_grid, kF2_grid)
        
        # Evaluate \delta U(k,q) \delta U^\dagger(q,k) for \tau and \tau'
        deltaU2_grid = deltaU2_func.ev(k_grid, q_grid)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 1
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        # 1/(4\pi) for averaging over \int d\Omega_q
        deltaU2_factor = (2/np.pi)**2 * (2*np.pi)**3/(4*np.pi)
        
        # Calculate k integrand (ntot_q, ntot_Q, ntot_R, ntot_k)
        integrand_k = (deltaU2_factor * deltaU2_grid * angle_avg_grid
                       * k_grid**2 * dk_grid * R_grid**2 * dR)
        
        # Integrate over k leaving R integrand (ntot_q, ntot_Q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R leaving a (ntot_q, ntot_Q) shape array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def compute_momentum_distribution(
            self, q_array, Q_array, pair, nucleus_name, Z, N, density='SLY4',
            save=False):
        """
        Pair momentum distribution for a specified nucleus.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
        pair : str
            Type of pair momentum distribution ('pp', 'pn', or 'nn').
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
        n_total : 2-D ndarray
            Pair momentum distribution n^{\tau \tau'}(q,Q) [fm^6].

        """
        
        self.ntot_Q = len(Q_array)
        
        # Split into cases based on pair
        if pair == 'pp':
            
            # Load R values and nucleonic densities
            R_array, rho_p_array = load_density('proton', nucleus_name, Z, N,
                                                density)
            
            # Evaluate kF values at each point in R_array
            kF1_array = (3*np.pi**2 * rho_p_array)**(1/3)
            kF2_array = None
            
        elif pair == 'nn':

            # Load R values and nucleonic densities
            R_array, rho_n_array = load_density('neutron', nucleus_name, Z, N,
                                                density)
            
            # Evaluate kF values at each point in R_array
            kF1_array = (3*np.pi**2 * rho_n_array)**(1/3)
            kF2_array = None
            
        else:
            
            R_array, rho_p_array = load_density('proton', nucleus_name, Z, N,
                                                density)
            R_array, rho_n_array = load_density('neutron', nucleus_name, Z, N,
                                                density)
            
            # Evaluate kF values at each point in R_array
            kF1_array = (3*np.pi**2 * rho_p_array)**(1/3)
            kF2_array = (3*np.pi**2 * rho_n_array)**(1/3)
        
        # Assuming linear spacing
        dR = R_array[2] - R_array[1]
        
        # Get each contribution with respect to q and Q (ntot_q, ntot_Q)
        n_I = self.get_I_term(q_array, Q_array, R_array, dR, kF1_array,
                              kF2_array)
        n_deltaU = self.get_deltaU_term(q_array, Q_array, R_array, dR,
                                        kF1_array, kF2_array)
        n_deltaU2 = self.get_deltaU2_term(q_array, Q_array, R_array, dR,
                                          kF1_array, kF2_array)

        n_total = n_I + n_deltaU + n_deltaU2
        
        if save:
            self.save_momentum_distribution(
                q_array, Q_array, n_total, n_I, n_deltaU, n_deltaU2, pair,
                nucleus_name, density)

        return n_total
    
    def save_momentum_distribution(
            self, q_array, Q_array, n_total, n_I, n_deltaU, n_deltaU2, pair,
            nucleus_name, density):
        """
        Saves momentum distribution in data/momentum_distributions. Here the
        specifications of the potential, SRG evolution, pair, and density
        dictate the file name.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        Q_array : 1-D ndarray
            C.o.M. momentum values [fm^-1].
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
        pair : str
            Type of pair momentum distribution ('pp', 'pn', or 'nn').
        nucleus_name : str
            Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
        density : str
            Name of nucleonic density (e.g., 'SLY4', 'Gogny').

        """
        
        # Directory for distributions data
        data_directory = f'../data/momentum_distributions/{nucleus_name}/'
        
        # Create file name
        file_name = (f'n_{pair}_{density}_kvnn_{self.kvnn}_kmax_{self.kmax}'
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
        
        # Split into cases on whether Q = 0
        if self.ntot_Q == 1:  # Q = 0
        
            file_name = replace_periods(file_name + '_Q0') + '.dat'
        
            # Open file and write header where we allocate 18 centered spaces
            # for each label
            f = open(data_directory + file_name, 'w')
        
            header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}'.format(
                'q', 'total', '1', '\delta U', '\delta U^2')
            f.write(header + '\n')
    
            # Loop over momenta q and write each contribution
            for iq, q in enumerate(q_array):

                # Write to data file following the format from the header
                line = (f'{q:^18.6f}{n_total[iq,0]:^18.6e}{n_I[iq,0]:^18.6e}'
                        f'{n_deltaU[iq,0]:^18.6e}{n_deltaU2[iq,0]:^18.6e}')
                f.write('\n' + line)

            f.close()
        
        else:  # Q != 0
        
            file_name = replace_periods(file_name) + '.dat'
        
            # Open file and write header where we allocate 18 centered spaces
            # for each label
            f = open(data_directory + file_name, 'w')
            
            header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}{:^18s}'.format(
                'q', 'Q', 'total', '1', '\delta U', '\delta U^2')
            f.write(header + '\n')
    
            # Loop over momenta q and write each contribution
            for iq, q in enumerate(q_array):
                for iQ, Q in enumerate(Q_array):

                    # Write to data file following the format from the header
                    line = (f'{q:^18.6f}{Q:^18.6f}{n_total[iq,iQ]:^18.6e}'
                            f'{n_I[iq,iQ]:^18.6e}{n_deltaU[iq,iQ]:^18.6e}'
                            f'{n_deltaU2[iq,iQ]:^18.6e}')
                    f.write('\n' + line)

            f.close()
        
    def compute_momentum_distribution_Q0(
            self, q_array, pair, nucleus_name, Z, N, density='SLY4',
            save=False):
        """
        Pair momentum distribution for a specified nucleus where the C.o.M.
        momentum iz zero Q=0. Note, this will save a separate file than the
        more general Q-dependent momentum distribution.

        Parameters
        ----------
        q_array : 1-D ndarray
            Relative momentum values [fm^-1].
        pair : str
            Type of pair momentum distribution ('pp', 'pn', or 'nn').
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
            Pair momentum distribution n^{\tau \tau'}(q,Q=0) [fm^6].

        """
        
        # Create Q_array with a single zero entry
        # Heaviside step functions simplify to give 1 or 0 in this case
        Q_array = np.array([0.0])
        
        # Call more general function
        n_total = self.compute_momentum_distribution(
            q_array, Q_array, pair, nucleus_name, Z, N, density, save)
        
        return n_total[:, 0]
    
    
# Run this script to compute and save momentum distributions
if __name__ == '__main__':
    
    # Specific imports
    import time
    from potentials import Potential
    
    # Default inputs
    kvnn = 6
    kmax, kmid, ntot = 15.0, 3.0, 120
    channels = ('1S0', '3S1')
    generator = 'Wegner'
    lamb = 1.35
    
    pmd = Pair(kvnn, kmax, kmid, ntot, channels, generator, lamb)
    
    # Get momentum values (channel argument doesn't matter here)
    potential = Potential(kvnn, '1S0', kmax, kmid, ntot)
    q_array, _ = potential.load_mesh()
    
    # Set C.o.M. momentum values
    Q_max = 2.0 # Starts to get wonky at Q_max > 2.3 fm^-1
    ntot_Q = 40
    Q_array, _ = gaussian_quadrature_mesh(Q_max, ntot_Q)
    
    # Calculate for Gogny nuclei
    density = 'Gogny'
    nuclei = (
        ('He4', 2, 2), ('Li7', 3, 4), ('Be9', 4, 5), ('C12', 6, 6),
        ('O16', 8, 8), ('Al27', 13, 14), ('Ca40', 20, 20), ('Ca48', 20, 28),
        ('Ti48', 22, 26), ('Fe56', 26, 30), ('Cu63', 29, 34),
        ('Ag107', 47, 60), ('Sn118', 50, 68), ('Ce140', 58, 82),
        ('Ta181', 73, 108), ('Au197', 79, 118), ('Pb208', 82, 126),
        ('U238', 92, 146)
    )
    
    # Compute for Q = 0 too?
    Q0_bool = True
    
    for nucleus in nuclei:
        
        nucleus_name = nucleus[0]
        Z = nucleus[1]
        N = nucleus[2]
        
        t0 = time.time()
        
        for pair in ('pn', 'pp', 'nn'):
            
            # Q > 0
            n_array = pmd.compute_momentum_distribution(
                q_array, Q_array, pair, nucleus_name, Z, N, density, save=True)
        
            # Q = 0
            if Q0_bool:
                n_Q0_array = pmd.compute_momentum_distribution_Q0(
                    q_array, pair, nucleus_name, Z, N, density, save=True)

        t1 = time.time()
        mins = (t1-t0)/60
        
        print(f'Done with {nucleus_name} after {mins:.2f} minutes.')