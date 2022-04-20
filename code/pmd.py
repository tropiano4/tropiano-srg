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
        
        

#------------------------------------------------------------------------------
# File: pmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 31, 2021
# 
# Calculates SRG-evolved pair momentum distributions for nuclei assuming the
# evolved wave function is given by HF treated in LDA. (Note, 'pmd' stands for
# pair momentum distribution.)
#
# Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
# screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!
#
# Revision history:
#   04/26/21 --- Correcting overall factors in front of \delta U and
#                \delta U^2 terms.
#   06/08/21 --- Corrected "I" term to have additional integration over R'.
#   06/22/21 --- Generalizing distribution to n(q, Q) from n(q, 0). Saved old
#                version as pmd_v1.py in Old_codes.
#   06/23/21 --- Speeding up code by switching from loops to np.sum() to do
#                integrations. Saved old version as pmd_v2.py in Old_codes. 
#   07/02/21 --- Added interpolation option for n_{\lambda}(q, Q).
#   07/15/21 --- Switched interpolation functions interp1d and 
#                RectBivariateSpline from cubic to linear since \delta U^2(k,q)
#                was returning negative values. Also, setting Q_max to 2 fm^-1
#                instead of max(kF1)+max(kF2) in file writing function.
#   08/24/21 --- Added SRG block-diagonal option.
#   02/15/22 --- Added optional lambda_init argument to class. This allows one
#                to use an SRG-evolved potential as the initial interaction.
#   02/24/22 --- Added optional kvnn_hard argument to class. This allows one
#                to SRG-evolve the initial potential back to a harder one
#                using transformations of the kvnn_hard potential.
#
#------------------------------------------------------------------------------


# import numpy as np
# from scipy.interpolate import interp1d, RectBivariateSpline
# # Scripts made by A.T.
# # from densities import load_density
# # from figures import figures_functions as ff
# # from misc.integration import gaussian_quadrature_mesh
# # from potentials.vsrg_macos import vnn
# # from srg.srg_unitary_transformation import SRG_unitary_transformation
# from .densities import load_density
# from .figure_labels import replace_periods
# from .integration import gaussian_quadrature_mesh
# from .vnn import Potential
# from .srg_transformation import get_transformation
# from .tools import channel_L_value, coupled_channel


# class pair_momentum_distributions(object):
    
    
#     def __init__(
#             self, kvnn, channels, lamb, kmax, kmid, ntot, generator='Wegner',
#             interp=False, lambda_init=np.inf, kvnn_hard=0):
#         """
#         Evaluates and saves the pp and pn matrix elements of \delta U and
#         \delta U^{\dagger} given the input potential and SRG \lambda.
        
#         Parameters
#         ----------
#         kvnn : int
#             This number specifies the potential.
#         channels : tuple
#             Partial wave channels to include in the calculation.
#         lamb : float
#             SRG evolution parameter \lambda [fm^-1].
#         kmax : float
#             Maximum value in the momentum mesh [fm^-1].
#         kmid : float
#             Mid-point value in the momentum mesh [fm^-1].
#         ntot : int
#             Number of momentum points in mesh.
#         interp : bool, optional
#             Option to use interpolated n_\lambda(q, Q) functions.
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
#         self.channels = channels
#         self.lamb = lamb
#         self.kmax = kmax
#         self.generator = generator

#         if interp == False:
            
#             # Save highest allowed L based on input channels
#             highest_L = 0
#             for channel in channels:
#                 next_L = channel_L_value(channel)
#                 if next_L > highest_L:
#                     highest_L = next_L
        
#             # Load and save momentum arrays
#             k_array, k_weights = Potential(kvnn, '1S0', kmax, kmid, ntot).load_mesh()
#             # Save k_array for writing files
#             self.k_array = k_array
        
#             # For dividing out momenta/weights
#             factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
#             # For coupled-channel matrices
#             factor_array_cc = np.concatenate( (factor_array, factor_array) )

#             # Initialize pp and pn matrix elements
#             deltaU_pp = np.zeros( (ntot, ntot) ) # \delta U linear term
#             deltaU_pn = np.zeros( (ntot, ntot) )
#             deltaU2_pp = np.zeros( (ntot, ntot) ) # \delta U \delta U^\dagger term
#             deltaU2_pn = np.zeros( (ntot, ntot) )
        
#             # Allowed channels for pp (and nn) up through the D-waves
#             pp_channels = ('1S0', '3P0', '3P1', '3P2', '1D2')
        
#             # Loop over channels and evaluate matrix elements
#             for channel in channels:
                
#                 potential = Potential(kvnn, channel, kmax, kmid, ntot)

#                 # Load SRG transformation
                
#                 # If lambda_init = np.inf, then take the initial Hamiltonian
#                 # of kvnn as the starting point
#                 if lambda_init == np.inf:
                    
#                     H_initial = potential.load_hamiltonian()
                    
#                 # Otherwise, this splits into two cases
#                 else:
                    
#                     # If kvnn_hard is zero, take an SRG-evolved Hamiltonian
#                     # corresponding to kvnn and lambda_init as the initial
#                     # Hamiltonian
#                     if kvnn_hard == 0:
                        
#                         # Block-diagonal
#                         if generator == 'Block-diag':
#                             H_initial = potential.load_hamiltonian(
#                                 'srg', generator, 1.0, lambda_init)
                            
#                         # Band-diagonal
#                         else:
#                             H_initial = potential.load_hamiltonian(
#                                 'srg', generator, lambda_init)
                    
#                     # If kvnn_hard is nonzero, SRG-evolve the initial 
#                     # Hamiltonian of kvnn back using transformations from
#                     # kvnn_hard at lambda_initial
#                     else:
                        
#                         potential_hard = Potential(kvnn_hard, channel, kmax, kmid, ntot)
                        
#                         H_hard_initial = potential_hard.load_hamiltonian()
                        
#                         # Block-diagonal
#                         if generator == 'Block-diag':
#                             H_hard_evolved = potential_hard.load_hamiltonian(
#                                 'srg', generator, 1.0, lambda_init)
                            
#                         # Band-diagonal
#                         else:
#                             H_hard_evolved = potential_hard.load_hamiltonian(
#                                 'srg', generator, lambda_init)
                            
#                         # Get SRG transformation from hard potential
#                         U_hard = get_transformation(H_hard_initial,
#                                                     H_hard_evolved)
                        
#                         # Get initial Hamiltonian for kvnn
#                         H_matrix = potential.load_hamiltonian()
                        
#                         # Do inverse transformation on softer Hamiltonian
#                         H_initial = U_hard.T @ H_matrix @ U_hard
                            
                            
#                 # Get evolved Hamiltonian
#                 if generator == 'Block-diag':
#                     # Take \lambda = 1 fm^-1 and set \Lambda_BD = input \lambda
#                     H_evolved = potential.load_hamiltonian('srg', generator,
#                                                            1.0, lamb)
#                 else:
#                     H_evolved = potential.load_hamiltonian('srg', generator,
#                                                            lamb)
#                 # Load U(k, k') [unitless]
#                 U_matrix_unitless = get_transformation(H_initial, H_evolved)

#                 # Isolate 2-body term and convert to fm^3
#                 if coupled_channel(channel):
#                     I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
#                     row, col = np.meshgrid(factor_array_cc, factor_array_cc)
#                 else:
#                     I_matrix_unitless = np.eye(ntot, ntot)
#                     row, col = np.meshgrid(factor_array, factor_array)
#                 delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
#                 delta_U_matrix = delta_U_matrix_unitless / row / col # fm^3
            
#                 # 2J+1 factor
#                 J = int( channel[-1] )
            
#                 # Add to the pp and pn terms
#                 # Coupled-channel
#                 if coupled_channel(channel):
                    
#                     # First L of coupled-channel
#                     # Isospin CG's=1/\sqrt(2) for pn
#                     deltaU_pn += (2*J+1)/2 * delta_U_matrix[:ntot, :ntot]
#                     deltaU2_pn += (2*J+1)/2 * ( \
#                                       delta_U_matrix[:ntot, :ntot]**2 + \
#                                       delta_U_matrix[:ntot, ntot:]**2 )

#                     # Isospin CG's=1 for pp
#                     if channel in pp_channels:
#                         deltaU_pp += (2*J+1) * delta_U_matrix[:ntot, :ntot]
#                         deltaU2_pp += (2*J+1) * ( \
#                                           delta_U_matrix[:ntot, :ntot]**2 \
#                                         + delta_U_matrix[:ntot, ntot:]**2 )
                    
#                     # Decide whether to add second L based on highest allowed
#                     # L value (e.g., 0 + 2 <= 2 meaning we include the 3D1-3D1
#                     # part of the coupled 3S1-3D1 channel if we input D-waves
#                     # in channels)
#                     if channel_L_value(channel) + 2 <= highest_L:
#                         deltaU_pn += (2*J+1)/2 * delta_U_matrix[ntot:, ntot:]
#                         deltaU2_pn += (2*J+1)/2 * ( \
#                                           delta_U_matrix[ntot:, :ntot]**2 + \
#                                           delta_U_matrix[ntot:, ntot:]**2 )
                        
#                         if channel in pp_channels:
#                             deltaU_pp += (2*J+1) * delta_U_matrix[ntot:, ntot:]
#                             deltaU2_pp += (2*J+1) * ( \
#                                               delta_U_matrix[ntot:, :ntot]**2 \
#                                             + delta_U_matrix[ntot:, ntot:]**2 )
            
#                 else:
                
#                     # Isospin CG's=1/\sqrt(2) for pn
#                     deltaU_pn += (2*J+1)/2 * delta_U_matrix
#                     deltaU2_pn += (2*J+1)/2 * delta_U_matrix**2
                
#                     # Isospin CG's=1 for pp
#                     if channel in pp_channels:
#                         deltaU_pp += (2*J+1) * delta_U_matrix
#                         deltaU2_pp += (2*J+1) * delta_U_matrix**2

#             # Interpolate pp and pn < k | \delta U | k >
#             self.deltaU_pp_func = interp1d( k_array, np.diag(deltaU_pp),
#                                             kind='linear', bounds_error=False,
#                                             fill_value='extrapolate' )
#             self.deltaU_pn_func = interp1d( k_array, np.diag(deltaU_pn),
#                                             kind='linear', bounds_error=False,
#                                             fill_value='extrapolate' )
        
#             # Interpolate pp and pn < k | \delta U \delta U^{\dagger} | k' > 
#             self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
#                                                        deltaU2_pp, kx=1, ky=1)
#             self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
#                                                        deltaU2_pn, kx=1, ky=1)


    # def theta_I(self, q_mesh, Q_mesh, kF1_mesh, kF2_mesh):
    #     """
    #     Evaluates angle-average of \theta( kF1(R) - \abs(Q/2 + q) ) x
    #     \theta( kF2(R') - \abs(Q/2 - q) ). This function appears in the I term.

    #     Parameters
    #     ----------
    #     q_mesh : 4-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_mesh : 4-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     kF1_mesh : 4-D ndarray
    #         Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
    #         with respect to R.
    #     kF2_mesh : 4-D ndarray
    #         Fermi momentum [fm^-1] for the second nucleon corresponding to
    #         \tau' with respect to R'.

    #     Returns
    #     -------
    #     theta_mesh : 4-D ndarray
    #         \theta function [unitless] evaluated for each q, Q, kF1(R), and
    #         kF2(R').

    #     """
        
    #     # Initialize 4-D array
    #     theta_mesh = np.zeros( (self.ntot_q, self.ntot_Q, self.ntot_R,
    #                             self.ntot_R) )
        
    #     # Evaluate each boolean case and use these to fill in the theta_mesh

    #     # Case 1: 2q+Q < 2kF1 and 2q+Q < 2kF2
    #     case_1 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
    #              ( 2*q_mesh + Q_mesh <= 2*kF2_mesh )
    #     theta_mesh[case_1] = 1
        
    #     # Q = 0 case simplifies \theta functions
    #     # If Q = 0, skip the following cases
    #     if self.ntot_Q == 1:
    #         return theta_mesh
                
    #     # Case 2: 2q+Q > 2kF1 and 2q+Q > 2kF2 and 4q^2+Q^2 < 2(kF1^2+kF2^2)
    #     case_2 = ( 2*q_mesh + Q_mesh > 2*kF1_mesh ) * \
    #              ( 2*q_mesh + Q_mesh > 2*kF2_mesh ) * \
    #              ( 4*q_mesh**2 + Q_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
    #     theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*q_mesh**2 - \
    #                            Q_mesh**2 ) / (4*q_mesh*Q_mesh) )[case_2]
                            
    #     # Case 3: 2q+Q < 2kF2 and -4 < (4q^2 - 4kF1^2 + Q^2)/(qQ) < 4
    #     case_3 = ( 2*q_mesh + Q_mesh <= 2*kF2_mesh ) * \
    #              ( -4 < (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) ) * \
    #              ( (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) <= 4 )
    #     theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
    #                            (8*q_mesh*Q_mesh) )[case_3]
                
    #     # Case 4: 2q+Q < 2kF1 and -4 < (4q^2 - 4kF2^2 + Q^2)/(qQ) < 4
    #     case_4 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
    #              ( -4 < (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) ) * \
    #              ( (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) <= 4 )
    #     theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
    #                            (8*q_mesh*Q_mesh) )[case_4]
        
    #     # This is a (ntot_q, ntot_Q, ntot_R, ntot_R) size array
    #     return theta_mesh
    
    
    # def theta_deltaU(self, q_mesh, Q_mesh, kF1_mesh, kF2_mesh):
    #     """
    #     Evaluates angle-average of \theta( kF1(R) - \abs(Q/2 + q) ) x
    #     \theta( kF2(R) - \abs(Q/2 - q) ). This function appears in the
    #     \delta U term.

    #     Parameters
    #     ----------
    #     q_mesh : 3-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_mesh : 3-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     kF1_mesh : 3-D ndarray
    #         Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
    #         with respect to R.
    #     kF2_mesh : 3-D ndarray
    #         Fermi momentum [fm^-1] for the second nucleon corresponding to
    #         \tau' with respect to R.

    #     Returns
    #     -------
    #     theta_mesh : 3-D ndarray
    #         \theta function [unitless] evaluated for each q, Q, kF1(R), and
    #         kF2(R).

    #     """
        
    #     # Initialize 3-D array
    #     theta_mesh = np.zeros( (self.ntot_q, self.ntot_Q, self.ntot_R) )
        
    #     # Evaluate each boolean case and use these to fill in the theta_mesh

    #     # Case 1: 2q+Q < 2kF1 and 2q+Q < 2kF2
    #     case_1 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
    #              ( 2*q_mesh + Q_mesh <= 2*kF2_mesh )
    #     theta_mesh[case_1] = 1
        
    #     # Q = 0 case simplifies \theta functions
    #     # If Q = 0, skip the following cases
    #     if self.ntot_Q == 1:
    #         return theta_mesh
                
    #     # Case 2: 2q+Q > 2kF1 and 2q+Q > 2kF2 and 4q^2+Q^2 < 2(kF1^2+kF2^2)
    #     case_2 = ( 2*q_mesh + Q_mesh > 2*kF1_mesh ) * \
    #              ( 2*q_mesh + Q_mesh > 2*kF2_mesh ) * \
    #              ( 4*q_mesh**2 + Q_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
    #     theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*q_mesh**2 - \
    #                            Q_mesh**2 ) / (4*q_mesh*Q_mesh) )[case_2]
                            
    #     # Case 3: 2q+Q < 2kF2 and -4 < (4q^2 - 4kF1^2 + Q^2)/(qQ) < 4
    #     case_3 = ( 2*q_mesh + Q_mesh <= 2*kF2_mesh ) * \
    #              ( -4 < (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) ) * \
    #              ( (4*q_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) <= 4 )
    #     theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
    #                            (8*q_mesh*Q_mesh) )[case_3]
                
    #     # Case 4: 2q+Q < 2kF1 and -4 < (4q^2 - 4kF2^2 + Q^2)/(qQ) < 4
    #     case_4 = ( 2*q_mesh + Q_mesh <= 2*kF1_mesh ) * \
    #              ( -4 < (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) ) * \
    #              ( (4*q_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
    #              (q_mesh*Q_mesh) <= 4 )
    #     theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (Q_mesh-2*q_mesh)**2 ) / \
    #                            (8*q_mesh*Q_mesh) )[case_4]
        
    #     # This is a (ntot_q, ntot_Q, ntot_R) size array
    #     return theta_mesh
    
    
    # def theta_deltaU2(self, Q_mesh, kF1_mesh, kF2_mesh, k_mesh):
    #     """
    #     Evaluates angle-average of \theta( kF1(R) - \abs(Q/2 + k) ) x
    #     \theta( kF2(R) - \abs(Q/2 - k) ). This function appears in the
    #     \delta U \delta U^\dagger term.

    #     Parameters
    #     ----------
    #     Q_mesh : 4-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     kF1_mesh : 4-D ndarray
    #         Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
    #         with respect to R.
    #     kF2_mesh : 4-D ndarray
    #         Fermi momentum [fm^-1] for the second nucleon corresponding to
    #         \tau' with respect to R.
    #     k_mesh : 4-D ndarray
    #         Relative momentum values [fm^-1].

    #     Returns
    #     -------
    #     theta_mesh : 4-D ndarray
    #         \theta function [unitless] evaluated for each q, Q, kF1(R), kF2(R),
    #         and k. Note, this function does not depend on q but we need the
    #         array to match the size of \delta U \delta U^\dagger matrix
    #         elements.

    #     """
        
    #     # Initialize 4-D array
    #     theta_mesh = np.zeros( (self.ntot_q, self.ntot_Q, self.ntot_R,
    #                             self.ntot_k) )
        
    #     # Evaluate each boolean case and use these to fill in the theta_mesh

    #     # Case 1: 2k+Q < 2kF1 and 2k+Q < 2kF2
    #     case_1 = ( 2*k_mesh + Q_mesh <= 2*kF1_mesh ) * \
    #              ( 2*k_mesh + Q_mesh <= 2*kF2_mesh )
    #     theta_mesh[case_1] = 1
        
    #     # Q = 0 case simplifies \theta functions
    #     # If Q = 0, skip the following cases
    #     if self.ntot_Q == 1:
    #         return theta_mesh
                
    #     # Case 2: 2k+Q > 2kF1 and 2k+Q > 2kF2 and 4k^2+Q^2 < 2(kF1^2+kF2^2)
    #     case_2 = ( 2*k_mesh + Q_mesh > 2*kF1_mesh ) * \
    #              ( 2*k_mesh + Q_mesh > 2*kF2_mesh ) * \
    #              ( 4*k_mesh**2 + Q_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
    #     theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*k_mesh**2 - \
    #                            Q_mesh**2 ) / (4*k_mesh*Q_mesh) )[case_2]
                            
    #     # Case 3: 2k+Q < 2kF2 and -4 < (4k^2 - 4kF1^2 + Q^2)/(kQ) < 4
    #     case_3 = ( 2*k_mesh + Q_mesh <= 2*kF2_mesh ) * \
    #              ( -4 < (4*k_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
    #              (k_mesh*Q_mesh) ) * \
    #              ( (4*k_mesh**2-4*kF1_mesh**2+Q_mesh**2) / \
    #              (k_mesh*Q_mesh) <= 4 )
    #     theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (Q_mesh-2*k_mesh)**2 ) / \
    #                            (8*k_mesh*Q_mesh) )[case_3]
                
    #     # Case 4: 2k+Q < 2kF1 and -4 < (4k^2 - 4kF2^2 + Q^2)/(kQ) < 4
    #     case_4 = ( 2*k_mesh + Q_mesh <= 2*kF1_mesh ) * \
    #              ( -4 < (4*k_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
    #              (k_mesh*Q_mesh) ) * \
    #              ( (4*k_mesh**2-4*kF2_mesh**2+Q_mesh**2) / \
    #              (k_mesh*Q_mesh) <= 4 )
    #     theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (Q_mesh-2*k_mesh)**2 ) / \
    #                            (8*k_mesh*Q_mesh) )[case_4]
        
    #     # This is a (ntot_q, ntot_Q, ntot_R, ntot_k) size array
    #     return theta_mesh
    
    
    # def n_I(self, q_array, Q_array, R_array, dR, kF1_array, kF2_array):
    #     """
    #     Evaluates the I term in U n(q) U^\dagger ~ I n(q) I.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_array : 1-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     kF1_array : 1-D ndarray
    #         Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
    #         with respect to R.
    #     kF2_array : 1-D ndarray
    #         Fermi momentum [fm^-1] for the second nucleon corresponding to
    #         \tau' with respect to R'.

    #     Returns
    #     -------
    #     output : 2-D ndarray
    #         Momentum distribution [fm^6] from I term as a function of q and Q.

    #     """
        
    #     # Initialize 4-D meshgrids (q, Q, R, R')
    #     q_mesh, Q_mesh, R_mesh, Rp_mesh = np.meshgrid(q_array, Q_array,
    #                                                   R_array, R_array,
    #                                                   indexing='ij')
    #     # Get 4-D kF1(R) and kF2(R') meshes
    #     _, _, kF1_mesh, kF2_mesh = np.meshgrid(q_array, Q_array, kF1_array,
    #                                            kF2_array, indexing='ij')
        
    #     # Evaluate angle-average of \theta-functions in I term
    #     theta_mesh = self.theta_I(q_mesh, Q_mesh, kF1_mesh, kF2_mesh)
        
    #     # Calculate R' integrand (ntot_q, ntot_Q, ntot_R, ntot_R)
    #     integrand_Rp = theta_mesh * Rp_mesh**2 * dR * R_mesh**2 * dR
        
    #     # Integrate over R' leaving R integrand (ntot_q, ntot_Q, ntot_R)
    #     integrand_R = 4*np.pi * np.sum(integrand_Rp, axis=-1)
        
    #     # Integrate over R
    #     # This is a (ntot_q, ntot_Q) size array
    #     # Factor of 2 is overall factor
    #     return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    
    # def n_deltaU(self, q_array, Q_array, R_array, dR, kF1_array, kF2_array,
    #              pair):
    #     """
    #     Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
    #     we are combining \delta U and \delta U^\dagger, hence the factor of 2.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_array : 1-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     kF1_array : 1-D ndarray
    #         Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
    #         with respect to R.
    #     kF2_array : 1-D ndarray
    #         Fermi momentum [fm^-1] for the second nucleon corresponding to
    #         \tau' with respect to R.
    #     pair : str
    #         Type of pair momentum distribution ('pp' or 'pn'). This will tell
    #         the function to either set isospin CGs to 1 or 1/\sqrt(2) (no need
    #         to specify 'nn' or 'np').

    #     Returns
    #     -------
    #     output : 2-D ndarray
    #         Momentum distribution [fm^6] from \delta U term as a function of q
    #         and Q.

    #     """
        
    #     # Initialize 3-D meshgrids (q, Q, R)
    #     q_mesh, Q_mesh, R_mesh = np.meshgrid(q_array, Q_array, R_array,
    #                                          indexing='ij')
        
    #     # Get 3-D kF1(R) and kF2(R) meshes
    #     _, _, kF1_mesh = np.meshgrid(q_array, Q_array, kF1_array,
    #                                  indexing='ij')
    #     _, _, kF2_mesh = np.meshgrid(q_array, Q_array, kF2_array,
    #                                  indexing='ij')
        
    #     # Evaluate angle-average of \theta-functions in \delta U term
    #     theta_mesh = self.theta_deltaU(q_mesh, Q_mesh, kF1_mesh, kF2_mesh)
        
    #     # Evaluate < q | \delta U | q >
    #     if pair == 'pp':
    #         deltaU_mesh = self.deltaU_pp_func(q_mesh)
    #     elif pair == 'pn':
    #         deltaU_mesh = self.deltaU_pn_func(q_mesh)

    #     # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
    #     # combine to give 1
    #     # Factor of 2 from \delta U + \delta U^\dagger
    #     # 2/\pi for two | k_vec > -> | k J L S ... > changes
    #     # 1/(4\pi) for averaging over \int d\Omega_q
    #     deltaU_factor =  2 * 2/np.pi * (2*np.pi)**3/(4*np.pi)
        
    #     # Calculate R integrand (ntot_q, ntot_Q, ntot_R)
    #     integrand_R = deltaU_factor * deltaU_mesh * theta_mesh * R_mesh**2 * dR
        
    #     # Integrate over R
    #     # This is a (ntot_q, ntot_Q) size array
    #     return 4*np.pi * np.sum(integrand_R, axis=-1)
        
    
    # def n_deltaU2(self, q_array, Q_array, R_array, dR, kF1_array, kF2_array,
    #               pair):
    #     """
    #     Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_array : 1-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     kF1_array : 1-D ndarray
    #         Fermi momentum [fm^-1] for the first nucleon corresponding to \tau
    #         with respect to R.
    #     kF2_array : 1-D ndarray
    #         Fermi momentum [fm^-1] for the second nucleon corresponding to
    #         \tau' with respect to R.
    #     pair : str
    #         Type of pair momentum distribution ('pp' or 'pn'). This will tell
    #         the function to either set isospin CGs to 1 or 1/\sqrt(2) (no need
    #         to specify 'nn' or 'np').

    #     Returns
    #     -------
    #     output : 2-D ndarray
    #         Momentum distribution [fm^6] from \delta U \delta U^\dagger term as
    #         a function of q and Q.

    #     """
        
    #     # Set number of k points for integration over k
    #     self.ntot_k = 40
        
    #     # Initialize 4-D meshgrids (q, Q, R, k) with k values equal to 0
    #     q_mesh, Q_mesh, R_mesh, k_mesh = np.meshgrid(q_array, Q_array, R_array,
    #                                                  np.zeros(self.ntot_k),
    #                                                  indexing='ij')
        
    #     # Get 4-D kF1(R) and kF2(R) meshes and initialize k weights mesh
    #     _, _, kF1_mesh, dk_mesh = np.meshgrid(q_array, Q_array, kF1_array,
    #                                           np.zeros(self.ntot_k),
    #                                           indexing='ij')
    #     _, _, kF2_mesh, _ = np.meshgrid(q_array, Q_array, kF2_array,
    #                                     np.zeros(self.ntot_k), indexing='ij')

    #     # Loop over q, Q, and R to find limits of k integration and then create
    #     # k_array using Gaussian quadrature
    #     for iQ, Q in enumerate(Q_array):
    #         for iR, R in enumerate(R_array):
                
    #             kF1, kF2 = kF1_array[iR], kF2_array[iR]
                
    #             # Minimum kF value
    #             kF_min = min(kF1, kF2)
                
    #             # Lower limit of integration
    #             k_min = max(Q/2 - kF_min, 0)
                
    #             # Upper limit of integration
    #             if Q**2/4 < (kF1**2 + kF2**2)/2:
    #                 k_max = min( np.sqrt( ( kF1**2 + kF2**2 )/2 - Q**2/4 ),
    #                              kF_min + Q/2 )
    #             else:
    #                 k_max = kF_min + Q/2
                    
    #             # Get Gaussian quadrature mesh
    #             k_array, k_weights = gaussian_quadrature_mesh(k_max,
    #                                                           self.ntot_k,
    #                                                           xmin=k_min)
                
    #             # Fill in k_mesh and dk_mesh given the specific k_array
    #             for iq in range(self.ntot_q):
    #                 k_mesh[iq, iQ, iR, :] = k_array
    #                 dk_mesh[iq, iQ, iR, :] = k_weights
                     
    #     # Evaluate angle-average of \theta-functions in \delta U^2 term
    #     theta_mesh = self.theta_deltaU2(Q_mesh, kF1_mesh, kF2_mesh, k_mesh)
        
    #     # Evaluate < k | \delta U | q > < q | \delta U^\dagger | k >
    #     if pair == 'pp':
    #         deltaU2_mesh = self.deltaU2_pp_func.ev(k_mesh, q_mesh)
    #     elif pair == 'pn':
    #         deltaU2_mesh = self.deltaU2_pn_func.ev(k_mesh, q_mesh)

    #     # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
    #     # combine to give 1
    #     # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
    #     # 1/(4\pi) for averaging over \int d\Omega_q
    #     deltaU2_factor =  (2/np.pi)**2 * (2*np.pi)**3/(4*np.pi)
        
    #     # Calculate k integrand (ntot_q, ntot_Q, ntot_R, ntot_k)
    #     integrand_k = deltaU2_factor * deltaU2_mesh * theta_mesh * \
    #                   k_mesh**2 * dk_mesh * R_mesh**2 * dR
        
    #     # Integrate over k leaving R integrand (ntot_q, ntot_Q, ntot_R)
    #     integrand_R = np.sum(integrand_k, axis=-1)

    #     # Integrate over R
    #     # This is a (ntot_q, ntot_Q) size array
    #     return 4*np.pi * np.sum(integrand_R, axis=-1)


    # def n_total(self, q_array, Q_array, R_array, dR, rho_1_array,
    #             rho_2_array=np.empty(0)):
    #     """
    #     Pair momentum distribution where the nucleonic densities specify the
    #     nucleus and distribution type (e.g., O16 and pn).

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_array : 1-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     rho_1_array : 1-D ndarray
    #         Densities as a function of R [fm^-3] for the first nucleon
    #         corresponding to \tau.
    #     rho_2_array : 1-D ndarray, optional
    #         Densities as a function of R [fm^-3] for the second nucleon
    #         corresponding to \tau'. If an empty array is input, the function
    #         assumes a proton-proton (or neutron-neutron) pair momentum
    #         distribution relying only on rho_1_array.

    #     Returns
    #     -------
    #     n_total : 2-D ndarray
    #         Pair momentum distribution [fm^6] for each q and Q.

    #     """
        
    #     # Save lengths of q_array, Q_array, and R_array
    #     self.ntot_q = len(q_array)
    #     self.ntot_Q = len(Q_array)
    #     self.ntot_R = len(R_array)

    #     # Evaluate kF values at each point in R_array
    #     kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
    #     # Calculate kF2 array if pn or np distribution
    #     if rho_2_array.any():
    #         kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
    #         pair = 'pn'
    #     # Otherwise set kF2=kF1 where previous functions will give pp or nn
    #     # distributions (depending on kF_1)
    #     else:
    #         kF2_array = kF1_array
    #         pair = 'pp'
            
    #     # Get each contribution with respect to q and Q
    #     n_I = self.n_I(q_array, Q_array, R_array, dR, kF1_array, kF2_array)
    #     n_deltaU = self.n_deltaU(q_array, Q_array, R_array, dR, kF1_array,
    #                              kF2_array, pair)
    #     n_deltaU2 = self.n_deltaU2(q_array, Q_array, R_array, dR, kF1_array,
    #                                kF2_array, pair)
        
    #     # Return total (ntot_q, ntot_Q)
    #     return n_I + n_deltaU + n_deltaU2
    
    
    # def n_contributions(self, q_array, Q_array, R_array, dR, rho_1_array,
    #                     rho_2_array=np.empty(0)):
    #     """
    #     Contributions to the pair momentum distribution where the nucleonic
    #     densities specify the nucleus and distribution type (e.g., O16 and pn).
    #     This function isolates the I, \delta U, and \delta U \delta U^\dagger
    #     terms.

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Relative momentum values [fm^-1].
    #     Q_array : 1-D ndarray
    #         C.o.M. momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     rho_1_array : 1-D ndarray
    #         Densities as a function of R [fm^-3] for the first nucleon
    #         corresponding to \tau.
    #     rho_2_array : 1-D ndarray, optional
    #         Densities as a function of R [fm^-3] for the second nucleon
    #         corresponding to \tau'. If an empty array is input, the function
    #         assumes a proton-proton (or neutron-neutron) pair momentum
    #         distribution relying only on rho_1_array.

    #     Returns
    #     -------
    #     n_contributions : tuple
    #         Tuple of 2-D ndarrays corresponding to the I, \delta U, and
    #         \delta U^\dagger terms of the pair momentum distribution [fm^6] for
    #         each q and Q.

    #     """
        
    #     # Save lengths of q_array, Q_array, and R_array
    #     self.ntot_q = len(q_array)
    #     self.ntot_Q = len(Q_array)
    #     self.ntot_R = len(R_array)

    #     # Evaluate kF values at each point in R_array
    #     kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
    #     # Calculate kF2 array if pn or np distribution
    #     if rho_2_array.any():
    #         kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
    #         pair = 'pn'
    #     # Otherwise set kF2=kF1 where previous functions will give pp or nn
    #     # distributions (depending on kF_1)
    #     else:
    #         kF2_array = kF1_array
    #         pair = 'pp'
            
    #     # Get each contribution with respect to q and Q
    #     n_I = self.n_I(q_array, Q_array, R_array, dR, kF1_array, kF2_array)
    #     n_deltaU = self.n_deltaU(q_array, Q_array, R_array, dR, kF1_array,
    #                              kF2_array, pair)
    #     n_deltaU2 = self.n_deltaU2(q_array, Q_array, R_array, dR, kF1_array,
    #                                kF2_array, pair)
        
    #     # Return tuple of contributions ( (ntot_q, ntot_Q), ... )
    #     return n_I, n_deltaU, n_deltaU2
    
    
    # def write_file(self, nucleus, pair, Z, N, edf='SLY4'):
    #     """
    #     Write pair momentum distribution file for interpolation purposes.
    #     Split things into total, I, \delta U, and \delta U^2 contributions.

    #     Parameters
    #     ----------
    #     nucleus : str
    #         Specify nucleus (e.g., 'O16', 'Ca40', etc.)
    #     pair : str
    #         Specify 'pp', 'pn', or 'nn'.
    #     Z : int
    #         Proton number of the nucleus.
    #     N : int
    #         Neutron number of the nucleus.
    #     edf : str, optional
    #         Name of EDF (e.g., 'SLY4').
            
    #     Notes
    #     -----
    #     'pn' means \tau = +1/2 and \tau' = -1/2 and does NOT account for the
    #     opposite case (np). Many references combine these two. Multiply by 2 to
    #     match those references.

    #     """
        
    #     # Get relative momentum values
    #     q_array = self.k_array
        
    #     # Directory for distributions data
    #     data_directory = 'data/pmd/kvnn_%d/%s' % (self.kvnn, edf)
        
    #     # Create file name
    #     file_name = '%s_%s_channels' % (nucleus, pair)
    #     # Add each channel to file name
    #     for channel in self.channels:
    #         file_name += '_%s' % channel
    #     if self.generator == 'Block-diag':
    #         file_name += '_LambdaBD_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
    #     else:
    #         file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
    #     file_name = replace_periods(file_name) + '.dat'
        
    #     # Load R values and nucleonic densities (the R_array's are the same)
    #     R_array, rho_p_array = load_density(nucleus, 'proton', Z, N, edf)
    #     if edf == 'AV18' and Z == N: # e.g., AV18 He4 densities
    #         rho_n_array = rho_p_array
    #     else: # e.g., AV18 He8 densities
    #         R_array, rho_n_array = load_density(nucleus, 'neutron', Z, N, edf)
    #     dR = R_array[2] - R_array[1] # Assuming linear spacing
        
    #     # Set C.o.M. momentum values
    #     Q_max = 2.0 # Starts to get wonky at Q_max > 2.3 fm^-1
    #     ntot_Q = 40
    #     Q_array, Q_weights = gaussian_quadrature_mesh(Q_max, ntot_Q)

    #     # Calculate n_\lambda^\tau(q, Q) for each q in q_array and Q in Q_array
    #     if pair == 'pp':
    #         n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
    #                                 q_array, Q_array, R_array, dR, rho_p_array)
    #     elif pair == 'nn':
    #         n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
    #                                 q_array, Q_array, R_array, dR, rho_n_array)
    #     else: # pn pair
    #         n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
    #                    q_array, Q_array, R_array, dR, rho_p_array, rho_n_array)
            
    #     # Total momentum distribution
    #     n_total_array = n_I_array + n_delU_array + n_delU2_array
    
    #     # Open file and write header where we allocate roughly 18 centered
    #     # spaces for each label
    #     f = open(data_directory + '/' + file_name, 'w')
    #     header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}{:^18s}'.format('q',
    #              'Q', 'total', '1', '\delta U', '\delta U^2')
    #     f.write(header + '\n')
    
    #     # Loop over momenta q and Q
    #     for iq, q in enumerate(q_array):
    #         for iQ, Q in enumerate(Q_array):

    #             # Write to data file following the format from the header
    #             line = '{:^18.6f}{:^18.6f}{:^18.6e}{:^18.6e}{:^18.6e}{:^18.6e}' \
    #                    .format( q, Q, n_total_array[iq, iQ], n_I_array[iq, iQ],
    #                             n_delU_array[iq, iQ], n_delU2_array[iq, iQ] )
    #             f.write('\n' + line)

    #     # Close file
    #     f.close()
        
        
    # def n_lambda_interp(self, nucleus, pair, Z, N, edf='SLY4'):
    #     """
    #     Interpolate the pair momentum distribution for the specified file.

    #     Parameters
    #     ----------
    #     nucleus : str
    #         Specify nucleus (e.g., 'O16', 'Ca40', etc.)
    #     pair : str
    #         Specify 'pp', 'pn', or 'nn'.
    #     Z : int
    #         Proton number of the nucleus.
    #     N : int
    #         Neutron number of the nucleus.
    #     edf : str, optional
    #         Name of EDF (e.g., 'SLY4').
            
    #     Notes
    #     -----
    #     'pn' means \tau = +1/2 and \tau' = -1/2 and does NOT account for the
    #     opposite case (np). Many references combine these two. Multiply by 2 to
    #     match those references.

    #     """
        
    #     # Directory for distributions data
    #     data_directory = 'data/pmd/kvnn_%d/%s' % (self.kvnn, edf)
        
    #     # Get file name
    #     file_name = '%s_%s_channels' % (nucleus, pair)
    #     # Add each channel to file name
    #     for channel in self.channels:
    #         file_name += '_%s' % channel
    #     if self.generator == 'Block-diag':
    #         file_name += '_LambdaBD_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
    #     else:
    #         file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
    #     file_name = replace_periods(file_name) + '.dat'
        
    #     # Load data which includes all contributions to n_\lambda(q)
    #     data = np.loadtxt(data_directory + '/' + file_name)

    #     # Get C.o.M. and relative momentum values
    #     ntot_Q = 40
    #     ntot_q = self.ntot
    #     q_array = np.reshape( data[:, 0], (ntot_q, ntot_Q) )[:, 0]
    #     Q_array = np.reshape( data[:, 1], (ntot_q, ntot_Q) )[0, :]
        
    #     # Split data into 1-D arrays for each column
    #     # Total distribution
    #     n_total_array = np.reshape( data[:, 2], (ntot_q, ntot_Q) )
    #     # 1 term
    #     n_I_array = np.reshape( data[:, 3], (ntot_q, ntot_Q) )
    #     # \delta U term
    #     n_delU_array = np.reshape( data[:, 4], (ntot_q, ntot_Q) )
    #     # \delta U^2 term
    #     n_delU2_array = np.reshape( data[:, 5], (ntot_q, ntot_Q) )
        
    #     # Interpolate each array
    #     n_total_func = RectBivariateSpline(q_array, Q_array, n_total_array,
    #                                        kx=1, ky=1)
    #     n_I_func = RectBivariateSpline(q_array, Q_array, n_I_array, kx=1, ky=1)
    #     n_delU_func = RectBivariateSpline(q_array, Q_array, n_delU_array, kx=1,
    #                                       ky=1)
    #     n_delU2_func = RectBivariateSpline(q_array, Q_array, n_delU2_array,
    #                                        kx=1, ky=1)
        
    #     # Return all contributions with total first
    #     # Note, these are functions of q and Q
    #     return n_total_func, n_I_func, n_delU_func, n_delU2_func
    
    
    # def n_Q0(self, q_array, R_array, dR, rho_1_array, rho_2_array=np.empty(0)):
    #     """
    #     Pair momentum distribution evaluated at Q = 0 where the nucleonic
    #     densities specify the nucleus and distribution type (e.g., O16 and pn).

    #     Parameters
    #     ----------
    #     q_array : 1-D ndarray
    #         Relative momentum values [fm^-1].
    #     R_array : 1-D ndarray
    #         C.o.M. coordinates [fm].
    #     dR : float
    #         C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
    #     rho_1_array : 1-D ndarray
    #         Densities as a function of R [fm^-3] for the first nucleon
    #         corresponding to \tau.
    #     rho_2_array : 1-D ndarray, optional
    #         Densities as a function of R [fm^-3] for the second nucleon
    #         corresponding to \tau'. If an empty array is input, the function
    #         assumes a proton-proton (or neutron-neutron) pair momentum
    #         distribution relying only on rho_1_array.

    #     Returns
    #     -------
    #     n_total : 1-D ndarray
    #         Pair momentum distribution [fm^6] for each q.

    #     """
        
    #     # Create Q_array with a single zero entry
    #     # \theta functions simplify to give 1 or 0 in this case
    #     Q_array = np.array([0.0])
        
    #     # Save lengths of q_array, Q_array, and R_array
    #     self.ntot_q = len(q_array)
    #     self.ntot_Q = len(Q_array)
    #     self.ntot_R = len(R_array)

    #     # Evaluate kF values at each point in R_array
    #     kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
    #     # Calculate kF2 array if pn or np distribution
    #     if rho_2_array.any():
    #         kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
    #         distribution = 'pn'
    #     # Otherwise set kF2=kF1 where previous functions will give pp or nn
    #     # distributions (depending on kF_1)
    #     else:
    #         kF2_array = kF1_array
    #         distribution = 'pp'
            
    #     # Get each contribution with respect to q and Q
    #     n_I = self.n_I(q_array, Q_array, R_array, dR, kF1_array, kF2_array)
    #     n_deltaU = self.n_deltaU(q_array, Q_array, R_array, dR, kF1_array,
    #                              kF2_array, distribution)
    #     n_deltaU2 = self.n_deltaU2(q_array, Q_array, R_array, dR, kF1_array,
    #                                kF2_array, distribution)
        
    #     # Return total (ntot_q, 1)
    #     return ( n_I + n_deltaU + n_deltaU2 )[:, 0]