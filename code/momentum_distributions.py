#!/usr/bin/env python3

"""
File: momentum_distributions.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 17, 2022

Handles momentum distribution codes and data. Momentum distributions are
computed using sub-classes from snmd.py, pmd.py, and dmd.py, and SRG-
transformed potentials from potentials.py. Momentum distribution data is
stored in data/momentum_distributions. These codes also rely on nucleonic
densities taken from external codes or data.
     
Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!

Last update: April 8, 2022

"""

# To-do: Probably want a function that gets k_array, k_weights independent of
# the channel.
# To-do: Make sure description above makes sense.

# Python imports
import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline

# Imports from A.T. codes
from modules.tools import channel_L_value, coupled_channel
from potentials import Potential
from srg import get_transformation


class MomentumDistribution:

    def __init__(self, kvnn, kmax, kmid, ntot):
        """
        Save the inputs of the potential excluding the channels argument.

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

        """

        self.kvnn = kvnn
        self.kmax = kmax
        self.kmid = kmid
        self.ntot = ntot

    def get_hamiltonians(
            self, channel, generator, lamb, lambda_initial=np.inf, kvnn_inv=0,
            delta_lambda=np.inf):
        """
        Get the initial and SRG-evolved Hamiltonians.
        
        Parameters
        ----------
        channel : str
            The partial wave channel (e.g. '1S0').
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
            
        Returns
        -------
        H_initial : 2-D ndarray
            Initial Hamiltonian matrix [MeV].
        H_evolved : 2-D ndarray
            Evolved Hamiltonian matrix [MeV].
            
        """

        # Set potential
        potential = Potential(self.kvnn, channel, self.kmax, self.kmid,
                              self.ntot)

        # Standard initial Hamiltonian
        if lambda_initial == np.inf:
            H_initial = potential.load_hamiltonian()

        # Forward-SRG-evolved Hamiltonian as starting point
        else:
            if generator == 'Block-diag':
                H_initial = potential.load_hamiltonian(
                    'srg', generator, 1.0, lambda_bd=lambda_initial)
            else:
                H_initial = potential.load_hamiltonian('srg', generator,
                                                       lambda_initial)

        # Backwards-SRG-evolved Hamiltonian as starting point
        if kvnn_inv:
            potential_hard = Potential(kvnn_inv, channel, self.kmax, self.kmid,
                                       self.ntot)
            H_hard_initial = potential_hard.load_hamiltonian()
            if generator == 'Block-diag':
                H_hard_evolved = potential_hard.load_hamiltonian(
                    'srg', generator, 1.0, lambda_bd=lambda_initial)
            else:
                H_hard_evolved = potential_hard.load_hamiltonian(
                    'srg', generator, lambda_initial)
            # Get SRG transformation for inverse transformation
            U_hard = get_transformation(H_hard_initial, H_hard_evolved)
            # Do inverse transformation on initial Hamiltonian
            H_initial = U_hard.T @ H_initial @ U_hard

        # Get SRG-evolved Hamiltonian at \lambda=lamb
        if generator == 'Block-diag':
            H_evolved = potential.load_hamiltonian('srg', generator, 1.0,
                                                   lambda_bd=lamb)
        else:
            H_evolved = potential.load_hamiltonian('srg', generator, lamb)

        return H_initial, H_evolved
    
    def get_deltaU_matrix_element(self, channel, delta_U_matrix):
        """
        Manipulates the \delta U(k,k') to return contributions up to the
        highest partial wave corresponding to the input channel.
        
        Parameters
        ----------
        channel : str
            The partial wave channel (e.g. '1S0').
        delta_U_matrix : 2-D ndarray
            \delta U matrix given some partial wave channel [fm^3].
        
        Returns
        -------
        deltaU : 2-D ndarray
            \delta U matrix [fm^3].
        deltaU_squared : 2-D
            \delta U \delta U^{\dagger} matrix [fm^6].
        
        """
        
        # This case corresponds to coupled-channel partial waves
        if coupled_channel(channel):
            
            # First L of coupled-channel
            deltaU = delta_U_matrix[:self.ntot, :self.ntot]
            deltaU_squared = (delta_U_matrix[:self.ntot, :self.ntot]**2
                              + delta_U_matrix[:self.ntot, self.ntot:]**2)

            # Decide whether to add second L based on highest allowed L value
            # (e.g., include the 3D1-3D1 part of the coupled 3S1-3D1 channel
            # if we input D-waves in channels)
            if channel_L_value(channel) + 2 <= self.highest_L:
                deltaU += delta_U_matrix[self.ntot:, self.ntot:]
                deltaU_squared += (delta_U_matrix[self.ntot:, :self.ntot]**2
                                   + delta_U_matrix[self.ntot:, self.ntot:]**2)

        else:

            deltaU = delta_U_matrix
            deltaU_squared = delta_U_matrix**2

        return deltaU, deltaU_squared

    def save_deltaU_funcs(
            self, channels, generator, lamb, lambda_initial=np.inf, kvnn_inv=0,
            delta_lambda=np.inf):
        """
        Save the function \delta U(k,k') and \delta U(k,k')^2 summing over all
        partial wave channels. We must define functions with \tau=\tau' and
        \tau=-\tau' separately since they have different isospin CG's.

        Parameters
        ----------
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

        Notes
        -----
        We are taking \lambda=1.0 fm^-1 for block-diagonal decoupling and
        assuming parameters lamb, lambda_initial, and delta_lambda correspond
        to \Lambda_BD.
        
        """
        
        # Save highest allowed L based on input channels
        highest_L = 0
        for channel in channels:
            next_L = channel_L_value(channel)
            if next_L > highest_L:
                highest_L = next_L
        self.highest_L = highest_L

        # Allowed channels for pp (and nn) up through the D-waves
        pp_channels = ('1S0', '3P0', '3P1', '3P2', '1D2')
        
        # Get momentum mesh (channel argument doesn't matter here)
        k_array, k_weights = Potential(
            self.kvnn, '1S0', self.kmax, self.kmid, self.ntot).load_mesh()

        # For dividing out momenta/weights
        factor_array = np.sqrt((2*k_weights)/np.pi) * k_array
        # For coupled-channel matrices
        factor_array_cc = np.concatenate((factor_array, factor_array))

        # Initialize \delta U linear term
        deltaU_pp = np.zeros((self.ntot, self.ntot))
        deltaU_pn = np.zeros((self.ntot, self.ntot))
        # Initialize \delta U \delta U^\dagger
        deltaU2_pp = np.zeros((self.ntot, self.ntot))
        deltaU2_pn = np.zeros((self.ntot, self.ntot))

        # Loop over channels and evaluate matrix elements
        for channel in channels:

            # Get initial and evolved Hamiltonians
            H_initial, H_evolved = self.get_hamiltonians(
                channel, generator, lamb, lambda_initial, kvnn_inv,
                delta_lambda)

            # Get SRG transformation U(k, k') [unitless]
            U_matrix_unitless = get_transformation(H_initial, H_evolved)

            # Isolate 2-body term and convert to fm^3
            if coupled_channel(channel):
                I_matrix_unitless = np.eye(2*self.ntot, 2*self.ntot)
                row, col = np.meshgrid(factor_array_cc, factor_array_cc)
            else:
                I_matrix_unitless = np.eye(self.ntot, self.ntot)
                row, col = np.meshgrid(factor_array, factor_array)
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            delta_U_matrix = delta_U_matrix_unitless / row / col  # fm^3
            
            # 2J+1 factor
            J = int(channel[-1])
            
            # Get matrix elements up to highest input partial wave channel
            deltaU, deltaU2 = self.get_deltaU_matrix_element(channel,
                                                             delta_U_matrix)
            
            # Isospin CG's=1/\sqrt(2) for pn
            deltaU_pn += (2*J+1)/2 * deltaU
            deltaU2_pn += (2*J+1)/2 * deltaU2
            
            # Isospin CG's=1 for pp
            if channel in pp_channels:
                deltaU_pp += (2*J+1) * deltaU
                deltaU2_pp += (2*J+1) * deltaU2

        # Interpolate pp and pn \delta U(k,k)
        self.deltaU_pp_func = interp1d(
            k_array, np.diag(deltaU_pp), kind='linear', bounds_error=False,
            fill_value='extrapolate')
        self.deltaU_pn_func = interp1d(
            k_array, np.diag(deltaU_pn), kind='linear', bounds_error=False,
            fill_value='extrapolate')

        # Interpolate pp and pn \delta U^2(k,k') 
        self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pp, kx=1, ky=1)
        self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pn, kx=1, ky=1)
        
    def save_deuteron_deltaU_funcs(
            self, generator, lamb, lambda_initial=np.inf, kvnn_inv=0,
            delta_lambda=np.inf):
        """
        Save the function \delta U(k,k') and \delta U(k,k')^2 for the deuteron
        momentum distribution (meaing 3S1-3D1 only). No 2J+1 factor since we
        are fixing M_J for deuteron.

        Parameters
        ----------
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

        Notes
        -----
        We are taking \lambda=1.0 fm^-1 for block-diagonal decoupling and
        assuming parameters lamb, lambda_initial, and delta_lambda correspond
        to \Lambda_BD.
        
        """
        
        # Channel is 3S1-3D1 for deuteron
        channel = '3S1'

        # Get momentum mesh (channel argument doesn't matter here)
        k_array, k_weights = Potential(
            self.kvnn, channel, self.kmax, self.kmid, self.ntot).load_mesh()

        # For dividing out momenta/weights
        factor_array = np.sqrt((2*k_weights)/np.pi) * k_array
        factor_array_cc = np.concatenate((factor_array, factor_array))

        # Get initial and evolved Hamiltonians
        H_initial, H_evolved = self.get_hamiltonians(
            channel, generator, lamb, lambda_initial, kvnn_inv, delta_lambda)

        # Get SRG transformation U(k, k') [unitless]
        U_matrix_unitless = get_transformation(H_initial, H_evolved)

        # Isolate 2-body term and convert to fm^3
        I_matrix_unitless = np.eye(2*self.ntot, 2*self.ntot)
        row, col = np.meshgrid(factor_array_cc, factor_array_cc)
        delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
        delta_U_matrix = delta_U_matrix_unitless / row / col  # fm^3

        # Get matrix elements for full 3S1-3D1 partial wave channel including
        # the isospin CG = 1/\sqrt(2)
        deltaU = 1/2 * (delta_U_matrix[:self.ntot, :self.ntot]
                        + delta_U_matrix[self.ntot:, self.ntot:])
        deltaU_squared = 1/2 * (delta_U_matrix[:self.ntot, :self.ntot]**2
                                + delta_U_matrix[:self.ntot, self.ntot:]**2
                                + delta_U_matrix[self.ntot:, :self.ntot]**2
                                + delta_U_matrix[self.ntot:, self.ntot:]**2)

        # Interpolate \delta U(k,k)
        self.deltaU_func = interp1d(
            k_array, np.diag(deltaU), kind='linear', bounds_error=False,
            fill_value='extrapolate')

        # Interpolate \delta U^2(k,k') 
        self.deltaU2_func = RectBivariateSpline(k_array, k_array,
                                                deltaU_squared, kx=1, ky=1)

# [sub-classes of MomentumDistributions] Not sure about this yet!
# There should probably be a distinction of multiplying by \theta(...) and
# multiplying by \int dx/2 \theta(...). Where should these functions be?
# 1. Separate .py file in code/modules as sub-class of momentum distributions?

# [sub-classes of MomentumDistributions]
# Do only n_contributions(meshgrids, kF_tau_array, kF_taup_array) which
# calculates the contributions to the momentum distribution (for each script:
# snmd.py, pmd.py, dmd.py).

# [sub-classes of MomentumDistributions]
# Add method that calculates n(q, Q=0) for pmd.py.

# [sub-classes of MomentumDistribution]
# We could do run_momentum_distribution() and save_momentum_distribution()
# where the latter (maybe) relies on np.savetxt().
# Include option for n(q, Q=0).
# Will rely on densities.py. Do not break pn and np convention.
# run_momentum_distribution() will call super().set_matrix_elements() which
# saves \delta U matrix elements (called using super()). Alternatively, can
# make the matrix elements an argument of n_contributions().

# [Main script]
# Add load_momentum_distribution() function which gives the
# momentum distribution given distribution_type = 'pn', 'pp', 'nn', 'p', 'n',
# or 'd', and optional arguments: nucleus, Z, N, density_type,
# interpolate=False, contributions=False, zero_total_momentum=False.
# * This could extend to three functions with an overheaf function if you want.
# * output : Either q_array, (and Q_array?), n_array or n_func.

# [Main script]
# Add method to integrate out Q-dependence of pmd's. See quasideuteron.ipynb
# for example.

# Copy the last function of dmd.py and add to src.ipynb for now. Pretty sure
# it's only used in one figure suggesting that it should be part of the
# function that plots that figure.

# Same idea for partial wave decomposition (see above).

# Add a note to try and understand the K integration in deuteron (does that
# even make sense?)
