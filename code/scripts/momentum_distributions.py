#!/usr/bin/env python3

"""
File: momentum_distributions.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 17, 2022

Initializes momentum distribution codes. Momentum distributions are computed
using subclasses from snmd.py, pmd.py, and dmd.py, and SRG-transformed
potentials from potentials.py. Momentum distribution data is stored in
data/momentum_distributions. These codes also rely on nucleonic densities
taken from external codes or data.

Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!

Last update: July 28, 2022

"""

# To-do: Further reduce load methods into smaller methods.
# To-do: Make save functions automatically create /nucleus_name subdirectories

# Python imports
import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline

# Imports from A.T. codes
from .integration import (
    unattach_weights_from_matrix, unattach_weights_from_vector
)
from .potentials import Potential
from .srg import get_transformation
from .tools import channel_L_value, coupled_channel, replace_periods
from .wave_function import wave_function


class MomentumDistribution:
    """
    Parent class of the SingleNucleon, Pair, and Deuteron momentum distribution
    classes. This sets up the calculation using SRG-evolved potentials from
    the Potentials class.
    
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

    def __init__(self, kvnn, kmax, kmid, ntot):
        """Save the inputs of the potential excluding the channels argument."""

        self.kvnn = kvnn
        self.kmax = kmax
        self.kmid = kmid
        self.ntot = ntot

    def get_hamiltonians(
            self, channel, generator, lamb, lambda_initial=None, kvnn_inv=None,
            lambda_m=None):
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
        lambda_m : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and lambda_m must be specified for
            this to run.
            
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
        if lambda_initial is None:
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
        if kvnn_inv is not None:
            potential_hard = Potential(kvnn_inv, channel, self.kmax, self.kmid,
                                       self.ntot)
            H_hard_initial = potential_hard.load_hamiltonian()
            if generator == 'Block-diag':
                H_hard_evolved = potential_hard.load_hamiltonian(
                    'srg', generator, 1.0, lambda_bd=lambda_m)
            else:
                H_hard_evolved = potential_hard.load_hamiltonian(
                    'srg', generator, lambda_m)
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
            deltaU_squared = (delta_U_matrix[:self.ntot, :self.ntot] ** 2
                              + delta_U_matrix[:self.ntot, self.ntot:] ** 2)

            # Decide whether to add second L based on highest allowed L value
            # (e.g., include the 3D1-3D1 part of the coupled 3S1-3D1 channel
            # if we input D-waves in channels)
            if channel_L_value(channel) + 2 <= self.highest_L:
                deltaU += delta_U_matrix[self.ntot:, self.ntot:]
                deltaU_squared += (
                    delta_U_matrix[self.ntot:, :self.ntot] ** 2
                    + delta_U_matrix[self.ntot:, self.ntot:] ** 2
                )

        else:

            deltaU = delta_U_matrix
            deltaU_squared = delta_U_matrix ** 2

        return deltaU, deltaU_squared

    def save_deltaU_funcs(
            self, channels, generator, lamb, lambda_initial=None,
            kvnn_inv=None, lambda_m=None):
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
        lambda_m : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and lambda_m must be specified for
            this to run.

        Notes
        -----
        We are taking \lambda=1.0 fm^-1 for block-diagonal decoupling and
        assuming parameters lamb, lambda_initial, and lambda_m correspond
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
                channel, generator, lamb, lambda_initial, kvnn_inv, lambda_m)

            # Get SRG transformation U(k, k') [unitless]
            U_matrix_unitless = get_transformation(H_initial, H_evolved)

            # Coupled-channel?
            cc_bool = coupled_channel(channel)

            # Isolate 2-body term and convert to fm^3
            if cc_bool:
                I_matrix_unitless = np.eye(2 * self.ntot, 2 * self.ntot)
            else:
                I_matrix_unitless = np.eye(self.ntot, self.ntot)
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            delta_U_matrix = unattach_weights_from_matrix(
                k_array, k_weights, delta_U_matrix_unitless, cc_bool)  # fm^3

            # 2J+1 factor
            J = int(channel[-1])

            # Get matrix elements up to the highest input partial wave channel
            deltaU, deltaU2 = self.get_deltaU_matrix_element(channel,
                                                             delta_U_matrix)

            # Isospin CG's=1/\sqrt(2) for pn
            deltaU_pn += (2 * J + 1) / 2 * deltaU
            deltaU2_pn += (2 * J + 1) / 2 * deltaU2

            # Isospin CG's=1 for pp
            if channel in pp_channels:
                deltaU_pp += (2 * J + 1) * deltaU
                deltaU2_pp += (2 * J + 1) * deltaU2

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
            self, generator, lamb, lambda_initial=None, kvnn_inv=None,
            lambda_m=None):
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
        lambda_m : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and lambda_m must be specified for
            this to run.

        Notes
        -----
        We are taking \lambda=1.0 fm^-1 for block-diagonal decoupling and
        assuming parameters lamb, lambda_initial, and lambda_m correspond
        to \Lambda_BD.
        
        """

        # Channel is 3S1-3D1 for deuteron
        channel = '3S1'

        # Get momentum mesh (channel argument doesn't matter here)
        k_array, k_weights = Potential(
            self.kvnn, channel, self.kmax, self.kmid, self.ntot).load_mesh()

        # Set mesh as instance attribute
        self.k_array, self.k_weights = k_array, k_weights

        # Get initial and evolved Hamiltonians
        H_initial, H_evolved = self.get_hamiltonians(
            channel, generator, lamb, lambda_initial, kvnn_inv, lambda_m)

        # Get SRG transformation U(k, k') [unitless]
        U_matrix_unitless = get_transformation(H_initial, H_evolved)

        # Need the deuteron wave function to get kF values for calculation
        psi_k_unitless = wave_function(H_initial, U_matrix=U_matrix_unitless)
        self.psi_k = unattach_weights_from_vector(
            k_array, k_weights, psi_k_unitless, coupled_channel=True)  # fm^3/2

        # Isolate 2-body term and convert to fm^3
        I_matrix_unitless = np.eye(2 * self.ntot, 2 * self.ntot)
        delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
        delta_U_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delta_U_matrix_unitless, coupled_channel=True)

        # Get matrix elements for full 3S1-3D1 partial wave channel including
        # the isospin CG = 1/\sqrt(2)
        deltaU = 1 / 2 * (delta_U_matrix[:self.ntot, :self.ntot]
                          + delta_U_matrix[self.ntot:, self.ntot:])
        deltaU_squared = 1 / 2 * (delta_U_matrix[:self.ntot, :self.ntot] ** 2
                                  + delta_U_matrix[:self.ntot, self.ntot:] ** 2
                                  + delta_U_matrix[self.ntot:, :self.ntot] ** 2
                                  + delta_U_matrix[self.ntot:, self.ntot:] ** 2)

        # Interpolate \delta U(k,k)
        self.deltaU_func = interp1d(
            k_array, np.diag(deltaU), kind='linear', bounds_error=False,
            fill_value='extrapolate')

        # Interpolate \delta U^2(k,k') 
        self.deltaU2_func = RectBivariateSpline(k_array, k_array,
                                                deltaU_squared, kx=1, ky=1)

    def get_single_nucleon_momentum_distribution(
            self, nucleon, nucleus_name, density, channels, generator, lamb,
            lambda_initial=None, kvnn_inv=None, lambda_m=None,
            contributions=False, interpolate=False):
        """
        Loads single-nucleon momentum distribution data from
        data/momentum_distributions. Here the arguments of the function
        specify the file name and directory. There is an option to return
        interpolated function(s) instead.

        Parameters
        ----------
        nucleon : str
            Specify 'proton' or 'neutron'.
        nucleus_name : str
            Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
        density : str
            Name of density (e.g., 'SLy4', 'Gogny').
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
        lambda_m : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and lambda_m must be specified
            for this to run.
        contributions : bool, optional
            Option to return isolated contributions to the momentum
            distribution in terms of the I, \delta U, and \delta U^2 terms.
        interpolate : bool, optional
            Option to return interpolated function instead of NumPy array.
            
        Returns
        -------
        output : tuple or func
            The output depends on the arguments of contributions and
            interpolate. The momentum distribution has units [fm^3] and the
            momentum values has units [fm^-1].

        """

        # Directory for distributions data
        data_directory = f'../data/momentum_distributions/{nucleus_name}/'

        # Get file name
        file_name = (f'n_{nucleon}_{density}_kvnn_{self.kvnn}_kmax_{self.kmax}'
                     f'_kmid_{self.kmid}_ntot_{self.ntot}')

        for channel in channels:
            file_name += f'_{channel}'

        file_name += f'_{generator}'
        if generator == 'Block-diag':
            file_name += f'_LambdaBD_{lamb}'
        else:
            file_name += f'_lambda_{lamb}'

        if lambda_initial is not None:
            file_name += f'_lambda_initial_{lambda_initial}'

        if kvnn_inv is not None:
            file_name += f'_kvnn_inv_{kvnn_inv}_lamb_m_{lambda_m}'

        file_name = replace_periods(file_name) + '.dat'

        # Load data which includes all contributions to n_\lambda(q)
        data = np.loadtxt(data_directory + file_name)

        # Momentum values and total distribution (1-D arrays)
        q_array = data[:, 0]
        n_total_array = data[:, 1]

        # Get isolated contributions too (1-D arrays)
        if contributions:
            n_I_array = data[:, 2]  # 1 term
            n_delU_array = data[:, 3]  # \delta U term
            n_delU2_array = data[:, 4]  # \delta U^2 term

        # Interpolate (UnivariateSpline is for smoothing whereas interp1d
        # gives closer value to the actual calculation)
        if interpolate:

            # Total distribution (function)
            n_total_func = interp1d(q_array, n_total_array, bounds_error=False,
                                    kind='linear', fill_value='extrapolate')

            # Isolated contributions too (functions)
            if contributions:
                n_I_func = interp1d(q_array, n_I_array, bounds_error=False,
                                    kind='linear', fill_value='extrapolate')
                n_delU_func = interp1d(
                    q_array, n_delU_array, bounds_error=False, kind='linear',
                    fill_value='extrapolate')
                n_delU2_func = interp1d(
                    q_array, n_delU2_array, bounds_error=False, kind='linear',
                    fill_value='extrapolate')

        # Return all contributions as functions
        if contributions and interpolate:
            return n_total_func, n_I_func, n_delU_func, n_delU2_func
        # Return all contributions as arrays of the data
        elif contributions:
            return (q_array, n_total_array, n_I_array, n_delU_array,
                    n_delU2_array)
        # Return the total momentum distribution as a function
        elif interpolate:
            return n_total_func
        # Return the total momentum distribution as an array of the data
        else:
            return q_array, n_total_array

    def get_pair_momentum_distribution(
            self, pair, nucleus_name, density, channels, generator, lamb,
            lambda_initial=None, kvnn_inv=None, lambda_m=None,
            Q_equals_zero=False, contributions=False, interpolate=False):
        """
        Loads pair momentum distribution data from data/momentum_distributions.
        Here the arguments of the function specify the file name and directory.
        There is an option to return interpolated function(s) instead.

        Parameters
        ----------
        pair : str
            Type of pair momentum distribution ('pp', 'pn', or 'nn').
        nucleus_name : str
            Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
        density : str
            Name of density (e.g., 'SLy4', 'Gogny').
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
        lambda_m : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and lambda_m must be specified
            for this to run.
        Q_equals_zero : bool, optional
            Option to get the pair momentum distribution at Q = 0: n(q, Q=0)
            which has units [fm^6].
        contributions : bool, optional
            Option to return isolated contributions to the momentum
            distribution in terms of the I, \delta U, and \delta U^2 terms.
        interpolate : bool, optional
            Option to return interpolated function instead of NumPy array.
            
        Returns
        -------
        output : tuple or func
            The output depends on the arguments of contributions and
            interpolate. The momentum distribution has units [fm^6] and the
            momentum values has units [fm^-1].

        """

        # Directory for distributions data
        data_directory = f'../data/momentum_distributions/{nucleus_name}/'

        # Get file name
        file_name = (f'n_{pair}_{density}_kvnn_{self.kvnn}_kmax_{self.kmax}'
                     f'_kmid_{self.kmid}_ntot_{self.ntot}')

        for channel in channels:
            file_name += f'_{channel}'

        file_name += f'_{generator}'
        if generator == 'Block-diag':
            file_name += f'_LambdaBD_{lamb}'
        else:
            file_name += f'_lambda_{lamb}'

        if lambda_initial is not None:
            file_name += f'_lambda_initial_{lambda_initial}'

        if kvnn_inv is not None:
            file_name += f'_kvnn_inv_{kvnn_inv}_lamb_m_{lambda_m}'

        # Split into cases on whether Q = 0
        if Q_equals_zero:

            file_name = replace_periods(file_name + '_Q0') + '.dat'

            # Load data which includes all contributions to n_\lambda(q,Q=0)
            data = np.loadtxt(data_directory + file_name)

            # Momentum values and total distribution (1-D arrays)
            q_array = data[:, 0]
            n_total_array = data[:, 1]

            # Get isolated contributions too (1-D arrays)
            if contributions:
                n_I_array = data[:, 2]  # 1 term
                n_delU_array = data[:, 3]  # \delta U term
                n_delU2_array = data[:, 4]  # \delta U^2 term

            # Interpolate (UnivariateSpline is for smoothing whereas interp1d
            # gives closer value to the actual calculation)
            if interpolate:

                # Total distribution (function)
                n_total_func = interp1d(
                    q_array, n_total_array, bounds_error=False, kind='linear',
                    fill_value='extrapolate')

                # Isolated contributions too (functions)
                if contributions:
                    n_I_func = interp1d(
                        q_array, n_I_array, bounds_error=False, kind='linear',
                        fill_value='extrapolate')
                    n_delU_func = interp1d(
                        q_array, n_delU_array, bounds_error=False,
                        kind='linear', fill_value='extrapolate')
                    n_delU2_func = interp1d(
                        q_array, n_delU2_array, bounds_error=False,
                        kind='linear', fill_value='extrapolate')

            # Return all contributions as functions
            if contributions and interpolate:
                return n_total_func, n_I_func, n_delU_func, n_delU2_func
            # Return all contributions as arrays of the data
            elif contributions:
                return (q_array, n_total_array, n_I_array, n_delU_array,
                        n_delU2_array)
            # Return the total momentum distribution as a function
            elif interpolate:
                return n_total_func
            # Return the total momentum distribution as an array of the data
            else:
                return q_array, n_total_array

        # Q > 0
        else:

            file_name = replace_periods(file_name) + '.dat'

            # Load data which includes all contributions to n_\lambda(q,Q)
            data = np.loadtxt(data_directory + file_name)

            # Get C.o.M. and relative momentum values
            q_array = np.unique(data[:, 0])
            Q_array = np.unique(data[:, 1])
            ntot_q, ntot_Q = len(q_array), len(Q_array)

            # Total distribution (2-D array)
            n_total_array = np.reshape(data[:, 2], (ntot_q, ntot_Q))

            # Get isolated contributions too (2-D arrays)
            if contributions:
                # I term
                n_I_array = np.reshape(data[:, 3], (ntot_q, ntot_Q))
                # \delta U term
                n_delU_array = np.reshape(data[:, 4], (ntot_q, ntot_Q))
                # \delta U^2 term
                n_delU2_array = np.reshape(data[:, 5], (ntot_q, ntot_Q))

            if interpolate:
                n_total_func = RectBivariateSpline(q_array, Q_array,
                                                   n_total_array, kx=1, ky=1)

                # Isolated contributions too (functions)
                if contributions:
                    n_I_func = RectBivariateSpline(q_array, Q_array,
                                                   n_I_array, kx=1, ky=1)
                    n_delU_func = RectBivariateSpline(q_array, Q_array,
                                                      n_delU_array, kx=1, ky=1)
                    n_delU2_func = RectBivariateSpline(
                        q_array, Q_array, n_delU2_array, kx=1, ky=1)

            # Return all contributions as functions
            if contributions and interpolate:
                return n_total_func, n_I_func, n_delU_func, n_delU2_func
            # Return all contributions as arrays of the data
            elif contributions:
                return (q_array, Q_array, n_total_array, n_I_array,
                        n_delU_array, n_delU2_array)
            # Return the total momentum distribution as a function
            elif interpolate:
                return n_total_func
            # Return the total momentum distribution as an array of the data
            else:
                return q_array, Q_array, n_total_array

    def get_deuteron_momentum_distribution(
            self, generator, lamb, lambda_initial=None, kvnn_inv=None,
            lambda_m=None, contributions=False, interpolate=False):
        """
        Loads deuteron momentum distribution data from
        data/momentum_distributions. Here the arguments of the function
        specify the file name and directory. There is an option to return
        interpolated function(s) instead.

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
        lambda_m : float, optional
            SRG evolution parameter \lambda for inverse-SRG transformations
            [fm^-1]. Note, both kvnn_inv and lambda_m must be specified for
            this to run.
        contributions : bool, optional
            Option to return isolated contributions to the momentum
            distribution in terms of the I, \delta U, and \delta U^2 terms.
        interpolate : bool, optional
            Option to return interpolated function instead of NumPy array.
            
        Returns
        -------
        output : tuple or func
            The output depends on the arguments of contributions and
            interpolate. The momentum distribution has units [fm^3] and the
            momentum values has units [fm^-1].

        """

        # Directory for distributions data
        data_directory = '../data/momentum_distributions/H2/'

        # Get file name
        file_name = (f'n_kvnn_{self.kvnn}_kmax_{self.kmax}_kmid_{self.kmid}'
                     f'_ntot_{self.ntot}_{generator}')

        if generator == 'Block-diag':
            file_name += f'_LambdaBD_{lamb}'
        else:
            file_name += f'_lambda_{lamb}'

        if lambda_initial is not None:
            file_name += f'_lambda_initial_{lambda_initial}'

        if kvnn_inv is not None:
            file_name += f'_kvnn_inv_{kvnn_inv}_lamb_m_{lambda_m}'

        file_name = replace_periods(file_name) + '.dat'

        # Load data which includes all contributions to n_\lambda(q)
        data = np.loadtxt(data_directory + file_name)

        # Momentum values and total distribution (1-D arrays)
        q_array = data[:, 0]
        n_total_array = data[:, 1]

        # Get isolated contributions too (1-D arrays)
        if contributions:
            n_I_array = data[:, 2]  # 1 term
            n_delU_array = data[:, 3]  # \delta U term
            n_delU2_array = data[:, 4]  # \delta U^2 term

        # Interpolate (UnivariateSpline is for smoothing whereas interp1d
        # gives closer value to the actual calculation)
        if interpolate:

            # Total distribution (function)
            n_total_func = interp1d(q_array, n_total_array, bounds_error=False,
                                    kind='linear', fill_value='extrapolate')

            # Isolated contributions too (functions)
            if contributions:
                n_I_func = interp1d(q_array, n_I_array, bounds_error=False,
                                    kind='linear', fill_value='extrapolate')
                n_delU_func = interp1d(
                    q_array, n_delU_array, bounds_error=False, kind='linear',
                    fill_value='extrapolate')
                n_delU2_func = interp1d(
                    q_array, n_delU2_array, bounds_error=False, kind='linear',
                    fill_value='extrapolate')

        # Return all contributions as functions
        if contributions and interpolate:
            return n_total_func, n_I_func, n_delU_func, n_delU2_func
        # Return all contributions as arrays of the data
        elif contributions:
            return (q_array, n_total_array, n_I_array, n_delU_array,
                    n_delU2_array)
        # Return the total momentum distribution as a function
        elif interpolate:
            return n_total_func
        # Return the total momentum distribution as an array of the data
        else:
            return q_array, n_total_array
