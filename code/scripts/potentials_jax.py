#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: potentials_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 12, 2025

Class for handling NN potentials using JAX arrays.

Last update: March 12, 2025

"""

# Python imports
from functools import partial
import numpy as np

# JAX imports
import jax.numpy as jnp
from jax import config, jit

# Imports from scripts
from .tools import channel_L_value, coupled_channel


# Enable double precision
config.update("jax_enable_x64", True)


class Potential:
    """Class that loads NN potentials as JAX arrays."""
    
    # Define class attribute for h-bar^2 / M [MeV fm^2]
    hbar_sq_over_m = 41.47
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb=jnp.inf, L_max=2):
        
        # Need kvnn as string (can cause error if kvnn < 10)
        if kvnn < 10:
            kvnn_str = '0' + str(kvnn)
        else:
            kvnn_str = str(kvnn)
        self.kvnn_str = kvnn_str
        
        # Get potential directory
        kmax_int = int(kmax)
        kmid_int = int(kmid)
        self.directory = (
            f'../data/potentials/vsrg_kvnn_{kvnn_str}_lam12.0_kmax{kmax_int:d}'
            f'_kmid{kmid_int:d}_ntot{ntot:d}/'
        )

        # Set momentum mesh in units fm^-1
        self.k_array, self.k_weights = self.get_momentum_mesh()
        
        # Set momentum mesh specifications as attributes
        self.kmax, self.kmid, self.ntot = kmax, kmid, ntot
        
        # Set potentials up to L_max in a big JAX arrays distinguishing
        # uncoupled- or coupled-channel
        self.uncoupled_potentials, self.coupled_potentials = (
            self.get_potentials(lamb, L_max)
        )
        
        # Set unevolved and evolved Hamiltonians in big JAX arrays
        # distinguishing uncoupled- or coupled-channel
        if lamb != jnp.inf:
            
            # Un-evolved Hamiltonians
            (
                self.uncoupled_hamiltonians_unevolved, 
                self.coupled_hamiltonians_unevolved
            ) = self.get_hamiltonians(jnp.inf, L_max)
            
            # Evolved Hamiltonians
            (
                self.uncoupled_hamiltonians_evolved, 
                self.coupled_hamiltonians_evolved
            ) = self.get_hamiltonians(lamb, L_max)
        
    def get_momentum_mesh(self):
        """Momentum mesh in units [fm^-1] as JAX arrays."""
        
        filename = f'vsrg_1S0_kvnn_{self.kvnn_str}_lam12.0_reg_0_3_0_mesh.out'
        momentum_mesh = np.loadtxt(self.directory + filename)
        k_array = jnp.asarray(momentum_mesh[:, 0])
        k_weights = jnp.asarray(momentum_mesh[:, 1])
        
        return k_array, k_weights
    
    def get_potentials(self, lamb, L_max):
        """JAX arrays of potentials up to L_max returning uncoupled- and
        coupled-channel separately."""
        
        # Initialize potentials as lists later to be converted to JAX arrays
        uncoupled_potentials_list = []
        coupled_potentials_list = []
        
        # Possible partial wave channels up to L = 5
        channels = [
            '1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3',
            '1F3', '3F3', '3F4', '1G4', '3G4', '3G5', '1H5', '3H5', '3H6'
        ]
        
        # Loop over channels
        for channel in channels:

            # Check that channel is within L_max
            L = channel_L_value(channel)
            if L <= L_max:
                
                # Load potential [fm]
                V_matrix = self.load_potential(channel, lamb)

                # Coupled-channel
                if coupled_channel(channel):
                    coupled_potentials_list.append(V_matrix)
                # Uncoupled-channel
                else:
                    uncoupled_potentials_list.append(V_matrix)
                    
        # Convert lists to JAX arrays with shape (# of channels, ntot, ntot)
        uncoupled_potentials = jnp.asarray(uncoupled_potentials_list)
        coupled_potentials = jnp.asarray(coupled_potentials_list)

        return uncoupled_potentials, coupled_potentials
    
    def get_hamiltonians(self, lamb, L_max):
        """JAX arrays of Hamiltonians up to L_max returning uncoupled- and
        coupled-channel separately."""
        
        # Initialize Hamiltonians as lists later to be converted to JAX arrays
        uncoupled_hamiltonians_list = []
        coupled_hamiltonians_list = []
        
        # Possible partial wave channels up to L = 5
        channels = [
            '1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3',
            '1F3', '3F3', '3F4', '1G4', '3G4', '3G5', '1H5', '3H5', '3H6'
        ]
        
        # Loop over channels
        for channel in channels:

            # Check that channel is within L_max
            L = channel_L_value(channel)
            if L <= L_max:
                
                # Load potential [MeV]
                H_matrix = self.load_hamiltonian(channel, lamb)

                # Coupled-channel
                if coupled_channel(channel):
                    coupled_hamiltonians_list.append(H_matrix)
                # Uncoupled-channel
                else:
                    uncoupled_hamiltonians_list.append(H_matrix)
                    
        # Convert lists to JAX arrays with shape (# of channels, ntot, ntot)
        uncoupled_hamiltonians = jnp.asarray(uncoupled_hamiltonians_list)
        coupled_hamiltonians = jnp.asarray(coupled_hamiltonians_list)

        return uncoupled_hamiltonians, coupled_hamiltonians
            
    def load_potential(self, channel, lamb):
        """Load the potential in the given partial wave channel."""
        
        # Unevolved potential
        if lamb == jnp.inf:

            filename = (f'vnn_{channel}_kvnn_{self.kvnn_str}_lam12.0_reg_0_3_0'
                        '.out')
        
        # Evolved potential
        else:
            
            # Get \lambda with correct number of decimals
            if lamb == round(lamb, 1):
                lamb_str = str(round(lamb, 1))
            else:
                lamb_str = str(round(lamb, 2))
            
            filename = (f'vnn_{channel}_kvnn_{self.kvnn_str}_srg_Wegner'
                        f'_lambda{lamb_str}.out')

        data = np.loadtxt(self.directory + filename)
        
        # Coupled-channel potential?
        if coupled_channel(channel):
        
            V11 = jnp.reshape(data[:, 2], (self.ntot, self.ntot))
            V12 = jnp.reshape(data[:, 3], (self.ntot, self.ntot))
            V21 = jnp.reshape(data[:, 4], (self.ntot, self.ntot))
            V22 = jnp.reshape(data[:, 5], (self.ntot, self.ntot))
            V_matrix = jnp.vstack(
                (jnp.hstack((V11, V12)), jnp.hstack((V21, V22)))
            )
        
        else:
        
            V_matrix = jnp.reshape(data[:, 2], (self.ntot, self.ntot))

        # Potential in units [fm] where the shape can vary
        return V_matrix
    
    def load_kinetic_energy(self, channel):
        """Loads relative kinetic energy."""

        # Matrix of (h-bar*k)^2 / M along diagonal [MeV]
        T_matrix = jnp.diag(self.k_array ** 2) * self.hbar_sq_over_m
    
        # Coupled-channel operator?
        if coupled_channel(channel):
        
            # Matrix of zeros (n x n)
            zeros = jnp.zeros((self.ntot, self.ntot))
        
            # Build coupled-channel kinetic energy matrix
            T_matrix = jnp.vstack(
                (jnp.hstack((T_matrix, zeros)), jnp.hstack((zeros, T_matrix)))
            )
        
        # Kinetic energy in units [MeV]
        return T_matrix
    
    def load_hamiltonian(self, channel, lamb):
        """Loads the Hamiltonian."""
        
        # Load relative kinetic energy [MeV]
        T_matrix = self.load_kinetic_energy(channel)
    
        # Load potential [fm]
        V_matrix_fm = self.load_potential(channel, lamb)

        # Convert potential from [fm] -> [MeV]
        V_matrix_MeV = self.hbar_sq_over_m * self.attach_integration_weights(
            V_matrix_fm, channel
        )
    
        H_matrix = T_matrix + V_matrix_MeV
    
        # Hamiltonian in units MeV
        return H_matrix
    
    def attach_integration_weights(self, V_matrix_fm, channel):
        """Attach integration k_i^2 w_i factor to potential."""
        
        # Get integration factor as a 1-D array
        factor_array = self.integration_measure(channel)
        
        # 1-D array to 2-D grids
        factor_i_grid, factor_j_grid = jnp.meshgrid(factor_array, factor_array,
                                                    indexing='ij')
    
        # Multiply the potential by the integration measure [fm^-2]
        return factor_i_grid * factor_j_grid * V_matrix_fm
    
    def integration_measure(self, channel):
        """Get factor k_i^2 w_i."""
        
        # Make array of integration factor [fm^-3/2]
        if coupled_channel(channel):
            factor_array = jnp.concatenate((
                jnp.sqrt(2 / jnp.pi * self.k_weights) * self.k_array,
                jnp.sqrt(2 / jnp.pi * self.k_weights) * self.k_array
            ))
        else:
            factor_array = jnp.sqrt(2 / jnp.pi * self.k_weights) * self.k_array
        
        # Return 1-D array with units [fm^-3/2]
        return factor_array

    @partial(jit, static_argnums=(0,))
    def get_uncoupled_potential(self, L, S, J):
        
        # Map L, S, and J onto index (special case for 3P0)
        cond = jnp.logical_and(L == 1, jnp.logical_and(S == 1, J == 0))
        index = jnp.where(cond, 1, L + S + J)
        
        # Units are [fm] and shape is (ntot, ntot)
        return self.uncoupled_potentials[index]
    
    @partial(jit, static_argnums=(0,))
    def get_coupled_potential(self, J):
        
        # Map J onto index
        index = J - 1
        
        # Units are [fm] and shape is (2*ntot, 2*ntot)
        return self.coupled_potentials[index]
    
    @partial(jit, static_argnums=(0,))
    def get_uncoupled_hamiltonian(self, L, S, J, lamb):
        
        # Map L, S, and J onto index (special case for 3P0)
        cond = jnp.logical_and(L == 1, jnp.logical_and(S == 1, J == 0))
        index = jnp.where(cond, 1, L + S + J)
        
        # Pick unevolved or evolved Hamiltonian
        cond = lamb == jnp.inf
        H_matrix = jnp.where(cond, self.uncoupled_hamiltonians_unevolved[index],
                             self.uncoupled_hamiltonians_evolved[index])
        
        # Units are [MeV] and shape is (ntot, ntot)
        return H_matrix
    
    @partial(jit, static_argnums=(0,))
    def get_coupled_hamiltonian(self, J, lamb):
        
        # Map J onto index
        index = J - 1
        
        # Pick unevolved or evolved Hamiltonian
        cond = lamb == jnp.inf
        H_matrix = jnp.where(cond, self.coupled_hamiltonians_unevolved[index],
                             self.coupled_hamiltonians_evolved[index])
        
        # Units are [fm] and shape is (2*ntot, 2*ntot)
        return H_matrix
    
    @partial(jit, static_argnums=(0,))
    def is_coupled(self, L, J):
        """Boolean value on whether the channel is coupled or not."""
        
        return jnp.logical_and(J > 0, J != L)
    
    @partial(jit, static_argnums=(0,))
    def channel_is_physical(self, L, Lp, S, J, T):
        """Check if the partial wave channel is physical."""

        # Make sure |L - S| <= J <= L + S and likewise with L'
        J_bool = jnp.logical_and(
            jnp.logical_and(
                jnp.less_equal(jnp.abs(L - S), J),
                jnp.less_equal(J, L + S)
            ),
            jnp.logical_and(
                jnp.less_equal(jnp.abs(Lp - S), J),
                jnp.less_equal(J, Lp + S)
            )
        )
        
        # Make sure L + S + T and L' + S + T are odd
        T_bool = jnp.logical_and(
            jnp.not_equal((L + S + T) % 2, 0),
            jnp.not_equal((Lp + S + T) % 2, 0)
        )
    
        return jnp.logical_and(J_bool, T_bool)