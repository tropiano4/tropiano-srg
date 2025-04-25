#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: delta_U_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

JAX class for evaluating the partial wave relative momentum matrix element of
\delta U where \delta U = I - U, and U is the SRG transformation.

Last update: April 25, 2025

"""

# Python imports
from functools import partial

# JAX imports
from jax import config, jit
from jax.lax import fori_loop
import jax.numpy as jnp
from jax.numpy.linalg import eigh
from jaxinterp2d import interp2d

# Imports from scripts
from .potentials_jax import Potential


# Enable double precision
config.update("jax_enable_x64", True)


class DeltaU(Potential):
    """Class for \delta U in partial waves."""
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max):
        
        # Initializes potentials up to L_max
        super().__init__(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Set SRG flow parameter as instance attribute
        self.lamb = lamb
        
    @partial(jit, static_argnums=(0,))
    def __call__(self, k, kp, J, L, Lp, S, T, hc=False):
        """Compute \delta U in a particular partial wave channel at the momenta
        k and k'. Units are [fm^3]. Option to return hermitian conjugate
        \delta U^\dagger.
        """
        
        delta_U_matrix = self.compute(J, L, Lp, S, T, hc)
        
        # Linear interpolation to the point (k, k')
        return interp2d(k, kp, self.k_array, self.k_array, delta_U_matrix)

    @partial(jit, static_argnums=(0,))
    def compute(self, J, L, Lp, S, T, hc=False):
        """Compute \delta U in a particular partial wave channel. Note, this
        method returns a matrix of zeros if the partial wave channel is 
        unphysical.
        """
        
        cond = self.channel_is_physical(J, L, Lp, S, T)
        
        # \delta U in units [fm^3] with shape (ntot, ntot)
        return jnp.where(
            cond, self.compute_coupled_or_uncoupled(J, L, Lp, S, hc),
            jnp.zeros((self.ntot, self.ntot))
        )
    
    @partial(jit, static_argnums=(0,))
    def compute_coupled_or_uncoupled(self, J, L, Lp, S, hc=False):
        """Compute \delta U in any partial wave channel."""
        
        # Condition on whether the channel is coupled or not
        cond = self.is_coupled(J, L)
        delta_U_matrix = jnp.where(
            cond, self.compute_coupled(J, L, Lp, hc),
            self.compute_uncoupled(J, L, S, hc)
        )
        
        # \delta U in units [fm^3] with shape (ntot, ntot)
        return delta_U_matrix
        
    @partial(jit, static_argnums=(0,))
    def compute_uncoupled(self, J, L, S, hc=False):
        """Compute \delta U in a uncoupled-channel."""
        
        # Get initial and evolved Hamiltonians
        H_initial = self.get_uncoupled_hamiltonian(J, L, S, jnp.inf)
        H_evolved = self.get_uncoupled_hamiltonian(J, L, S, self.lamb)
        
        # Compute SRG transformation U(k, k') which is unitless
        U_matrix_unitless = self.compute_srg_transformation(H_initial,
                                                            H_evolved, hc)
        
        # Subtract out identity matrix
        I_matrix = jnp.identity(self.ntot)
        delta_U_matrix_unitless = U_matrix_unitless - I_matrix
        
        # Convert to units [fm^3]
        delta_U_matrix = self.unattach_weights_uncoupled(
            delta_U_matrix_unitless)

        # \delta U in units [fm^3] with shape (ntot, ntot)
        return delta_U_matrix
        
    @partial(jit, static_argnums=(0,))
    def compute_coupled(self, J, L, Lp, hc=False):
        """Compute \delta U in a coupled-channel."""
        
        # Get initial and evolved Hamiltonians
        H_initial = self.get_coupled_hamiltonian(J, jnp.inf)
        H_evolved = self.get_coupled_hamiltonian(J, self.lamb)
        
        # Compute SRG transformation U(k, k') which is unitless
        U_matrix_unitless = self.compute_srg_transformation(H_initial,
                                                            H_evolved, hc)
        
        # Subtract out identity matrix
        I_matrix = jnp.identity(2 * self.ntot)
        delta_U_matrix_unitless = U_matrix_unitless - I_matrix
        
        # Convert to units [fm^3]
        delta_U_matrix = self.unattach_weights_coupled(delta_U_matrix_unitless)
        
        # Select particular sub-block of coupled-channel matrix
        delta_U_subblock = self.select_subblock(J, L, Lp, delta_U_matrix)

        # \delta U in units [fm^3] with shape (ntot, ntot)
        return delta_U_subblock
    
    @partial(jit, static_argnums=(0,))
    def compute_srg_transformation(self, H_initial, H_evolved, hc=False):
        """SRG unitary transformation built out of eigenvectors of the initial
        and evolved Hamiltonians.
            U = \sum_i | \psi_i(\lambda) > < \psi_i(\inf) |
        """
        
        @jit
        def step(i, U_matrix):
            """Compute outer product of eigenvectors."""
            
            # Individual eigenvectors (sorted correctly from eigh)
            psi_initial = self.vecs_initial[:, i]
            psi_evolved = self.vecs_evolved[:, i]
            
            # Make sure the phases match
            cond = psi_initial.T @ psi_evolved < 0
            psi_evolved = jnp.where(cond, -psi_evolved, psi_evolved)
            
            # Outer product of eigenvectors
            U_matrix += jnp.outer(psi_evolved, psi_initial)
            
            return U_matrix

        Ntot = H_initial.shape[0]

        # Get the eigenvectors of the initial and SRG-evolved Hamiltonians
        _, self.vecs_initial = eigh(H_initial)
        _, self.vecs_evolved = eigh(H_evolved)

        # Initialize unitary transformation U with same size as Hamiltonians
        U_matrix = jnp.zeros((Ntot, Ntot))
        
        # Sum over states i of the Hamiltonian
        U_matrix = fori_loop(0, Ntot, step, U_matrix)
        
        # Take hermitian conjugate?
        return jnp.where(hc, jnp.conj(U_matrix.T), U_matrix)
    
    @partial(jit, static_argnums=(0,))
    def unattach_weights_uncoupled(self, delta_U_matrix_unitless):
        """Divide out integration measure k_i^2 w_i factor from uncoupled
        SRG transformation.
        """
        
        # Get integration factor as a 1-D array [fm^-3/2]
        factor_array = jnp.sqrt(2 / jnp.pi * self.k_weights) * self.k_array
        
        # 1-D array to 2-D grids
        factor_i_grid, factor_j_grid = jnp.meshgrid(factor_array, factor_array,
                                                    indexing='ij')
    
        # Divide the transformation by the integration measure [fm^3]
        return delta_U_matrix_unitless / factor_i_grid / factor_j_grid
    
    @partial(jit, static_argnums=(0,))
    def unattach_weights_coupled(self, delta_U_matrix_unitless):
        """Divide out integration measure k_i^2 w_i factor from coupled SRG
        transformation.
        """
        
        # Get integration factor as a 1-D array [fm^-3/2]
        factor_array = jnp.concatenate(
            (jnp.sqrt(2 / jnp.pi * self.k_weights) * self.k_array,
             jnp.sqrt(2 / jnp.pi * self.k_weights) * self.k_array)
        )
        
        # 1-D array to 2-D grids
        factor_i_grid, factor_j_grid = jnp.meshgrid(factor_array, factor_array,
                                                    indexing='ij')
    
        # Divide the transformation by the integration measure [fm^3]
        return delta_U_matrix_unitless / factor_i_grid / factor_j_grid
    
    @partial(jit, static_argnums=(0,))
    def select_subblock(self, J, L, Lp, delta_U_matrix):
        """Select a particular sub-block (L and L') of a coupled-channel
        \delta U matrix.
        """
        
        # List of conditions for each sub-block
        condlist = [
            jnp.logical_and(L == Lp, J > L),   # 0-0 sub-block
            jnp.logical_and(L != Lp, J > L),   # 0-2 sub-block
            jnp.logical_and(L != Lp, J > Lp),  # 2-0 sub-block
            jnp.logical_and(L == Lp, J < L)    # 2-2 sub-block
        ]
        
        # List of the T-matrix in each sub-block
        choicelist = [
            delta_U_matrix[:self.ntot, :self.ntot],
            delta_U_matrix[:self.ntot, self.ntot:],
            delta_U_matrix[self.ntot:, :self.ntot],
            delta_U_matrix[self.ntot:, self.ntot:]
        ]
        
        # Get particular sub-block with shape (ntot+1, ntot+1)
        delta_U_subblock = jnp.select(condlist, choicelist)
        
        # Units remain [fm]
        return delta_U_subblock

