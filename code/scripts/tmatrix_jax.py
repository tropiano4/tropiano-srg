#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: tmatrix_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Class for evaluating the NN scattering T-matrix written in JAX.

Last update: April 25, 2025

"""


# Python imports
from functools import partial

# JAX imports
from jax import config, jit
import jax.numpy as jnp
from jax.numpy.linalg import solve
from jaxinterp2d import interp2d

# Imports from scripts
from .potentials_jax import Potential


# Enable double precision
config.update("jax_enable_x64", True)


class TMatrix(Potential):
    """Class that computes the NN T-matrix."""
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max):
        
        # Initializes potentials up to L_max
        super().__init__(kvnn, kmax, kmid, ntot, lamb, L_max)

        # Maximum momentum [fm^-1]
        self.Lamb = kmax
        
    @partial(jit, static_argnums=(0,))
    def half_offshell(self, pp, J, L, Lp, S, T):
        """Compute the half off-shell T-matrix T_{L, L'}(k_i, p'; E'=p'^2/M)
        for a particular partial wave channel.
        """
        
        # T-matrix in units [fm] with shape (ntot+1, ntot+1)
        T_matrix = self.compute(pp, J, L, Lp, S, T)
        
        # Half off-shell T-matrix in units [fm] with shape (ntot, 1)
        return T_matrix[:self.ntot, self.ntot]
    
    @partial(jit, static_argnums=(0,))
    def onshell(self, pp, J, L, Lp, S, T):
        """Compute the on-shell T-matrix T_{L, L'}(p', p'; E'=p'^2/M) for a
        particular partial wave channel.
        """
        
        # T-matrix in units [fm] with shape (ntot+1, ntot+1)
        T_matrix = self.compute(pp, J, L, Lp, S, T)
        
        # Half off-shell T-matrix in units [fm] which is a complex scalar
        return T_matrix[self.ntot, self.ntot]
    
    @partial(jit, static_argnums=(0,))
    def compute(self, pp, J, L, Lp, S, T):
        """Compute the T-matrix for given momentum p' with energy E' = p'^2 / M
        for a particular partial wave channel. Note, this method returns a
        matrix of zeros if the partial wave channel is unphysical.
        """
        
        cond = self.channel_is_physical(J, L, Lp, S, T)
        
        # T-matrix in units [fm] with shape (ntot+1, ntot+1)
        return jnp.where(
            cond,
            self.compute_coupled_or_uncoupled(pp, J, L, Lp, S),
            jnp.zeros((self.ntot + 1, self.ntot + 1))
        )
    
    @partial(jit, static_argnums=(0,))
    def compute_coupled_or_uncoupled(self, pp, J, L, Lp, S):
        """Compute the T-matrix for given momentum p' with energy E' = p'^2 / M
        for any partial wave channel.
        """
        
        # Condition on whether the channel is coupled or not
        cond = self.is_coupled(J, L)
        T_matrix = jnp.where(
            cond, self.compute_coupled(pp, J, L, Lp),
            self.compute_uncoupled(pp, J, L, S)
        )
        
        # T-matrix in units [fm] with shape (ntot+1, ntot+1)
        return T_matrix
    
    @partial(jit, static_argnums=(0,))
    def compute_uncoupled(self, pp, J, L, S):
        """Compute the half off-shell T-matrix in a uncoupled-channel."""

        # Evaluate D-vector for solving matrix inversion problem
        D_vector = self.D_vector(pp)
        
        # Append p' to end of mesh (ntot+1,)
        k_full = jnp.append(self.k_array, pp)
        
        # Create meshes for interpolation (ntot+1, ntot+1)
        k_grid, kp_grid = jnp.meshgrid(k_full, k_full, indexing='ij')
        
        # Load potential in units [fm] with shape (ntot, ntot)
        V_matrix = self.get_uncoupled_potential(J, L, S)

        # Append p' points by linear interpolation (ntot+1, ntot+1)
        V_interpolated = interp2d(k_grid, kp_grid, self.k_array, self.k_array,
                                  V_matrix)
        
        # Build F-matrix [unitless] where F_ij = \delta_ij + D_j V_ij
        F_matrix = (jnp.identity(self.ntot + 1)
                    + jnp.tile(D_vector, (self.ntot + 1, 1)) * V_interpolated)

        # Solve for T-matrix by matrix inversion
        T_matrix = solve(F_matrix, V_interpolated)

        # Units are [fm] and shape is (ntot+1, ntot+1)
        return T_matrix
    
    @partial(jit, static_argnums=(0,))
    def compute_coupled(self, pp, J, L, Lp):
        """Compute the half off-shell T-matrix in a coupled-channel."""
        
        # Evaluate D-vector for solving matrix inversion problem
        D_vector = self.D_vector(pp)
        
        # Append p' to end of mesh (ntot+1,)
        k_full = jnp.append(self.k_array, pp)
        
        # Create meshes for interpolation (ntot+1, ntot+1)
        k_grid, kp_grid = jnp.meshgrid(k_full, k_full, indexing='ij')
        
        # Load potential in units [fm] with shape (2*ntot, 2*ntot)
        V_matrix = self.get_coupled_potential(J)

        # Append p' points by linear interpolation (2*ntot+2, 2*ntot+2)
        V_interpolated = self.interpolate_coupled_potential(k_grid, kp_grid,
                                                            V_matrix)
        
        # Build F-matrix [unitless] where F_ij = \delta_ij + D_j V_ij
        F_matrix = (
            jnp.identity(2 * (self.ntot + 1))
            + jnp.tile(D_vector, (2 * (self.ntot + 1), 2)) * V_interpolated
        )

        # Solve for T-matrix by matrix inversion
        T_matrix = solve(F_matrix, V_interpolated)
        
        # Select particular sub-block of coupled-channel matrix
        T_subblock = self.select_subblock(J, L, Lp, T_matrix)

        # Units are [fm] and shape is (ntot+1, ntot+1)
        return T_subblock
    
    @partial(jit, static_argnums=(0,))
    def D_vector(self, pp):
        """Compute D-vector used for solving matrix inversion problem:
            T = F^-1 V
        where F_ij = \delta_ij + D_j V_ij.
        """
        
        # First ntot elements of D_vector [fm^-1]
        D_vector = (2.0 / jnp.pi * (self.k_weights * self.k_array ** 2)
                    / (self.k_array ** 2 - pp ** 2))
        
        # ntot + 1 element of D_vector [fm^-1]
        D_last = (
            -2.0 / jnp.pi * pp ** 2 * (
                jnp.sum(self.k_weights / (self.k_array ** 2 - pp ** 2))
                + jnp.log((self.Lamb + pp) / (self.Lamb - pp)) / (2.0 * pp)
            )
        ) + 1j * pp
        
        # Append ntot + 1 element to D_vector
        return jnp.append(D_vector, D_last)

    @partial(jit, static_argnums=(0,))
    def interpolate_coupled_potential(self, k_grid, kp_grid, V_matrix):
        """Interpolate half off-shell T-matrix matrix in a coupled-channel."""
        
        # Get each sub-block separately with shapes (ntot, ntot)
        V11, V12, V21, V22 = self.get_subblocks(V_matrix)
        
        # Append p' points by linear interpolation (ntot+1, ntot+1)
        V11_interpolated = interp2d(k_grid, kp_grid, self.k_array, self.k_array,
                                    V11)
        V12_interpolated = interp2d(k_grid, kp_grid, self.k_array, self.k_array,
                                    V12)
        V21_interpolated = interp2d(k_grid, kp_grid, self.k_array, self.k_array,
                                    V21)
        V22_interpolated = interp2d(k_grid, kp_grid, self.k_array, self.k_array,
                                    V22)
        
        # Recombine sub-blocks for coupled-channel matrix
        V_interpolated = jnp.vstack((
            jnp.hstack((V11_interpolated, V12_interpolated)),
            jnp.hstack((V21_interpolated, V22_interpolated))
        ))

        # Shape is (2*ntot+2, 2*ntot+2)
        return V_interpolated
    
    @partial(jit, static_argnums=(0,))
    def select_subblock(self, J, L, Lp, T_matrix):
        """Select a particular sub-block (L and L') of a coupled-channel T-
        matrix.
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
            T_matrix[:self.ntot+1, :self.ntot+1],
            T_matrix[:self.ntot+1, self.ntot+1:],
            T_matrix[self.ntot+1:, :self.ntot+1],
            T_matrix[self.ntot+1:, self.ntot+1:]
        ]
        
        # Get particular sub-block with shape (ntot+1, ntot+1)
        T_subblock = jnp.select(condlist, choicelist)
        
        # Units remain [fm]
        return T_subblock
    
    @partial(jit, static_argnums=(0,))
    def get_subblocks(self, V_matrix):
        """Gets each sub-block (L and L') from a coupled-channel potential."""
        
        V11 = V_matrix[:self.ntot, :self.ntot]
        V12 = V_matrix[:self.ntot, self.ntot:]
        V21 = V_matrix[self.ntot:, :self.ntot]
        V22 = V_matrix[self.ntot:, self.ntot:]

        return V11, V12, V21, V22