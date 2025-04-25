#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: deuteron_wavefunction_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Class for computing the deuteron wave function in JAX.

Last update: April 25, 2025

"""

# Python imports
from functools import partial

# JAX imports
from jax import config, jit
import jax.numpy as jnp
from jax.numpy.linalg import eigh

# Imports from scripts
from .potentials_jax import Potential


# Enable double precision
config.update("jax_enable_x64", True)


class DeuteronWaveFunction(Potential):
    """Class that computes the deuteron wave function."""
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb=jnp.inf):
        
        # Initializes potentials up to L_max = 2
        super().__init__(kvnn, kmax, kmid, ntot, lamb=lamb, L_max=2)
        
        # Set deuteron wave function as instance attributes
        self.psi_0_array, self.psi_2_array = self.set_deuteron_wf(lamb)
        
    def set_deuteron_wf(self, lamb):
        """Get deuteron S- and D-wave functions."""

        # Hamiltonian in 3S1-3D1 coupled-channel [MeV]
        H_matrix = self.load_hamiltonian('3S1', lamb)

        # Diagonalize Hamiltonian for deuteron wave function [unitless]
        eigenvalues, eigenvectors = eigh(H_matrix)
        psi_d_unitless = eigenvectors[:, 0]  # Shape is (2*ntot, 1)
        
        # Remove factor k_i^2 w_i from deuteron wave function [fm^3/2]
        factor_array = self.integration_measure('3S1')
        psi_d_units = psi_d_unitless / factor_array

        # Split into S- and D-waves
        psi_0_array = psi_d_units[:self.ntot]  # S-wave [fm^3/2]
        psi_2_array = psi_d_units[self.ntot:]  # D-wave [fm^3/2]
        
        # Units are [fm^3/2] and shapes are (ntot, 1)
        return psi_0_array, psi_2_array
    
    @partial(jit, static_argnums=(0,))
    def compute_wf(self, k, L):
        """Compute the wave function at relative momenta k given orbital
        angular momentum L.
        """
        
        # List of conditions for S- or D-wave
        condlist = [L == 0, L == 2]
        
        # List of deuteron wave function in S- and D-wave
        choicelist = [
            jnp.interp(k, self.k_array, self.psi_0_array),
            jnp.interp(k, self.k_array, self.psi_2_array)
        ]
        
        # Deuteron wave function in units [fm^3/2] with shape (ntot, 1)
        return jnp.select(condlist, choicelist)