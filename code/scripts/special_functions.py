#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: special_functions.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 11, 2025

Class for computing special functions (e.g., spherical harmonics) in JAX.

Last update: March 10, 2025

"""

# Python imports
from functools import partial

# JAX imports
import jax.numpy as jnp
from jax import config, jit
from jax.scipy.special import factorial


# Enable double precision
config.update("jax_enable_x64", True)


# TODO: Add spherical Bessel functions of the second kind n_l(\rho)
class SpecialFunctions:
    """Class for computing the first few associated Legendre polynomials and
    spherical harmonics.
    """

    @partial(jit, static_argnums=(0,))
    def ylm(self, theta, phi, L, M_L):
        """Spherical harmonic Y_{L M_L}(\theta, \phi) assuming the quantum
        mechanics sign convention. Note, if |M_L| > L this method returns 0.
        """
        
        sign = (-1) ** M_L
        plm = self.plm(jnp.cos(theta), L, M_L)
        exp = jnp.exp(1j * M_L * phi)
        
        return sign * plm * exp

    @partial(jit, static_argnums=(0,))
    def plm(self, x, L, M_L):
        """Associated Legendre polynomial P_{L M_L}(x). Note, if |M_L| > L
        this method returns 0.
        """
        
        cond = jnp.abs(M_L) <= L
        
        return jnp.where(cond, self.plm_all_m(x, L, M_L), 0)
    
    @partial(jit, static_argnums=(0,))
    def plm_all_m(self, x, L, M_L):
        """Associated Legendre polynomial P_{L M_L}(x) for all M_L."""
        
        cond = M_L < 0
        
        # Include normalization factor in Legendre polynomial
        factor = jnp.sqrt(
            (2 * L + 1) / (4 * jnp.pi) * factorial(L - M_L) / factorial(L + M_L)
        )
        
        return factor * jnp.where(cond, self.plm_negative_m(x, L, M_L),
                                  self.plm_positive_m(x, L, M_L))
        
    @partial(jit, static_argnums=(0,))
    def plm_negative_m(self, x, L, M_L):
        """Convert P_{L M_L}(x) to case where m is negative."""
        
        M_L_positive = -M_L
        sign = (-1) ** (M_L_positive)
        factor = factorial(L - M_L_positive) / factorial(L + M_L_positive)
        plm = self.plm_positive_m(x, L, M_L_positive)
        
        return sign * factor * plm
        
    @partial(jit, static_argnums=(0,))
    def plm_positive_m(self, x, L, M_L):
        """P_{L M_L}(x) for positive m."""
        
        condlist = [
            jnp.logical_and(L == 0, M_L == 0),
            jnp.logical_and(L == 1, M_L == 0),
            jnp.logical_and(L == 1, M_L == 1),
            jnp.logical_and(L == 2, M_L == 0),
            jnp.logical_and(L == 2, M_L == 1),
            jnp.logical_and(L == 2, M_L == 2),
        ]
        choicelist = [
            self.p00(x), self.p10(x), self.p11(x), self.p20(x), self.p21(x),
            self.p22(x)
        ]
        
        return jnp.select(condlist, choicelist)

    @partial(jit, static_argnums=(0,))
    def p00(self, x):
        """Associated Legendre polynomial with L = 0 and M_L = 0."""
        
        return jnp.ones_like(x)
    
    @partial(jit, static_argnums=(0,))
    def p10(self, x):
        """Associated Legendre polynomial with L = 1 and M_L = 0."""
        
        return x
    
    @partial(jit, static_argnums=(0,))
    def p11(self, x):
        """Associated Legendre polynomial with L = 1 and M_L = 1."""
        
        return -jnp.sqrt(1 - x ** 2)
    
    @partial(jit, static_argnums=(0,))
    def p20(self, x):
        """Associated Legendre polynomial with L = 2 and M_L = 0."""
        
        return (3 * x ** 2 - 1) / 2
    
    @partial(jit, static_argnums=(0,))
    def p21(self, x):
        """Associated Legendre polynomial with L = 2 and M_L = 1."""
        
        return -3 * x * jnp.sqrt(1 - x ** 2)
    
    @partial(jit, static_argnums=(0,))
    def p22(self, x):
        """Associated Legendre polynomial with L = 2 and M_L = 2."""
        
        return 3 * (1 - x ** 2)
    
    @partial(jit, static_argnums=(0,))
    def jL(self, rho, L):
        """Spherical Bessel function of the first kind j_L(\rho)."""
        
        condlist = [L == 0, L == 1, L == 2, L == 3]
        choicelist = [self.j0(rho), self.j1(rho), self.j2(rho), self.j3(rho)]
        
        return jnp.select(condlist, choicelist)
    
    @partial(jit, static_argnums=(0,))
    def j0(self, rho):
        """Spherical Bessel function of the first kind with L = 0."""
        
        return jnp.sin(rho) / rho
    
    @partial(jit, static_argnums=(0,))
    def j1(self, rho):
        """Spherical Bessel function of the first kind with L = 1."""
        
        return jnp.sin(rho) / rho ** 2 - jnp.cos(rho) / rho
    
    @partial(jit, static_argnums=(0,))
    def j2(self, rho):
        """Spherical Bessel function of the first kind with L = 2."""
        
        return ((3 / rho ** 2 - 1) * jnp.sin(rho) / rho
                - 3 * jnp.cos(rho) / rho ** 2)
    
    @partial(jit, static_argnums=(0,))
    def j3(self, rho):
        """Spherical Bessel function of the first kind with L = 3."""
        
        return ((15 / rho ** 3 - 6 / rho) * jnp.sin(rho) / rho
                - (15 / rho ** 2 - 1) * jnp.cos(rho) / rho)