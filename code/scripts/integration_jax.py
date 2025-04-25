#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: integration_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Class for integrating with Gaussian quadrature written in JAX.

Last update: April 25, 2025

"""

# Python imports
from functools import partial

# JAX imports
from jax import config, jit
from jax.lax import fori_loop, while_loop
import jax.numpy as jnp


# Enable double precision
config.update("jax_enable_x64", True)


class GaussQuadrature:
    """Gauss-Legendre quadrature integration in JAX."""
    
    def __init__(self, n):
        """Save number of points as an instance attribute. Note, make sure n is
        an odd valued integer.
        """
        
        # Check that n is odd
        if n % 2 == 0:
            raise ValueError("n must be an odd integer!")
        
        # Number of points within (-1, 1)
        self.n = n
        
        # Number of points within (-1, 0)
        self.m = int((n + 1) / 2)
        
        # Tolerance in accepting root
        self.eps = 3e-10
        
    @partial(jit, static_argnums=(0,))
    def __call__(self, a, b):
        """Creates a Gaussian quadrature mesh for the interval (a, b) with n
        points.
        """

        # Initialize points and weights
        x_init = jnp.zeros((2, self.m))

        # Interval from (-1, 0)
        x_array, x_weights = fori_loop(1, self.m + 1, self.step_i, x_init)

        # Get full interval on (-1, 1)
        y_array = jnp.append(x_array, -jnp.flip(x_array[:-1]))
        y_weights = jnp.append(x_weights, jnp.flip(x_weights[:-1]))

        # Rescale between (a, b)
        y_array = y_array * (b - a) / 2 + (b + a) / 2
        y_weights = y_weights * (b - a) / 2
    
        return y_array, y_weights
    
    @partial(jit, static_argnums=(0,))
    def step_i(self, i, x):
        """Perform one step in the loop over i."""

        t = jnp.cos(jnp.pi * (i - 1 / 4) / (self.n + 1 / 2))
        t1 = 1

        # While loop on |t - t_1| >= \eps
        init_val = jnp.array([t1, t, 0])
        t1, t, pp = while_loop(self.cond_func, self.body_func, init_val)

        # Points x_i
        x = x.at[0, i - 1].set(-t)
        # Weights w_i
        x = x.at[1, i - 1].set(2 / ((1 - t ** 2) * pp ** 2))

        return x
    
    @partial(jit, static_argnums=(0,))
    def cond_func(self, loop_carry):
        """Evaluate the while loop condition |t - t_1| >= \eps."""
    
        t1, t, _ = loop_carry
    
        return jnp.abs(t - t1) >= self.eps
    
    @partial(jit, static_argnums=(0,))
    def body_func(self, loop_carry):
        """Perform one step in the while loop."""
    
        # Unpack t1, t, and pp
        t1, t, pp = loop_carry
    
        # Loop over j
        init = jnp.array([1, 0, t])
        p1, p2, t = fori_loop(1, self.n + 1, self.step_j, init)
    
        # Update pp, t1, and t
        pp = self.n * (t * p1 - p2) / (t ** 2 - 1)
        t1 = t
        t = t1 - p1 / pp

        return jnp.array([t1, t, pp])

    @partial(jit, static_argnums=(0,))
    def step_j(self, j, loop_carry):
        """Perform one step in the loop over j."""
    
        p2, p3, t = loop_carry
        p1 = ((2 * j - 1) * t * p2 - (j - 1) * p3) / j

        return jnp.array([p1, p2, t])