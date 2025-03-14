#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: clebsch_gordan.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: March 7, 2025

Class for computing Clebsch-Gordan coefficients using JAX and SymPy.

Last update: March 7, 2025

"""

# Python imports
from functools import partial
import numpy as np
from sympy.physics.quantum.cg import CG

# JAX imports
import jax.numpy as jnp
from jax import config, jit


# Enable double precision
config.update("jax_enable_x64", True)


# TODO: Enable half-integer values of j and m_j
class ClebschGordan:
    """
    Class for computing Clebsch-Gordan coefficients using SymPy to precompute
    and storing values in a JAX array for look-up. Note, we do not need to
    worry about half-integer values here.
    """
    
    def __init__(self, j_max):
        
        # Precompute the JAX array of Clebsch-Gordan coefficients
        self.cg_array = self.compute_clebsch_gordan_array(j_max)
        
        # Set j_max as an instance attribute
        self.j_max = j_max
    
    def compute_clebsch_gordan_array(self, j_max):
        """Calculate Clebsch-Gordan coefficients for combinations of j and m_j
        up to j_max.
        """
        
        # Initialize big JAX array
        ntot_j = j_max + 1
        ntot_m = 2 * ntot_j + 1
        cg_array = jnp.zeros((ntot_j, ntot_m, ntot_j, ntot_m, ntot_j, ntot_m))
        
        # Loop over quantum numbers and compute/store coefficients
        j_array = np.arange(0, j_max + 1, 1)
        m_array = np.arange(-j_max, j_max + 1, 1)
        for ij1, j_1 in enumerate(j_array):
            for im1, m_1 in enumerate(m_array):
                for ij2, j_2 in enumerate(j_array):
                    for im2, m_2 in enumerate(m_array):
                        for ij3, j_3 in enumerate(j_array):
                            for im3, m_3 in enumerate(m_array):
                                
                                # Check if the quantum numbers make sense
                                cond = self.is_physical(j_1, m_1, j_2, m_2, j_3,
                                                        m_3)
                                
                                # Value of coefficient will be zero otherwise
                                if cond:
                                    
                                    cg_array = (
                                        cg_array.at[
                                            ij1, im1, ij2, im2, ij3, im3
                                        ].set(
                                            float(
                                                CG(
                                                    j_1, m_1, j_2, m_2, j_3, m_3
                                                ).doit()
                                            )
                                        )
                                    )
                                    
        return cg_array
    
    @partial(jit, static_argnums=(0,))
    def is_physical(self, j_1, m_1, j_2, m_2, j_3, m_3):
        """Check if quantum numbers make sense."""
        
        m_1_bool = self.check_projection(j_1, m_1)
        m_2_bool = self.check_projection(j_2, m_2)
        # j_1 - j_2 <= j_3 <= j_1 + j_2 condition
        j_3_bool = jnp.logical_and(
            jnp.less_equal(jnp.abs(j_1 - j_2), j_3),
            jnp.less_equal(j_3, j_1 + j_2)
        )
        # m_3 = m_1 + m_2 and |m_3| <= j_3 conditions
        m_3_bool = jnp.logical_and(
            jnp.equal(m_3, m_1 + m_2),
            self.check_projection(j_3, m_3)
        )
        
        return jnp.logical_and(
            m_1_bool, jnp.logical_and(
                m_2_bool, jnp.logical_and(j_3_bool, m_3_bool)
            )
        )
    
    @partial(jit, static_argnums=(0,))
    def check_projection(self, j, m_j):
        """Check if angular momentum projection is physical."""
        
        # |m_j| <= j condition
        return jnp.less_equal(jnp.abs(m_j), j)

    @partial(jit, static_argnums=(0,))
    def get_coefficient(self, j_1, m_1, j_2, m_2, j_3, m_3):
        """Return the CG coefficient given input angular momentum."""
        
        # Values of j are the same as its index
        ij1, ij2, ij3 = j_1, j_2, j_3
        
        # Mapping from m_j to indices
        im1, im2, im3 = m_1 + self.j_max, m_2 + self.j_max, m_3 + self.j_max

        # Make sure quantum numbers make sense
        cond = self.is_physical(j_1, m_1, j_2, m_2, j_3, m_3)
        return jnp.where(cond, self.cg_array[ij1, im1, ij2, im2, ij3, im3], 0)