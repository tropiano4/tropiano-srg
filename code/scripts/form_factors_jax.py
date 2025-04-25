#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: form_factors_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Class for loading and interpolating electric proton and neutron form factors
written in JAX.

Last update: April 25, 2025

"""

# Python imports
from functools import partial
import numpy as np

# JAX imports
from jax import config, jit
import jax.numpy as jnp


# Enable double precision
config.update("jax_enable_x64", True)


class FormFactors:
    """Form factors from Sushant's data files."""

    def __init__(self):
        """Interpolate data files."""
        
        gep_data = np.loadtxt("../data/form_factors/gep.dat")
        gen_data = np.loadtxt("../data/form_factors/gen.dat")
        
        # Set data as instance attributes
        self.Q2_gep_GeV = jnp.asarray(gep_data[:, 0])
        self.gep_array = jnp.asarray(gep_data[:, 1])
        self.Q2_gen_GeV = jnp.asarray(gen_data[:, 0])
        self.gen_array = jnp.asarray(gen_data[:, 1])

    @partial(jit, static_argnums=(0,))
    def GEp(self, Q2):
        """Electric proton form factor w.r.t. Q^2 in GeV^2."""

        return jnp.interp(Q2, self.Q2_gep_GeV, self.gep_array) * self.GD(Q2)
    
    @partial(jit, static_argnums=(0,))
    def GEn(self, Q2):
        """Electric neutron form factor w.r.t. Q^2 in GeV^2."""
        
        return jnp.interp(Q2, self.Q2_gen_GeV, self.gen_array)
    
    @partial(jit, static_argnums=(0,))
    def GD(self, Q2):
        """Dipole form factor w.r.t. Q^2 in GeV^2."""
        
        mD2 = 0.71  # GeV^2
        
        return (1 + Q2 / mD2) ** -2

