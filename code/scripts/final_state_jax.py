#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: final_state_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Class for computing the final state wavefunction in NN scattering written in 
JAX.

Last update: April 25, 2025

"""

# Python imports
from functools import partial
import numpy as np

# JAX imports
from jax import config, jit
import jax.numpy as jnp

# Imports from scripts
from .clebsch_gordan import ClebschGordan
from .tmatrix_jax import TMatrix


# Enable double precision
config.update("jax_enable_x64", True)


class FinalState:
    """
    Final-state wave function without plane-wave term evaluated at k_i != p'.
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb):
        
        # Set table of Clebsch-Gordan coefficients with j_max = 2
        self.cg = ClebschGordan(2)
        
        # Initialize T-matrix class
        self.tmatrix = TMatrix(kvnn, kmax, kmid, ntot, lamb, L_max=2)

    @partial(jit, static_argnums=(0,))
    def delta_psi(self, pp):

        # Not sure about these!
        m_J_d = 0
        m_S_f = 0
        
        k_array = self.tmatrix.k_array
        delta_psi = jnp.zeros(self.tmatrix.ntot, dtype=complex)
        
        # Denominator of Green's function without imaginary part
        greens_func = 1 / (pp ** 2 - k_array ** 2)
        
        # Sum over L' = 0, 2 with L, S, and J fixed
        J, L, S, T = 1, 0, 1, 0

        # Half off-shell T-matrix in units [fm]
        T_0_array = self.tmatrix.half_offshell(pp, J, L, 0, S, T)  # L' = 0
        T_2_array = self.tmatrix.half_offshell(pp, J, L, 2, S, T)  # L' = 2
            
        # Look-up Clebsch-Gordan coefficients
        cg_0 = self.cg.get_coefficient(0, m_J_d - m_S_f, S, m_S_f, J, m_J_d)
        cg_2 = self.cg.get_coefficient(2, m_J_d - m_S_f, S, m_S_f, J, m_J_d)
            
        # Extra factor of 2 from [1 + (-1)^T (-1) L']
        delta_psi = 0.5 * jnp.sqrt(2/np.pi) * greens_func * 2 * (
            T_0_array * cg_0 + T_2_array * cg_2
        )
   
        return k_array, delta_psi

