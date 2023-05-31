#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: clebsch_gordan.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 23, 2023

This script contains functions that compute Clebsch-Gordan coefficients.

Last update: May 25, 2023

"""

# Python imports
from functools import partial
from jax import config, jit, vmap
import jax.numpy as jnp
import numpy as np
from sympy.physics.quantum.cg import CG


# Enable double-precision with JAX arrays
config.update("jax_enable_x64", True)


def compute_clebsch_gordan_table(j_max):
    """
    Calculate Clebsch-Gordan coefficients for combinations of j and m_j up
    to j_max.
    
    Parameters
    ----------
    j_max : int
        Maximum j value for j_1, j_2, and j_3. This also constrains m_j.
    
    Returns
    -------
    cg_table : dict
        Table of Clebsch-Gordan coefficients <j_1 m_j_1 j_2 m_j_2|j_3 m_j_3>
        for each combination of angular momenta.
        
    """
        
    cg_table = {}
        
    j_array = np.arange(0, j_max+1/2, 1/2)
    
    for j_1 in j_array:
        for j_2 in j_array:
            j_3_array = np.arange(abs(j_1-j_2), j_1+j_2+1)
            for j_3 in j_3_array:
                for m_1 in np.arange(-j_1, j_1+1, 1):
                    for m_2 in np.arange(-j_2, j_2+1, 1):
                        
                        m_3 = m_1 + m_2
                        
                        if abs(m_3) <= j_3:
                            cg_table[(j_1,m_1,j_2,m_2,j_3,m_3)] = float(
                                CG(j_1,m_1,j_2,m_2,j_3,m_3).doit()
                            )
                                
    return cg_table


def clebsch_gordan_coefficient_v1(j1, m1, j2, m2, j3, m3, cg_table):
    """Clebsch-Gordan coefficient < j1 m1 j2 m2 | j3 m3 >."""
        
    # This try/except automatically enforces selection rules
    try:
        
        return cg_table[(j1, m1, j2, m2, j3, m3)]
                
    except KeyError:
        
        return 0
    

def compute_clebsch_gordan_array(j_max):
    """
    Calculate Clebsch-Gordan coefficients for combinations of j and m_j up
    to j_max.
    
    Parameters
    ----------
    j_max : int
        Maximum j value for j_1, j_2, and j_3. This also constrains m_j.
    
    Returns
    -------
    cg_table : 6-D JAX array
        Array of Clebsch-Gordan coefficients <j_1 m_j_1 j_2 m_j_2|j_3 m_j_3>
        for each combination of angular momenta.
        
    """
    
    j_array = np.arange(0, j_max+1/2, 1/2)
    N_j = j_array.size
    
    # 0, 1/2, 1, ..., J, -1/2, -1, ..., -J
    m_array = np.concatenate((j_array, -j_array[1:]))
    N_m = m_array.size
    
    cg_array = np.zeros((N_j, N_j, N_j, N_m, N_m, N_m))
    
    for i, j_1 in enumerate(j_array):
        m_1_array = np.arange(-j_1, j_1+1)
        for j, j_2 in enumerate(j_array):
            m_2_array = np.arange(-j_2, j_2+1)
            j_3_array = np.arange(np.abs(j_1-j_2), j_1+j_2+1)
            for k, j_3 in enumerate(j_array):
                m_3_array = np.arange(-j_3, j_3+1)
                for l, m_1 in enumerate(m_array):
                    for m, m_2 in enumerate(m_array):
                        for n, m_3 in enumerate(m_array):
                            
                            selection_rules = (
                                np.any(j_3 == j_3_array)
                                and np.any(m_1 == m_1_array)
                                and np.any(m_2 == m_2_array)
                                and np.any(m_3 == m_3_array)
                                and m_1 + m_2 == m_3
                            )

                            if selection_rules:

                                cg_array[i, j, k, l, m, n] = float(
                                    CG(j_1,m_1,j_2,m_2,j_3,m_3).doit()
                                )
    
    # return jnp.array(cg_array, dtype=jnp.float64), N_j
    return cg_array, N_j


@partial(jit, static_argnames=['N_j'])
def cg_mapping(j, m, N_j):
    """Return the indices of the input angular momentum and projection for the
    array of Clebsch-Gordan coefficients.
    """
    
    j_index = jnp.array(j / 0.5, dtype=int)
    m_index = jnp.array(jnp.abs(m/0.5) + jnp.heaviside(-m, 0) * (N_j-1),
                        dtype=int)

    return j_index, m_index


@partial(jit, static_argnames=['N_j'])
def clebsch_gordan_coefficient(j1, m1, j2, m2, j3, m3, cg_array, N_j):
    """Clebsch-Gordan coefficient < j1 m1 j2 m2 | j3 m3 >."""
    
    ij, im = cg_mapping(j1, m1, N_j)
    jj, jm = cg_mapping(j2, m2, N_j)
    kj, km = cg_mapping(j3, m3, N_j)
    
    return cg_array[ij, jj, kj, im, jm, km]


@partial(jit, static_argnames=['N_j'])
def clebsch_gordan_coefficient_vmap(j1, m1, j2, m2, j3, m3, cg_array, N_j):
    return vmap(
        clebsch_gordan_coefficient, in_axes=(0, 0, 0, 0, 0, 0, None, None),
        out_axes=(0)
    )(j1, m1, j2, m2, j3, m3, cg_array, N_j)