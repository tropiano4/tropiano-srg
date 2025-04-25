#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: current_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Class for evaluating the current operator written in JAX.

Last update: April 25, 2025

"""

# Python imports
from functools import partial

# JAX imports
from jax import config, jit, vmap
import jax.numpy as jnp

# Imports from scripts
from .clebsch_gordan import ClebschGordan
from .form_factors_jax import FormFactors
from .delta_U_jax import DeltaU
from .integration_jax import GaussQuadrature
from .special_functions import SpecialFunctions


# Enable double precision
config.update("jax_enable_x64", True)


class Current:
    """Momentum-space matrix element of the current operator J_0, where SRG
    evolution is split according to
        U_\lambda J_0 U_\lambda^\dagger = J_0
            + J_0 \delta U^\dagger
            + \delta U J_0
            + \delta U J_0 \delta U^\dagger
        = J_4 + J_3 + J_2 + J_1.
    Note, we use the J_4, J_3, J_2, and J_1 naming convention in this class.
    """
    
    def __init__(
            self, kvnn=6, kmax=25.0, kmid=4.0, ntot=120, lamb=jnp.inf, L_max=2
    ):
        """Initialize other classes and set possible m_s values."""
        
        # Load form factors from data files
        self.ff = FormFactors()
        
        # Special functions for Legendre polynomials
        self.sf = SpecialFunctions()
        
        # Clebsch-Gordan coefficients
        self.cg = ClebschGordan(L_max + 1)  # Set table up to j_max = L_max + 1
        
        # Possible m_s values
        self.m_s_array = jnp.array([-1, 0, 1])
        
        # Set possible L values for intermediate sums
        self.L_array = jnp.arange(0, L_max + 1, 1)
        
        # Set-up Gauss-quadrature integration mesh with ptot points
        ptot = 101
        self.gq = GaussQuadrature(ptot)
        
        # Set SRG \lambda as instance attribute
        self.lamb = lamb
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, k1_array, k2_array, q, J_1, m_J_d, L_1, L_2, T_1):
        """Matrix element of current operator assuming the form factors
        G_E^p = 1 and G_E^n = 0.
        """

        # Form factor part \pi^2 / 2 [G_E^p(Q^2) + (-1)^T_1 G_E^n(Q^2)]
        # Common in J_1, J_2, J_3, and J_4
        form_factors = jnp.pi ** 2 / 2

        # First vectorize over k_1
        matrix_element_vmap_k1 = vmap(
            self.matrix_element_wrt_k1_k2,
            in_axes=(0, None, None, None, None, None, None, None)
        )

        # Then vectorize over k_2 and plug in arrays of k_1 and k_2
        matrix_element_grid = vmap(
            matrix_element_vmap_k1,
            in_axes=(None, 0, None, None, None, None, None, None)
        )(k1_array, k2_array, q, J_1, m_J_d, L_1, L_2, T_1)

        # Units are fm^3 with shape (k_array.size, k_array.size)
        return form_factors * matrix_element_grid
    
    @partial(jit, static_argnums=(0,))
    def matrix_element_wrt_k1_k2(self, k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1):
        """Method for calling the unevolved or evolved current operator with
        respect to k_1 and k_2.
        """
        
        # Condition for no SRG evolution
        cond = self.lamb == jnp.inf
        
        # Return unevolved or evolved matrix element [fm^3]
        return jnp.where(
            cond,
            self.J4_wrt_k1_k2(k_1, k_2, q, J_1, m_J_d, L_1, L_2),
            self.J4_wrt_k1_k2(k_1, k_2, q, J_1, m_J_d, L_1, L_2)
            + self.J3_wrt_k1_k2(k_1, k_2, q, J_1, m_J_d, L_1, L_2)
            + self.J2_wrt_k1_k2(k_1, k_2, q, J_1, m_J_d, L_1, L_2)
            + self.J1_wrt_k1_k2(k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1)
        )
    
    @partial(jit, static_argnums=(0,))
    def J4_wrt_k1_k2(self, k_1, k_2, q, J_1, m_J_d, L_1, L_2):
        """Matrix element of the unevolved current operator for particular
        values of k_1 and k_2.
        """
        
        # Vectorize J_4 over m_s
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2
        matrix_element_array = vmap(
            self.J4_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J4_wrt_ms(self, m_s, args):
        """Matrix element of the unevolved current operator for particular
        values of k_1, k_2, and m_s.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2 = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # CG < J_1 m_J_d | L_1 m_J_d - m_s S = 1 m_s >
        cg_L1 = self.cg.get_coefficient(L_1, m_L, 1, m_s, J_1, m_J_d)
        
        # CG < L_2 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L2 = self.cg.get_coefficient(L_2, m_L, 1, m_s, 1, m_J_d)
        
        # Legendre polynomials
        k1_k2_plus_q = (k_1 ** 2 - k_2 ** 2 + q ** 2 / 4) / (k_1 * q)
        k1_k2_minus_q = (k_1 ** 2 - k_2 ** 2 - q ** 2 / 4) / (k_2 * q)
        P_L1 = self.sf.plm(k1_k2_plus_q, L_1, m_L)
        P_L2 = self.sf.plm(k1_k2_minus_q, L_2, m_L)
        
        # Units of matrix element are fm^3
        matrix_element = cg_L1 * cg_L2 * P_L1 * P_L2 * 2 / (k_1 * k_2 * q)
        
        # Non-zero only if |k_1 - q / 2 | < k_2 < k_1 + q / 2
        cond = jnp.logical_and(jnp.less(jnp.abs(k_1 - q / 2), k_2),
                               jnp.less(k_2, k_1 + q / 2))
        return jnp.where(cond, matrix_element, 0)
    
    @partial(jit, static_argnums=(0,))
    def J3_wrt_k1_k2(self, k_1, k_2, q, J_1, m_J_d, L_1, L_2):
        """Matrix element J_0 \delta U^\dagger for particular values of k_1 and
        k_2.
        """
        
        # Vectorize J_3 over m_s
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2
        matrix_element_array = vmap(
            self.J3_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J3_wrt_ms(self, m_s, args):
        """Matrix element J_0 \delta U^\dagger for particular values of k_1,
        k_2, and m_s.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2 = args
        
        # Vectorize J_3 over L_p
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2, m_s
        matrix_element_array = vmap(
            self.J3_wrt_Lp,
            in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_p
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J3_wrt_Lp(self, L_p, args):
        """Matrix element J_0 \delta U^\dagger for particular values of k_1,
        k_2, m_s, and L_p.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2, m_s = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # CG < J_1 m_J_d | L_1 m_J_d - m_s S = 1 m_s >
        cg_L1 = self.cg.get_coefficient(L_1, m_L, 1, m_s, J_1, m_J_d)
        
        # CG < L_p m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_Lp = self.cg.get_coefficient(L_p, m_L, 1, m_s, 1, m_J_d)
        
        # Integrate over p from [|k_1 - q / 2], k_1 + q / 2]
        p_min = jnp.abs(k_1 - q / 2)
        p_max = k_1 + q / 2
        # p_array = jnp.linspace(p_min, p_max, self.ptot)
        p_array, p_weights = self.gq(p_min, p_max)
        
        # Legendre polynomials
        k1_p_plus_q = (k_1 ** 2 - p_array ** 2 + q ** 2 / 4) / (k_1 * q)
        # Make sure (k_1 ^ 2 - p^2 - q^2 / 4) / (p q) is within [-1, 1]
        k1_p_minus_q = jnp.clip(
            (k_1 ** 2 - p_array ** 2 - q ** 2 / 4) / (p_array * q), -1, 1
        )
        P_L1_array = self.sf.plm(k1_p_plus_q, L_1, m_L)
        P_Lp_array = self.sf.plm(k1_p_minus_q, L_p, m_L)

        # \delta U^\dagger matrix element where S = 1, J = 1, and T = 0
        delta_U_dagger = self.deltaU(p_array, k_2, 1, L_p, L_2, 1, 0, hc=True)
        
        # Set integrand over p
        integrand = (p_array ** 2 * P_L1_array * P_Lp_array
                     * 2 / (k_1 * p_array * q) * delta_U_dagger)
        
        # Integrate over p
        matrix_element = 2 / jnp.pi * cg_L1 * cg_Lp * jnp.sum(p_weights
                                                              * integrand)

        # Units of matrix element are fm^3
        return matrix_element
    
    @partial(jit, static_argnums=(0,))
    def J2_wrt_k1_k2(self, k_1, k_2, q, J_1, m_J_d, L_1, L_2):
        """Matrix element \delta U J_0 for particular values of k_1 and k_2."""

        # Vectorize J_2 over m_s
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2
        matrix_element_array = vmap(
            self.J2_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J2_wrt_ms(self, m_s, args):
        """Matrix element \delta U J_0 for particular values of k_1, k_2, and
        m_s.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2 = args
        
        # Vectorize J_2 over L_p
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2, m_s
        matrix_element_array = vmap(
            self.J2_wrt_Lp, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_p
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J2_wrt_Lp(self, L_p, args):
        """Matrix element \delta U J_0 for particular values of k_1, k_2, m_s,
        and L_p.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2, m_s = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # CG < J_1 m_J_d | L_p m_J_d - m_s S = 1 m_s >
        cg_Lp = self.cg.get_coefficient(L_p, m_L, 1, m_s, J_1, m_J_d)
        
        # CG < L_2 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L2 = self.cg.get_coefficient(L_2, m_L, 1, m_s, 1, m_J_d)
        
        # Integrate over p from [|k_2 - q / 2], k_2 + q / 2]
        p_min = jnp.abs(k_2 - q / 2)
        p_max = k_2 + q / 2
        # p_array = jnp.linspace(p_min, p_max, self.ptot)
        p_array, p_weights = self.gq(p_min, p_max)
        
        # Legendre polynomials
        # Make sure (p^2 - k_2^2 + q^2 / 4) / (p q) is within [-1, 1]
        k2_p_plus_q = jnp.clip(
            (p_array ** 2 - k_2 ** 2 + q ** 2 / 4) / (p_array * q), -1, 1
        )
        k2_p_minus_q = (p_array ** 2 - k_2 ** 2 - q ** 2 / 4) / (k_2 * q)
        P_Lp_array = self.sf.plm(k2_p_plus_q, L_p, m_L)
        P_L2_array = self.sf.plm(k2_p_minus_q, L_2, m_L)

        # \delta U matrix element where S = 1, J = 1, and T = 0
        delta_U = self.deltaU(k_1, p_array, 1, L_1, L_p, 1, 0)
        
        # Set integrand over p
        integrand = (p_array ** 2 * P_Lp_array * P_L2_array
                     * 2 / (k_2 * p_array * q) * delta_U)
        
        # Integrate over p
        matrix_element = 2 / jnp.pi * cg_Lp * cg_L2 * jnp.sum(p_weights
                                                              * integrand)

        # Units of matrix element are fm^3
        return matrix_element
    
    @partial(jit, static_argnums=(0,))
    def J1_wrt_k1_k2(self, k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1):
        """Matrix element \delta U J_0 \delta U^\dagger for particular values
        of k_1 and k_2.
        """

        # Vectorize J_1 over m_s
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1
        matrix_element_array = vmap(
            self.J1_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J1_wrt_ms(self, m_s, args):
        """Matrix element \delta U J_0 \delta U^\dagger for particular values
        of k_1, k_2, and m_s.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1 = args
        
        # Vectorize J_1 over L_p_1
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1, m_s
        matrix_element_array = vmap(
            self.J1_wrt_Lp1, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_p_1
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J1_wrt_Lp1(self, L_p_1, args):
        """Matrix element \delta U J_0 \delta U^\dagger for particular values
        of k_1, k_2, m_s, and L_p_1.
        """
        
        #  Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1, m_s = args
        
        # Vectorize J_1 over L_p_2
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1, m_s, L_p_1
        matrix_element_array = vmap(
            self.J1_wrt_Lp2, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_p_2
        return jnp.sum(matrix_element_array)
    
    @partial(jit, static_argnums=(0,))
    def J1_wrt_Lp2(self, L_p_2, args):
        """Matrix element \delta U J_0 \delta U^\dagger for particular values
        of k_1, k_2, m_s, L_p_1, and L_p_2.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1, m_s, L_p_1 = args
        
        # CG < J_1 m_J_d | L_p_1 m_J_d - m_s S = 1 m_s >
        cg_Lp1 = self.cg.get_coefficient(L_p_1, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_p_2 m_J_d - m_s S = 1 m_s | 1 m_J_d >
        cg_Lp2 = self.cg.get_coefficient(L_p_2, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Integrate p_1 from 0 to kmax
        p_1_array, p_1_weights = self.gq(0.0, self.deltaU.kmax)
        
        # Vectorize J_1 over p_1
        args = k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1, m_s, L_p_1, L_p_2
        matrix_element_array = vmap(
            self.J1_wrt_p1, in_axes=(0, None)
        )(p_1_array, args)
        
        # Integrand w.r.t. p_1
        integrand = ((2 / jnp.pi) ** 2 * cg_Lp1 * cg_Lp2 * p_1_array ** 2
                     * matrix_element_array)
        
        # Integrate over p_1
        return jnp.sum(p_1_weights * integrand)
    
    @partial(jit, static_argnums=(0,))
    def J1_wrt_p1(self, p_1, args):
        """Matrix element \delta U J_0 \delta U^\dagger for particular values
        of k_1, k_2, m_s, L_p_1, L_p_2, and p_1.
        """
        
        # Unpack other arguments
        k_1, k_2, q, J_1, m_J_d, L_1, L_2, T_1, m_s, L_p_1, L_p_2 = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s

        # Integrate over p_2 from [|p_1 - q / 2], p_1 + q / 2]
        p_2_min = jnp.abs(p_1 - q / 2)
        p_2_max = p_1 + q / 2
        p_2_array, p_2_weights = self.gq(p_2_min, p_2_max)
        
        # Legendre polynomials
        # Make sure momenta arguments are within [-1, 1]
        p1_p2_plus_q = jnp.clip(
            (p_1 ** 2 - p_2_array ** 2 + q ** 2 / 4) / (p_1 * q), -1, 1
        )
        p1_p2_minus_q = jnp.clip(
            (p_1 ** 2 - p_2_array ** 2 - q ** 2 / 4) / (p_2_array * q), -1, 1
        )
        P_Lp1_array = self.sf.plm(p1_p2_plus_q, L_p_1, m_L)
        P_Lp2_array = self.sf.plm(p1_p2_minus_q, L_p_2, m_L)
        
        # \delta U matrix elements
        delta_U = self.deltaU(k_1, p_1, J_1, L_1, L_p_1, 1, T_1)
        delta_U_dagger = self.deltaU(p_2_array, k_2, 1, L_p_2, L_2, 1, 0,
                                     hc=True)
        
        # Set integrand over p_2
        integrand = (
            p_2_array ** 2 * P_Lp1_array * P_Lp2_array
            * 2 / (p_1 * p_2_array * q) * delta_U * delta_U_dagger
        )
        
        # Integrate over p_2
        return jnp.sum(p_2_weights * integrand)

