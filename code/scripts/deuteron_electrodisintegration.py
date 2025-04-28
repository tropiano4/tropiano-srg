#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: deuteron_electrodisintegration.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Several classes for computing the deuteron electrodisintegration longitudinal
structure function in JAX.

Last update: April 25, 2025

"""

# Python imports
from functools import partial

# JAX imports
from jax import config, jit, vmap
from jax.lax import fori_loop
import jax.numpy as jnp
from jax.scipy.special import sph_harm

# Imports from scripts
from .clebsch_gordan import ClebschGordan
from .delta_U_jax import DeltaU
from .deuteron_wf_jax import DeuteronWaveFunction
from .form_factors_jax import FormFactors
from .integration_jax import GaussQuadrature
from .special_functions import SpecialFunctions
from .tmatrix_jax import TMatrix


# Enable double precision
config.update("jax_enable_x64", True)


# TODO: Get rid of option 5 when things are working
# TODO: Try vmap over QFR again?
class DeuteronElectrodisintegration:
    """Class that calculates the longitudinal structure function for deuteron
    electrodisintegration.
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb=jnp.inf, L_max=2, option=1):

        # Constants will be used as instance attributes
        self.hbar_c = 197.32696  # \hbar c [MeV fm]
        M_n = 939.56563  # Neutron mass [MeV]
        M_p = 938.27231  # Proton mass [MeV]
        self.M = (M_p + M_n) / 2  # Nucleon mass [MeV]
        B_d = 2.224  # Binding energy of deuteron [MeV]
        self.M_d = 2 * self.M - B_d  # Deuteron mass [MeV]
        self.alpha = 1 / 137.03599  # Fine structure constant
        
        # Quantum numbers for longitudinal structure function S_f, m_S_f, m_J_d
        self.fL_quantum_numbers = self.fL_sum()

        # Load form factors from data files
        self.ff = FormFactors()

        # Classes required for overlaps: special functions and Clebsch-Gordan
        # coefficients
        sf = SpecialFunctions()
        cg = ClebschGordan(L_max + 1)  # Set table up to j_max = L_max + 1
        
        # Set overlap function based on option
        self.set_overlap(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg, option)
        
    def set_overlap(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg, option):
        """This method sets-up the overlap function based on the option.
            Option 1: Impulse approximation with no SRG evolution -> B_4
            Option 2: Impulse approximation with SRG evolution -> A + B
            Option 3: Include FSI with no SRG evolution -> B_4 + F_4
            Option 4: Include FSI with SRG evolution -> B + F
        """
        
        # All options include B_4
        b4 = B4(kvnn, kmax, kmid, ntot, lamb, L_max, cg)
        
        # Option 1: B_4 only
        if option == 1:
            
            self.overlap = b4
        
        # Option 2: A + B
        elif option == 2:
            
            b3 = B3(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            b2 = B2(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            b1 = B1(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            a4 = A4(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            a3 = A3(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            a2 = A2(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            a1 = A1(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            
            self.overlap = (
                lambda pp, thetap, q, gep, gen, fL_quantum_numbers:
                    b4(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + b3(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + b2(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + b1(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + a4(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + a3(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + a2(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + a1(pp, thetap, q, gep, gen, fL_quantum_numbers)
            )
            
        # Option 3: B_4 + F_4
        elif option == 3:
            
            f4 = F4(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            
            self.overlap = (
                lambda pp, thetap, q, gep, gen, fL_quantum_numbers:
                    b4(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + f4(pp, thetap, q, gep, gen, fL_quantum_numbers)
            )
            
        # Option 4: B + F
        elif option == 4:
            
            b3 = B3(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            b2 = B2(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            b1 = B1(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            f4 = F4(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            f3 = F3(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            f2 = F2(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            f1 = F1(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            
            self.overlap = (
                lambda pp, thetap, q, gep, gen, fL_quantum_numbers:
                    b4(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + b3(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + b2(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + b1(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + f4(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + f3(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + f2(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + f1(pp, thetap, q, gep, gen, fL_quantum_numbers)
            )
            
        # TESTING
        elif option == 5:
            
            self.f1 = F1(kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg)
            self.overlap = (
                lambda pp, thetap, q, gep, gen, fL_quantum_numbers:
                    b4(pp, thetap, q, gep, gen, fL_quantum_numbers)
                    + self.f1(pp, thetap, q, gep, gen, fL_quantum_numbers)
            )

    def fL_sum(self):
        """Repackage sum over S_f, m_S_f, and m_J_d into a JAX array."""
        
        fL_list = []
        
        # Possible values of S_f and m_J_d
        S_f_array = jnp.array([0, 1])
        m_J_d_array = jnp.array([-1, 0, 1])
        
        # Sum over S_f, m_S_f, and m_J_d
        for S_f in S_f_array:
            
            m_S_f_array = jnp.arange(-S_f, S_f + 1, 1)
            for m_S_f in m_S_f_array:
                
                for m_J_d in m_J_d_array:
                    
                    # Nothing depends on S_f!
                    fL_list.append([m_S_f, m_J_d])

        # Return as JAX array
        return jnp.asarray(fL_list)
    
    @partial(jit, static_argnums=(0,))
    def fL_vmap_theta(self, Ep, thetap_array, q):
        """Longitudinal structure function f_L [fm] with respect to scalars E'
        [MeV] and q [fm^-1] and array \theta' [deg].
        """

        # Vectorize fL method over theta'
        return vmap(self.fL, in_axes=(None, 0, None))(Ep, thetap_array, q)

    @partial(jit, static_argnums=(0,))
    def fL_quasifree_ridge(self, Ep_array, thetap):
        """Longitudinal structure function f_L [fm] with respect to scalar
        \theta' [deg] and array E' [MeV] where \omega = 0.
        """
        
        # Initialize f_L array over E'
        fL_init_array = jnp.zeros_like(Ep_array)
        
        # Set E' and \theta' as instance attributes
        self.Ep_array = Ep_array
        self.thetap = thetap
    
        # Loop over E' to compute f_L
        fL_array = fori_loop(0, Ep_array.size, self.quasifree_ridge_step,
                             fL_init_array)
        
        return fL_array
    
    @partial(jit, static_argnums=(0,))
    def quasifree_ridge_step(self, i, fL_array):
        """Perform one step over E' to compute f_L along quasifree ridge."""
            
        Ep = self.Ep_array[i]
        
        # Quasifree ridge where \omega = 0 determines q [fm^-1]
        q = jnp.sqrt(
            (Ep + 2 * self.M) ** 2 - self.M_d ** 2
        ) / self.hbar_c

        # Add f_L(E'=E[i]) to array
        fL_array = fL_array.at[i].set(self.fL(Ep, self.thetap, q))
            
        return fL_array
    
    @partial(jit, static_argnums=(0,))
    def fL(self, Ep, thetap, q):
        """Longitudinal structure function f_L [fm] where the arguments are E'
        [MeV], \theta' [deg], and q [fm^-1].
        """
        
        # Convert \theta' to radians
        thetap_radians = jnp.radians(thetap)

        # Outgoing proton momentum [fm^-1]
        pp = jnp.sqrt(self.M * Ep + Ep ** 2 / 4) / self.hbar_c
        
        # Nucleon energy [MeV]
        E_N = jnp.sqrt((self.hbar_c * pp) ** 2 + self.M ** 2)
        
        # Virtual momentum q in units [MeV]
        q_MeV = self.hbar_c * q
        
        # Deuteron energy [MeV]
        E_d = jnp.sqrt(q_MeV ** 2 + self.M_d ** 2)
        
        # Energy of virtual photon [MeV]
        omega = Ep - E_d + 2 * self.M
        
        # Four momentum transfer squared [GeV^2]
        Q2 = (q_MeV ** 2 - omega ** 2) * 1e-6
        
        # Compute form factors [unitless]
        gep = self.ff.GEp(Q2)
        gen = self.ff.GEn(Q2)

        # Kinematic factor [fm^-1]
        factor = -jnp.pi * jnp.sqrt(2 * self.alpha * pp * E_N * E_d / self.M_d
                                    / self.hbar_c)

        # Vectorize amplitude over all S_f, m_S_f, and m_J_d combinations
        amplitude_array = factor * vmap(
            self.overlap, in_axes=(None, None, None, None, None, 0)
        )(pp, thetap_radians, q, gep, gen, self.fL_quantum_numbers)
        
        # Sum amplitude absolute value squared over S_f, m_S_f, and m_J_d
        f_L = jnp.sum(jnp.abs(amplitude_array) ** 2)

        # Return longitudinal structure function [fm] which is a scalar
        return f_L
    
    
class B4:
    """Class for calculating the overlap matrix element in the impulse
    approximation.
        < \phi | J_0 | \psi_i(\lambda) >
    Note, \lambda sets the SRG evolution of the initial state, where
    \lambda = \infty is unevolved.
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, cg):

        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Array for sum over L_d
        self.L_d_array = jnp.array([0, 2])

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers

        # Vectorize overlap method over L_d
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_Ld, in_axes=(0, None)
        )(self.L_d_array, args)
        
        # Sum over L_d
        overlap = jnp.sum(overlap_array)
        
        # Return overlap [fm^3/2] which is a scalar
        return overlap

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_Ld(self, L_d, args):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        m_J_d, and L_d.
        """
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_S_f

        # Dot product p'_vector \dot q_vector [fm^-2]
        pqx = pp * q * jnp.cos(thetap)
        
        # Magnitude of momenta |p'_vector +/- q_vector / 2| [fm^-1]
        pq_minus = jnp.sqrt(pp ** 2 + q ** 2 / 4 - pqx)
        pq_plus = jnp.sqrt(pp ** 2 + q ** 2 / 4 + pqx)
        
        # Angles between the unit vector z^\hat and p'_vector -/+ q_vector / 2
        theta_minus = jnp.arccos((pp * jnp.cos(thetap) - q / 2) / pq_minus)
        theta_plus = jnp.arccos((pp * jnp.cos(thetap) + q / 2) / pq_plus)
        
        # CG < L_d m_J_d - m_S_f S = 1 m_S_f | J = 1 m_J_d >
        cg = self.cg.get_coefficient(L_d, m_L, 1, m_S_f, 1, m_J_d)
                
        # Spherical harmonics
        theta_array = jnp.array([theta_minus, theta_plus])
        phip_array = jnp.zeros_like(theta_array)
        # Quantum numbers must have same shape as \theta and \phi
        L_d_array = jnp.repeat(L_d, 2)
        m_L_array = jnp.repeat(m_L, 2)
        Ylm_array = sph_harm(m_L_array, L_d_array, phip_array, theta_array,
                             n_max=self.L_max)
        Ylm_minus, Ylm_plus = Ylm_array
 
        # Deuteron wave function [fm^3/2]
        psi_minus = self.dwf.compute_wf(pq_minus, L_d)
        psi_plus = self.dwf.compute_wf(pq_plus, L_d)

        # Compute overlap [fm^3/2] which is a scalar
        overlap = jnp.sqrt(2 / jnp.pi) * cg * (
            gep * psi_minus * Ylm_minus + gen * psi_plus * Ylm_plus
        )
                
        return overlap
    

class B3:
    """Class for calculating the overlap matrix element:
    < \phi | J_0 \delta U^\dagger | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):
        
        # Set L_max as instance attribute
        self.L_max = L_max

        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])

        # Set cos(\theta) integration mesh
        ntot_theta = 21
        theta_array, theta_weights = GaussQuadrature(ntot_theta)(0, jnp.pi)

        # Meshgrids and Jacobian for \theta and k_3
        self.theta_grid, self.k3_grid = jnp.meshgrid(
            theta_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk3_grid = jnp.meshgrid(
            theta_weights, self.dwf.k_weights, indexing='ij'
        )
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk3_grid
                         * self.k3_grid ** 2)
        
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for m_s in m_s_array:
                        for L_2 in L_array:
                            for L_d in L_d_array:
                                        
                                # Check if the partial wave channels are
                                # physical
                                deltaU_bool = self.deltaU.channel_is_physical(
                                    1, L_d, L_2, 1, 0
                                )

                                # Check if T_1 and L_1 factor is 0
                                if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                    lt_bool = False
                                else:
                                    lt_bool = True
                                        
                                # Append combination to quantum numbers
                                if deltaU_bool and lt_bool:
                                    quantum_numbers.append([T_1, L_1, J_1, m_s,
                                                            L_2, L_d])
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, m_s, L_2, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # CG < J_1 m_J_d | L_1 m_J_d - m_s S = 1 m_s >
        cg_L1_ms = self.cg.get_coefficient(L_1, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_2 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Dot product p'_vector \dot q_vector
        ppqx_grid = pp * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |p'_vector - q_vector / 2|
        ppq_minus_grid = jnp.sqrt(pp ** 2 + q ** 2 / 4 - ppqx_grid)
        # Angle between the unit vector z^\hat and p'_vector - q_vector / 2
        cos_alphap_pp = ((pp * jnp.cos(self.theta_grid) - q / 2)
                         / ppq_minus_grid)
        
        # Legendre polynomials
        P_L1_grid = self.sf.plm(jnp.cos(self.theta_grid), L_1, m_J_d - m_s)
        P_L2_grid = self.sf.plm(cos_alphap_pp, L_2, m_J_d - m_s)
        
        # Deuteron wave function [fm^3/2]
        psi_k3_grid = self.dwf.compute_wf(self.k3_grid, L_d)
        
        # \delta U_{J = 1, L_d, L_2, S = 1, T = 0}(k_3, |p' - q / 2|)
        delta_U_grid = self.deltaU(self.k3_grid, ppq_minus_grid, 1, L_d, L_2, 1,
                                   0)

        # Integrate over \theta and k_3
        integrand = P_L1_grid * P_L2_grid * psi_k3_grid * delta_U_grid
        integral = jnp.sum(self.jacobian * integrand)

        # Return overlap
        return 2 * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L1_ms
            * cg_L2_ms * integral
        )

    
class B2:
    """Class for calculating the overlap matrix element:
    < \phi | \delta U J_0 | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):

        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])

        # Set cos(\theta) integration mesh
        ntot_theta = 21
        theta_array, theta_weights = GaussQuadrature(ntot_theta)(0, jnp.pi)

        # Meshgrids and Jacobian for \theta and k_2
        self.theta_grid, self.k2_grid = jnp.meshgrid(
            theta_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk2_grid = jnp.meshgrid(
            theta_weights, self.dwf.k_weights, indexing='ij'
        )
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk2_grid
                         * self.k2_grid ** 2)
        
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for m_s in m_s_array:
                        for L_2 in L_array:
                            for L_d in L_d_array:
                                        
                                # Check if the partial wave channels are
                                # physical
                                deltaU_bool = self.deltaU.channel_is_physical(
                                    J_1, L_1, L_2, 1, T_1
                                )

                                # Check if T_1 and L_1 factor is 0
                                if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                    lt_bool = False
                                else:
                                    lt_bool = True
                                        
                                # Append combination to quantum numbers
                                if deltaU_bool and lt_bool:
                                    quantum_numbers.append([T_1, L_1, J_1, m_s,
                                                            L_2, L_d])
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, m_s, L_2, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # CG < J_1 m_J_d | L_2 m_J_d - m_s S = 1 m_s >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_d m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_Ld_ms = self.cg.get_coefficient(L_d, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Dot product k2_vector \dot q_vector
        k2qx_grid = self.k2_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k2_vector - q_vector / 2|
        k2q_minus_grid = jnp.sqrt(self.k2_grid ** 2 + q ** 2 / 4 - k2qx_grid)
        # Angle between the unit vector z^\hat and k2_vector - q_vector / 2
        cos_alphap_k2 = ((self.k2_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k2q_minus_grid)
        
        # Legendre polynomials
        P_L2_grid = self.sf.plm(jnp.cos(self.theta_grid), L_2, m_J_d - m_s)
        P_Ld_grid = self.sf.plm(cos_alphap_k2, L_d, m_J_d - m_s)
        
        # Deuteron wave function [fm^3/2]
        psi_k2q_grid = self.dwf.compute_wf(k2q_minus_grid, L_d)
        
        # \delta U_{J_1, L_1, L_2, J_1, S = 1, T_1}(p', k_2)
        delta_U_grid = self.deltaU(pp, self.k2_grid, J_1, L_1, L_2, 1, T_1)

        # Integrate over \theta and k_2
        integrand = P_L2_grid * P_Ld_grid * psi_k2q_grid * delta_U_grid
        integral = jnp.sum(self.jacobian * integrand)

        # Return overlap
        return 2 * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L2_ms
            * cg_Ld_ms * integral
        )


class B1:
    """Class for calculating the overlap matrix element:
    < \phi | \delta U J_0 \delta U^\dagger | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):

        # Set L_max as instance attribute
        self.L_max = L_max

        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])

        # Set cos(\theta) integration mesh
        ntot_theta = 21
        theta_array, theta_weights = GaussQuadrature(ntot_theta)(0, jnp.pi)

        # Meshgrids and Jacobian for \theta, k_2, and k_4
        self.theta_grid, self.k2_grid, self.k4_grid = jnp.meshgrid(
            theta_array, self.dwf.k_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk2_grid, dk4_grid = jnp.meshgrid(
            theta_weights, self.dwf.k_weights, self.dwf.k_weights, indexing='ij'
        )
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk2_grid
                         * self.k2_grid ** 2 * dk4_grid * self.k4_grid ** 2)
        
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for m_s in m_s_array:
                        for L_2 in L_array:
                            for L_3 in L_array:
                                for L_d in L_d_array:
                                        
                                    # Check if the partial wave channels are
                                    # physical
                                    deltaU_12_bool = (
                                        self.deltaU.channel_is_physical(
                                            J_1, L_1, L_2, 1, T_1
                                        )
                                    )
                                    deltaU_d3_bool = (
                                        self.deltaU.channel_is_physical(
                                            1, L_d, L_3, 1, 0
                                        )
                                    )
                                    channel_bool = (deltaU_12_bool
                                                    and deltaU_d3_bool)
                                        
                                    # Check if T_1 and L_1 factor is 0
                                    if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                        lt_bool = False
                                    else:
                                        lt_bool = True
                                        
                                    # Append combination to quantum numbers
                                    if channel_bool and lt_bool:
                                        quantum_numbers.append(
                                            [T_1, L_1, J_1, m_s, L_2, L_3, L_d]
                                        )
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, m_s, L_2, L_3, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # CG < J_1 m_J_d | L_2 m_J_d - m_s S = 1 m_s >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_3 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, 1, m_J_d)

        # Dot product k2_vector \dot q_vector
        k2qx_grid = self.k2_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k2_vector - q_vector / 2|
        k2q_minus_grid = jnp.sqrt(self.k2_grid ** 2 + q ** 2 / 4 - k2qx_grid)
        # Angle between the unit vector z^\hat and k2_vector - q_vector / 2
        cos_alphap_k2 = ((self.k2_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k2q_minus_grid)
        
        # Legendre polynomials
        P_L2_grid = self.sf.plm(jnp.cos(self.theta_grid), L_2, m_J_d - m_s)
        P_L3_grid = self.sf.plm(cos_alphap_k2, L_3, m_J_d - m_s)
        
        # Deuteron wave function [fm^3/2]
        psi_k4_grid = self.dwf.compute_wf(self.k4_grid, L_d)
        
        # \delta U_{J_1, L_1, L_2, S = 1, T_1}(p', k_2)
        delta_U_L1_L2_grid = self.deltaU(pp, self.k2_grid, J_1, L_1, L_2, 1,
                                         T_1)
        
        # \delta U_{J = 1, L_d, L_3, S = 1, T = 0}(k_4, |k_2 - q / 2|)
        delta_U_Ld_L3_grid = self.deltaU(self.k4_grid, k2q_minus_grid, 1, L_d,
                                         L_3, 1, 0)
        
        # Integrate over \theta, k_2, and k_4
        integrand = (P_L2_grid * P_L3_grid * psi_k4_grid * delta_U_L1_L2_grid
                     * delta_U_Ld_L3_grid)
        integral = jnp.sum(self.jacobian * integrand)

        # Return overlap
        return 4 / jnp.pi * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L2_ms
            * cg_L3_ms * integral
        )
    

### TODO: Apply same changes as in F_1?
class A4:
    """Class for calculating the overlap matrix element:
    < \phi | \delta U^\dagger J_0 | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):

        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])
        
        # Set cos(\theta) integration mesh
        self.ntot_theta = 21
        gq = GaussQuadrature(self.ntot_theta)
        self.theta_array, self.theta_weights = gq(0, jnp.pi)
        
        # Meshgrids for integrations over \theta and k_2
        self.theta_grid, self.k2_grid = jnp.meshgrid(
            self.theta_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk2_grid = jnp.meshgrid(
            self.theta_weights, self.dwf.k_weights, indexing='ij'
        )
        
        # Jacobian for integrations over \theta and k_2
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk2_grid
                         * self.k2_grid ** 2)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over T_1
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_T1, in_axes=(0, None)
        )(self.T_1_array, args)
        
        # Sum over T_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_T1(self, T_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        and T_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # Vectorize overlap method over L_1
        args = pp, thetap, q, m_S_f, m_J_d, T_1
        overlap_array = form_factors * vmap(
            self.overlap_wrt_L1, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L1(self, L_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, and L_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, m_S_f, m_J_d, T_1 = args
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # Vectorize overlap method over J_1
        args = pp, q, m_S_f, m_J_d, T_1, L_1
        overlap_array = factor_T1_L1 * Y_L_1 * vmap(
            self.overlap_wrt_J1, in_axes=(0, None)
        )(self.J_1_array, args)
        
        # Sum over J_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_J1(self, J_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, and J_1.
        """
        
        # Unpack other arguments
        pp, q, m_S_f, m_J_d, T_1, L_1 = args
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # Vectorize overlap method over m_s
        args = pp, q, m_J_d, T_1, L_1, J_1
        overlap_array = cg_L1_mSf * vmap(
            self.overlap_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_ms(self, m_s, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, and m_s.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1 = args
        
        # Vectorize overlap method over L_2
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s
        overlap_array = vmap(
            self.overlap_wrt_L2, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_2
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L2(self, L_2, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, and L_2.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s = args
        
        # CG < J_1 m_J_d | L_2 m_J_d - m_s S = 1 m_s >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # Vectorize overlap method over L_d
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_2
        overlap_array = cg_L2_ms * vmap(
            self.overlap_wrt_Ld, in_axes=(0, None)
        )(self.L_d_array, args)
        
        # Sum over L_d
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_Ld(self, L_d, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_2, and L_d.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_2 = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # CG < L_d m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_Ld_ms = self.cg.get_coefficient(L_d, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Dot product k2_vector \dot q_vector
        k2qx_grid = self.k2_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k2_vector - q_vector / 2|
        k2q_minus_grid = jnp.sqrt(self.k2_grid ** 2 + q ** 2 / 4 - k2qx_grid)
        # Angle between the unit vector z^\hat and k2_vector - q_vector / 2
        cos_alphap_k2 = ((self.k2_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k2q_minus_grid)
        
        # Legendre polynomials
        P_L2_grid = self.sf.plm(jnp.cos(self.theta_grid), L_2, m_L)
        P_Ld_grid = self.sf.plm(cos_alphap_k2, L_d, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k2q_grid = self.dwf.compute_wf(k2q_minus_grid, L_d)
        
        # \delta U^\dagger_{J_1, L_1, L_2, J_1, S = 1, T_1}(p', k_2)
        delta_U_grid = self.deltaU(pp, self.k2_grid, J_1, L_1, L_2, 1, T_1,
                                   hc=True)

        # Integrate over \theta and k_2
        integrand = P_L2_grid * P_Ld_grid * psi_k2q_grid * delta_U_grid
        integral = jnp.sum(self.jacobian * integrand)

        # Return overlap
        return 2 * jnp.sqrt(2 / jnp.pi) * cg_Ld_ms * integral
    
    
### TODO: Apply same changes as in F_1?
class A3:
    """Class for calculating the overlap matrix element:
    < \phi | \delta U^\dagger J_0 \delta U^\dagger | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):

        # Set L_max as instance attribute
        self.L_max = L_max

        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])
        
        # Set cos(\theta) integration mesh
        self.ntot_theta = 21
        gq = GaussQuadrature(self.ntot_theta)
        self.theta_array, self.theta_weights = gq(0, jnp.pi)
        
        # Meshgrids for integrations over \theta, k_2, and k_4
        self.theta_grid, self.k2_grid, self.k4_grid = jnp.meshgrid(
            self.theta_array, self.dwf.k_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk2_grid, dk4_grid = jnp.meshgrid(
            self.theta_weights, self.dwf.k_weights, self.dwf.k_weights,
            indexing='ij'
        )
        
        # Jacobian for integrations over \theta, k_2, and k_4
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk2_grid
                         * self.k2_grid ** 2 * dk4_grid * self.k4_grid ** 2)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over T_1
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_T1, in_axes=(0, None)
        )(self.T_1_array, args)
        
        # Sum over T_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_T1(self, T_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        and T_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # Vectorize overlap method over L_1
        args = pp, thetap, q, m_S_f, m_J_d, T_1
        overlap_array = form_factors * vmap(
            self.overlap_wrt_L1, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L1(self, L_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, and L_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, m_S_f, m_J_d, T_1 = args
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # Vectorize overlap method over J_1
        args = pp, q, m_S_f, m_J_d, T_1, L_1
        overlap_array = factor_T1_L1 * Y_L_1 * vmap(
            self.overlap_wrt_J1, in_axes=(0, None)
        )(self.J_1_array, args)
        
        # Sum over J_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_J1(self, J_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, and J_1.
        """
        
        # Unpack other arguments
        pp, q, m_S_f, m_J_d, T_1, L_1 = args
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # Vectorize overlap method over m_s
        args = pp, q, m_J_d, T_1, L_1, J_1
        overlap_array = cg_L1_mSf * vmap(
            self.overlap_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_ms(self, m_s, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, and m_s.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1 = args
        
        # Vectorize overlap method over L_2
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s
        overlap_array = vmap(
            self.overlap_wrt_L2, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_2
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L2(self, L_2, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, and L_2.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s = args
        
        # CG < J_1 m_J_d | L_2 m_J_d - m_s S = 1 m_s >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # Vectorize overlap method over L_3
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_2
        overlap_array = cg_L2_ms * vmap(
            self.overlap_wrt_L3, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_3
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L3(self, L_3, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_2, and L_3.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_2 = args
        
        # CG < L_3 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Vectorize overlap method over L_d
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_2, L_3
        overlap_array = cg_L3_ms * vmap(
            self.overlap_wrt_Ld, in_axes=(0, None)
        )(self.L_d_array, args)
        
        # Sum over L_d
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_Ld(self, L_d, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_2, L_3, and L_d.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_2, L_3 = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # Dot product k2_vector \dot q_vector
        k2qx_grid = self.k2_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k2_vector - q_vector / 2|
        k2q_minus_grid = jnp.sqrt(self.k2_grid ** 2 + q ** 2 / 4 - k2qx_grid)
        # Angle between the unit vector z^\hat and k2_vector - q_vector / 2
        cos_alphap_k2 = ((self.k2_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k2q_minus_grid)
        
        # Legendre polynomials
        P_L2_grid = self.sf.plm(jnp.cos(self.theta_grid), L_2, m_L)
        P_L3_grid = self.sf.plm(cos_alphap_k2, L_3, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k4_grid = self.dwf.compute_wf(self.k4_grid, L_d)
        
        # \delta U^\dagger_{J_1, L_1, L_2, S = 1, T_1}(p', k_2)
        delta_U_L1_L2_grid = self.deltaU(pp, self.k2_grid, J_1, L_1, L_2, 1,
                                         T_1, hc=True)
        
        # \delta U_{J = 1, L_d, L_3, S = 1, T = 0}(k_4, |k_2 - q / 2|)
        delta_U_Ld_L3_grid = self.deltaU(self.k4_grid, k2q_minus_grid, 1, L_d,
                                         L_3, 1, 0)
        
        # Integrate over \theta, k_2, and k_4
        integrand = (P_L2_grid * P_L3_grid * psi_k4_grid * delta_U_L1_L2_grid
                     * delta_U_Ld_L3_grid)
        integral = jnp.sum(self.jacobian * integrand)

        # Return overlap
        return 4 / jnp.pi * jnp.sqrt(2 / jnp.pi) * integral
    
    
### TODO: Apply same changes as in F_1?
class A2:
    """Class for calculating the overlap matrix element:
    < \phi | \delta U^\dagger \delta U J_0 | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):

        # Set L_max as instance attribute
        self.L_max = L_max

        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])
        
        # Set cos(\theta) integration mesh
        self.ntot_theta = 21
        gq = GaussQuadrature(self.ntot_theta)
        self.theta_array, self.theta_weights = gq(0, jnp.pi)
        
        # Meshgrids for integrations over \theta, k_2, and k_3
        self.theta_grid, self.k2_grid, self.k3_grid = jnp.meshgrid(
            self.theta_array, self.dwf.k_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk2_grid, dk3_grid = jnp.meshgrid(
            self.theta_weights, self.dwf.k_weights, self.dwf.k_weights,
            indexing='ij'
        )
        
        # Jacobian for integrations over \theta and k_3
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk2_grid
                         * self.k2_grid ** 2 * dk3_grid * self.k3_grid ** 2)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over T_1
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_T1, in_axes=(0, None)
        )(self.T_1_array, args)
        
        # Sum over T_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_T1(self, T_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        and T_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # Vectorize overlap method over L_1
        args = pp, thetap, q, m_S_f, m_J_d, T_1
        overlap_array = form_factors * vmap(
            self.overlap_wrt_L1, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L1(self, L_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, and L_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, m_S_f, m_J_d, T_1 = args
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # Vectorize overlap method over J_1
        args = pp, q, m_S_f, m_J_d, T_1, L_1
        overlap_array = factor_T1_L1 * Y_L_1 * vmap(
            self.overlap_wrt_J1, in_axes=(0, None)
        )(self.J_1_array, args)
        
        # Sum over J_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_J1(self, J_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, and J_1.
        """
        
        # Unpack other arguments
        pp, q, m_S_f, m_J_d, T_1, L_1 = args
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # Vectorize overlap method over m_s
        args = pp, q, m_J_d, T_1, L_1, J_1
        overlap_array = cg_L1_mSf * vmap(
            self.overlap_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_ms(self, m_s, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, and m_s.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1 = args
        
        # Vectorize overlap method over L_3
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s
        overlap_array = vmap(
            self.overlap_wrt_L3, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_3
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L3(self, L_3, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, and L_3.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s = args
        
        # CG < J_1 m_J_d | L_3 m_J_d - m_s S = 1 m_s >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, J_1,
                                           m_J_d)

        # Vectorize overlap method over L_2
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3
        overlap_array = cg_L3_ms * vmap(
            self.overlap_wrt_L2, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_2
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L2(self, L_2, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, and L_2.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3 = args
        
        # Vectorize overlap method over L_d
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_2
        overlap_array = vmap(
            self.overlap_wrt_Ld, in_axes=(0, None)
        )(self.L_d_array, args)
        
        # Sum over L_d
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_Ld(self, L_d, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_2, and L_d.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_2 = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # CG < L_d m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_Ld_ms = self.cg.get_coefficient(L_d, m_L, 1, m_s, 1, m_J_d)
        
        # Dot product k3_vector \dot q_vector
        k3qx_grid = self.k3_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k3_vector - q_vector / 2|
        k3q_minus_grid = jnp.sqrt(self.k3_grid ** 2 + q ** 2 / 4 - k3qx_grid)
        # Angle between the unit vector z^\hat and k3_vector - q_vector / 2
        cos_alphap_k3 = ((self.k3_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k3q_minus_grid)
        
        # Legendre polynomials
        P_L3_grid = self.sf.plm(jnp.cos(self.theta_grid), L_3, m_L)
        P_Ld_grid = self.sf.plm(cos_alphap_k3, L_d, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k3q_grid = self.dwf.compute_wf(k3q_minus_grid, L_d)
        
        # \delta U_{J_1, L_2, L_1, J_1, S = 1, T_1}(k_2, p')
        delta_U_L2_L1_grid = self.deltaU(self.k2_grid, pp, J_1, L_2, L_1, 1,
                                         T_1)
        
        # \delta U_{J_1, L_2, L_3, S = 1, T_1}(k_2, k_3)
        delta_U_L2_L3_grid = self.deltaU(self.k2_grid, self.k3_grid, J_1, L_2,
                                         L_3, 1, T_1)
        
        # Integrate over \theta, k_2, and k_3
        integrand = (P_L3_grid * P_Ld_grid * psi_k3q_grid * delta_U_L2_L1_grid
                     * delta_U_L2_L3_grid)
        integral = jnp.sum(self.jacobian * integrand)

        # Return overlap
        return 4 / jnp.pi * jnp.sqrt(2 / jnp.pi) * cg_Ld_ms * integral
    
    
### TODO: Apply same changes as in F_1?
class A1:
    """Class for calculating the overlap matrix element:
    < \phi | \delta U^\dagger \delta U J_0 \delta U^\dagger | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):

        # Set L_max as instance attribute
        self.L_max = L_max

        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf

        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])
        
        # Set cos(\theta) integration mesh
        self.ntot_theta = 21
        self.theta_array, self.theta_weights = GaussQuadrature(self.ntot_theta)(
            0, jnp.pi)
        
        # Set integration mesh for momenta
        self.k_array, self.k_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Jacobians
        self.k_jacobian = self.k_weights * self.k_array ** 2
        self.theta_jacobian = self.theta_weights * jnp.sin(self.theta_array)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over T_1
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_T1, in_axes=(0, None)
        )(self.T_1_array, args)
        
        # Sum over T_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_T1(self, T_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        and T_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # Vectorize overlap method over L_1
        args = pp, thetap, q, m_S_f, m_J_d, T_1
        overlap_array = form_factors * vmap(
            self.overlap_wrt_L1, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L1(self, L_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, and L_1.
        """
        
        # Unpack other arguments
        pp, thetap, q, m_S_f, m_J_d, T_1 = args
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # Vectorize overlap method over J_1
        args = pp, q, m_S_f, m_J_d, T_1, L_1
        overlap_array = factor_T1_L1 * Y_L_1 * vmap(
            self.overlap_wrt_J1, in_axes=(0, None)
        )(self.J_1_array, args)
        
        # Sum over J_1
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_J1(self, J_1, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, and J_1.
        """
        
        # Unpack other arguments
        pp, q, m_S_f, m_J_d, T_1, L_1 = args
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # Vectorize overlap method over m_s
        args = pp, q, m_J_d, T_1, L_1, J_1
        overlap_array = cg_L1_mSf * vmap(
            self.overlap_wrt_ms, in_axes=(0, None)
        )(self.m_s_array, args)
        
        # Sum over m_s
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_ms(self, m_s, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, and m_s.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1 = args
        
        # Vectorize overlap method over L_3
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s
        overlap_array = vmap(
            self.overlap_wrt_L3, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_3
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L3(self, L_3, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, and L_3.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s = args
        
        # CG < J_1 m_J_d | L_3 m_J_d - m_s S = 1 m_s >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # Vectorize overlap method over L_4
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3
        overlap_array = cg_L3_ms * vmap(
            self.overlap_wrt_L4, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_4
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L4(self, L_4, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, and L_4.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3 = args
        
        # CG < L_4 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L4_ms = self.cg.get_coefficient(L_4, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Vectorize overlap method over L_2
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4
        overlap_array = cg_L4_ms * vmap(
            self.overlap_wrt_L2, in_axes=(0, None)
        )(self.L_array, args)
        
        # Sum over L_2
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_L2(self, L_2, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_4, and L_2.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4 = args
        
        # Vectorize overlap method over L_d
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2
        overlap_array = vmap(
            self.overlap_wrt_Ld, in_axes=(0, None)
        )(self.L_d_array, args)
        
        # Sum over L_d
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_Ld(self, L_d, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_4, L_2, and L_d.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2 = args
        
        # Vectorize overlap method over \theta
        args = pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d
        overlap_array = vmap(
            self.overlap_wrt_theta, in_axes=(0, None)
        )(self.theta_array, args)
        
        # Integrate over \theta
        integral_theta = jnp.sum(self.theta_jacobian * overlap_array)

        return integral_theta

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_theta(self, theta, args):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d, and \theta.
        """
        
        # Unpack other arguments
        pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d = args

        # Legendre polynomial P_{L_3 m_J_d - m_s}(cos(\theta))
        P_L3 = self.sf.plm(jnp.cos(theta), L_3, m_J_d - m_s)
        
        # Integrate over k_2
        init_val_k2 = (0.0, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2,
                       L_d, theta)
        val_k2 = fori_loop(0, self.dwf.ntot, self.overlap_wrt_k2, init_val_k2)
        integral_k2 = val_k2[0]
        
        # Integrand over \theta
        integrand_theta = P_L3 * integral_k2

        return integrand_theta

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k2(self, i, val_k2):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d, \theta, and k_2.
        """
        
        # Integration variable
        k_2 = self.k_array[i]
        jacobian_k2 = self.k_jacobian[i]
        
        # Unpack val
        (integral_k2, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d,
         theta) = val_k2
        
        # \delta U_{J_1, L_2, L_1, S = 1, T_1}(k_2, p')
        delta_U_L2_L1 = self.deltaU(k_2, pp, J_1, L_2, L_1, 1, T_1)
        
        # Integrate over k_3
        init_val_k3 = (0.0, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2,
                       L_d, theta, k_2)
        val_k3 = fori_loop(0, self.dwf.ntot, self.overlap_wrt_k3, init_val_k3)
        integral_k3 = val_k3[0]
        
        # Integrand over k_2
        integrand_k2 = delta_U_L2_L1 * integral_k3
        
        # Sum up contribution to compute integral over k_2
        integral_k2 += integrand_k2 * jacobian_k2
        
        # Re-pack val
        val_k2 = (integral_k2, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2,
                  L_d, theta)
        
        return val_k2
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k3(self, i, val_k3):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d, \theta, k_2, and k_3.
        """
        
        # Integration variable
        k_3 = self.k_array[i]
        jacobian_k3 = self.k_jacobian[i]
        
        # Unpack val
        (integral_k3, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d,
         theta, k_2) = val_k3

        # \delta U_{J_1, L_2, L_3, S = 1, T_1}(k_2, k_3)
        delta_U_L2_L3 = self.deltaU(k_2, k_3, J_1, L_2, L_3, 1, T_1)
        
        # Integrate over k_5
        init_val_k5 = (0.0, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2,
                       L_d, theta, k_2, k_3)
        val_k5 = fori_loop(0, self.dwf.ntot, self.overlap_wrt_k5, init_val_k5)
        integral_k5 = val_k5[0]
        
        # Integrand over k_3
        integrand_k3 = delta_U_L2_L3 * integral_k5

        # Sum up contribution to compute integral over k_3
        integral_k3 += integrand_k3 * jacobian_k3
        
        # Re-pack val
        val_k3 = (integral_k3, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2,
                  L_d, theta, k_2)
        
        return val_k3

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k5(self, i, val_k5):
        """Overlap matrix element for particular values of S_f, m_S_f, m_J_d,
        T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d, \theta, k_2, k_3, and k_5.
        """
        
        # Integration variable
        k_5 = self.k_array[i]
        jacobian_k5 = self.k_jacobian[i]
        
        # Unpack other arguments
        (integral_k5, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2, L_d,
         theta, k_2, k_3) = val_k5
        
        # Dot product k3_vector \dot q_vector
        k3qx = k_3 * q * jnp.cos(theta)
        # Magnitude of momenta |k3_vector - q_vector / 2|
        k3q_minus = jnp.sqrt(k_3 ** 2 + q ** 2 / 4 - k3qx)
        # Angle between the unit vector z^\hat and k3_vector - q_vector / 2
        cos_alphap_k3 = (k_3 * jnp.cos(theta) - q / 2) / k3q_minus
        
        # Legendre polynomial P_{L_4 m_J_d - m_s}(cos(\alpha'(k_3, \theta)))
        P_L4 = self.sf.plm(cos_alphap_k3, L_4, m_J_d - m_s)
        
        # Deuteron wave function [fm^3/2]
        psi_k5 = self.dwf.compute_wf(k_5, L_d)
        
        # \delta U_{J = 1, L_d, L_4, S = 1, T = 0}(k_5, |k_3 - q / 2|)
        delta_U_Ld_L4 = self.deltaU(k_5, k3q_minus, 1, L_d, L_4, 1, 0)
        
        # Integrand over k_5
        integrand_k5 = P_L4 * psi_k5 * delta_U_Ld_L4
        
        # Sum up contribution to compute integral over k_5
        integral_k5 += (8 / jnp.pi ** 2 * jnp.sqrt(2 / jnp.pi) * integrand_k5
                        * jacobian_k5)
        
        # Re-pack val
        val_k5 = (integral_k5, pp, q, m_J_d, T_1, L_1, J_1, m_s, L_3, L_4, L_2,
                  L_d, theta, k_2, k_3)
        
        return val_k5
    
    
class F4:
    """Class for calculating final state interactions of the overlap matrix
    element.
        < \phi | T^\dagger(\lambda) G_0^\dagger J_0 | \psi_i(\lambda) >
    Note, \lambda sets the SRG evolution of the initial state and final state,
    where \lambda = \infty is unevolved.
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):
        
        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf
        
        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)

        # Initialize T-matrix class
        self.tmatrix = TMatrix(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])
        
        # Set cos(\theta) integration mesh and its Jacobian
        ntot_theta = 21
        self.theta_array, self.theta_weights = GaussQuadrature(ntot_theta)(
            0, jnp.pi
        )
        self.jacobian = self.theta_weights * jnp.sin(self.theta_array)
        
        # Set integration mesh for momenta
        self.k_array, self.k_weights = self.dwf.k_array, self.dwf.k_weights
        
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for L_2 in L_array:
                        for m_s in m_s_array:
                            for L_d in L_d_array:
                                        
                                # Check if the partial wave channels are
                                # physical
                                tmatrix_bool = self.tmatrix.channel_is_physical(
                                    J_1, L_2, L_1, 1, T_1
                                )

                                # Check if T_1 and L_1 factor is 0
                                if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                    lt_bool = False
                                else:
                                    lt_bool = True
                                        
                                # Append combination to quantum numbers
                                if tmatrix_bool and lt_bool:
                                    quantum_numbers.append([T_1, L_1, J_1, L_2,
                                                            m_s, L_d])
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, L_2, m_s, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                           m_J_d)
        
        # CG < J_1 m_J_d | L_2 m_J_d - m_s S=1 m_s >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_d m_J_d - m_s S=1 m_s | J=1 m_J_d >
        cg_Ld_ms = self.cg.get_coefficient(L_d, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Call T-matrix for given partial wave channel [fm]
        T_matrix = jnp.conj(self.tmatrix.compute(pp, J_1, L_2, L_1, 1, T_1))
        
        # Half off-shell T-matrix with shape (ntot_k, 1)
        ntot_k = self.tmatrix.ntot
        T_half_offshell_array = T_matrix[:ntot_k, ntot_k]
        
        # On-shell T-matrix (scalar)
        T_onshell = T_matrix[ntot_k, ntot_k]
        
        # k_2 momentum mesh
        k2_array, k2_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Denominator of Green's function
        greens_func_array = 1 / (pp ** 2 - k2_array ** 2)
        
        # Repackage arguments with L_d
        args = q, m_J_d, T_1, L_1, J_1, L_2, m_s, L_d

        # Vectorize overlap method over k_2 and attach k_2 ** 2 factor
        f_k2_array = (
            k2_array ** 2 * greens_func_array * T_half_offshell_array
            * self.vmap_overlap_wrt_k2(k2_array, args)
        )
        # Integrate over k_2
        integral_f_k2 = jnp.sum(k2_weights * f_k2_array)
        
        # Evaluate part with k_2 = p'
        f_pp = pp ** 2 * T_onshell * self.overlap_wrt_k2(pp, args)
        Lamb = self.tmatrix.Lamb
        # \int_0^\Lamb dk_2 / (p'^2 - k_2^2)
        integral_k2 = jnp.sum(k2_weights * greens_func_array)
        
        # Factor of 2 for < J_0 > = 2 < J_0^- > (Eq. 20 in More 2015)
        overlap = 2 * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L2_ms 
            * cg_Ld_ms * (
                integral_f_k2 - f_pp * (
                    integral_k2 - 1 / (2 * pp) * (
                        jnp.log((Lamb + pp) / (Lamb - pp)) + 1j * jnp.pi
                    )
                )
            )
        )
        
        return overlap
    
    @partial(jit, static_argnums=(0,))
    def vmap_overlap_wrt_k2(self, k2_array, args):
        """vmap of the method below w.r.t. k_2."""
        
        return vmap(self.overlap_wrt_k2, in_axes=(0, None))(k2_array, args)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k2(self, k_2, args):
        """Overlap matrix element for particular values k_2."""
        
        # Unpack other arguments
        q, m_J_d, T_1, L_1, J_1, L_2, m_s, L_d = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # Dot product k2_vector \dot q_vector
        k2qx_array = k_2 * q * jnp.cos(self.theta_array)
        # Magnitude of momenta |k2_vector - q_vector / 2|
        k2q_minus_array = jnp.sqrt(k_2 ** 2 + q ** 2 / 4 - k2qx_array)
        # Angle between the unit vector z^\hat and k2_vector - q_vector / 2
        cos_alphap_k2 = ((k_2 * jnp.cos(self.theta_array) - q / 2)
                         / k2q_minus_array)

        # Legendre polynomials
        P_Ld_array = self.sf.plm(cos_alphap_k2, L_d, m_L)
        P_L2_array = self.sf.plm(jnp.cos(self.theta_array), L_2, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k2_array = self.dwf.compute_wf(k2q_minus_array, L_d)

        # Integrand over \theta
        integrand = P_Ld_array * P_L2_array * psi_k2_array
        
        # Return integrand over k_2 by integrating over \theta
        return jnp.sum(self.jacobian * integrand)
    
  
class F3:
    """Class for calculating the overlap matrix element:
    < \phi |
        T^\dagger(\lambda) G_0^\dagger J_0 \delta U^\dagger 
    | \psi_i(\lambda) >
    """

    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):
        
        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf
        
        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)

        # Initialize T-matrix class
        self.tmatrix = TMatrix(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])

        # Set cos(\theta) integration mesh
        ntot_theta = 21
        theta_array, theta_weights = GaussQuadrature(ntot_theta)(0, jnp.pi)
        
        # Set integration mesh for momenta
        self.k_array, self.k_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Meshgrids and Jacobian for \theta and k_5
        self.theta_grid, self.k5_grid = jnp.meshgrid(
            theta_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk5_grid = jnp.meshgrid(
            theta_weights, self.dwf.k_weights, indexing='ij'
        )
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk5_grid
                         * self.k5_grid ** 2)
        
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for L_2 in L_array:
                        for m_s in m_s_array:
                            for L_3 in L_array:
                                for L_d in L_d_array:
                                        
                                    # Check if the partial wave channels are
                                    # physical
                                    tmatrix_bool = (
                                        self.tmatrix.channel_is_physical(
                                            J_1, L_2, L_1, 1, T_1
                                        )
                                    )
                                    deltaU_bool = (
                                        self.deltaU.channel_is_physical(
                                            1, L_d, L_3, 1, 0
                                        )
                                    )
                                    channel_bool = tmatrix_bool and deltaU_bool
                                        
                                    # Check if T_1 and L_1 factor is 0
                                    if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                        lt_bool = False
                                    else:
                                        lt_bool = True
                                        
                                    # Append combination to quantum numbers
                                    if channel_bool and lt_bool:
                                        quantum_numbers.append(
                                            [T_1, L_1, J_1, L_2, m_s, L_3, L_d]
                                        )
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, L_2, m_s, L_3, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # CG < J_1 m_J_d | L_2 m_J_d - m_s S = 1 m_s >
        cg_L2_ms = self.cg.get_coefficient(L_2, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_3 m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Call T-matrix for given partial wave channel [fm]
        T_matrix = jnp.conj(self.tmatrix.compute(pp, J_1, L_2, L_1, 1, T_1))
        
        # Half off-shell T-matrix with shape (ntot_k, 1)
        ntot_k = self.tmatrix.ntot
        T_half_offshell_array = T_matrix[:ntot_k, ntot_k]
        
        # On-shell T-matrix (scalar)
        T_onshell = T_matrix[ntot_k, ntot_k]
        
        # k_2 momentum mesh
        k2_array, k2_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Denominator of Green's function
        greens_func_array = 1 / (pp ** 2 - k2_array ** 2)
        
        # Repackage arguments with L_d
        args = q, m_J_d, T_1, L_1, J_1, L_2, m_s, L_3, L_d

        # Vectorize overlap method over k_2 and attach k_2 ** 2 factor
        f_k2_array = (
            k2_array ** 2 * greens_func_array * T_half_offshell_array
            * self.vmap_overlap_wrt_k2(k2_array, args)
        )
        # Integrate over k_2
        integral_f_k2 = jnp.sum(k2_weights * f_k2_array)
        
        # Evaluate part with k_2 = p'
        f_pp = pp ** 2 * T_onshell * self.overlap_wrt_k2(pp, args)
        Lamb = self.tmatrix.Lamb
        # \int_0^\Lamb dk_2 / (p'^2 - k_2^2)
        integral_k2 = jnp.sum(k2_weights * greens_func_array)
        
        # Overlap matrix element
        overlap = 4 / jnp.pi * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L2_ms 
            * cg_L3_ms * (
                integral_f_k2 - f_pp * (
                    integral_k2 - 1 / (2 * pp) * (
                        jnp.log((Lamb + pp) / (Lamb - pp)) + 1j * jnp.pi
                    )
                )
            )
        )
        
        return overlap
    
    @partial(jit, static_argnums=(0,))
    def vmap_overlap_wrt_k2(self, k2_array, args):
        """vmap of the method below w.r.t. k_2."""
        
        return vmap(self.overlap_wrt_k2, in_axes=(0, None))(k2_array, args)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k2(self, k_2, args):
        """Overlap matrix element for particular values k_2."""
        
        # Unpack other arguments
        q, m_J_d, T_1, L_1, J_1, L_2, m_s, L_3, L_d = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # Dot product k2_vector \dot q_vector
        k2qx_grid = k_2 * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k2_vector - q_vector / 2|
        k2q_minus_grid = jnp.sqrt(k_2 ** 2 + q ** 2 / 4 - k2qx_grid)
        # Angle between the unit vector z^\hat and k2_vector - q_vector / 2
        cos_alphap_k2 = ((k_2 * jnp.cos(self.theta_grid) - q / 2)
                         / k2q_minus_grid)

        # Legendre polynomials
        P_L2_grid = self.sf.plm(jnp.cos(self.theta_grid), L_2, m_L)
        P_L3_grid = self.sf.plm(cos_alphap_k2, L_3, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k5_grid = self.dwf.compute_wf(self.k5_grid, L_d)
        
        # \delta U_{J = 1, L_d, L_3, S = 1, T = 0}(k_5, |k_2 - q / 2|)
        delta_U_Ld_L3_grid = self.deltaU(self.k5_grid, k2q_minus_grid, 1, L_d,
                                         L_3, 1, 0)

        # Integrand over k_2, \theta, and k_4
        integrand = P_L2_grid * P_L3_grid * psi_k5_grid * delta_U_Ld_L3_grid
        
        # Return integrand over k_2 by integrating over \theta and k_5
        return jnp.sum(self.jacobian * integrand)

class F2:
    """Class for calculating the overlap matrix element:
    < \phi | T^\dagger(\lambda) G_0^\dagger \delta U J_0 | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):
        
        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf
        
        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)

        # Initialize T-matrix class
        self.tmatrix = TMatrix(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)
        
        # Arrays for sums over quantum numbers
        self.T_1_array = jnp.array([0, 1])
        self.L_array = jnp.arange(0, L_max + 1, 1)
        self.m_s_array = jnp.array([-1, 0, 1])
        self.J_1_array = jnp.arange(0, L_max + 2, 1)
        self.L_d_array = jnp.array([0, 2])
        
        # Set cos(\theta) integration mesh
        ntot_theta = 21
        theta_array, theta_weights = GaussQuadrature(ntot_theta)(0, jnp.pi)
        
        # Set integration mesh for momenta
        self.k_array, self.k_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Meshgrids and Jacobian for \theta and k_4
        self.theta_grid, self.k4_grid = jnp.meshgrid(
            theta_array, self.dwf.k_array, indexing='ij'
        )
        dtheta_grid, dk4_grid = jnp.meshgrid(
            theta_weights, self.dwf.k_weights, indexing='ij'
        )
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk4_grid
                         * self.k4_grid ** 2)
        
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for L_2 in L_array:
                        for L_3 in L_array:
                            for m_s in m_s_array:
                                for L_d in L_d_array:
                                        
                                    # Check if the partial wave channels are
                                    # physical
                                    tmatrix_bool = (
                                        self.tmatrix.channel_is_physical(
                                            J_1, L_2, L_1, 1, T_1
                                        )
                                    )
                                    deltaU_bool = (
                                        self.deltaU.channel_is_physical(
                                            J_1, L_2, L_3, 1, T_1
                                        )
                                    )
                                    channel_bool = tmatrix_bool and deltaU_bool
                                        
                                    # Check if T_1 and L_1 factor is 0
                                    if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                        lt_bool = False
                                    else:
                                        lt_bool = True
                                        
                                    # Append combination to quantum numbers
                                    if channel_bool and lt_bool:
                                        quantum_numbers.append(
                                            [T_1, L_1, J_1, L_2, L_3, m_s, L_d]
                                        )
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, L_2, L_3, m_s, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # CG < J_1 m_J_d | L_3 m_J_d - m_s S = 1 m_s >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_d m_J_d - m_s S = 1 m_s | J = 1 m_J_d >
        cg_Ld_ms = self.cg.get_coefficient(L_d, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Call T-matrix for given partial wave channel [fm]
        T_matrix = jnp.conj(self.tmatrix.compute(pp, J_1, L_2, L_1, 1, T_1))

        # Half off-shell T-matrix with shape (ntot_k, 1)
        ntot_k = self.tmatrix.ntot
        T_half_offshell_array = T_matrix[:ntot_k, ntot_k]
        
        # On-shell T-matrix (scalar)
        T_onshell = T_matrix[ntot_k, ntot_k]
        
        # k_2 momentum mesh
        k2_array, k2_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Denominator of Green's function
        greens_func_array = 1 / (pp ** 2 - k2_array ** 2)
        
        # Repackage arguments with L_d
        args = q, m_J_d, T_1, L_1, J_1, L_2, L_3, m_s, L_d

        # Vectorize overlap method over k_2 and attach k_2 ** 2 factor
        f_k2_array = (
            k2_array ** 2 * greens_func_array * T_half_offshell_array
            * self.vmap_overlap_wrt_k2(k2_array, args)
        )
        # Integrate over k_2
        integral_f_k2 = jnp.sum(k2_weights * f_k2_array)
        
        # Evaluate part with k_2 = p'
        f_pp = pp ** 2 * T_onshell * self.overlap_wrt_k2(pp, args)
        Lamb = self.tmatrix.Lamb
        # \int_0^\Lamb dk_2 / (p'^2 - k_2^2)
        integral_k2 = jnp.sum(k2_weights * greens_func_array)
        
        # Overlap matrix element
        overlap = 4 / jnp.pi * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L3_ms 
            * cg_Ld_ms * (
                integral_f_k2 - f_pp * (
                    integral_k2 - 1 / (2 * pp) * (
                        jnp.log((Lamb + pp) / (Lamb - pp)) + 1j * jnp.pi
                    )
                )
            )
        )
        
        return overlap
    
    @partial(jit, static_argnums=(0,))
    def vmap_overlap_wrt_k2(self, k2_array, args):
        """vmap of the method below w.r.t. k_2."""
        
        return vmap(self.overlap_wrt_k2, in_axes=(0, None))(k2_array, args)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k2(self, k_2, args):
        """Overlap matrix element for particular values k_2."""
        
        # Unpack other arguments
        q, m_J_d, T_1, L_1, J_1, L_2, L_3, m_s, L_d = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # Dot product k4_vector \dot q_vector
        k4qx_grid = self.k4_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k4_vector - q_vector / 2|
        k4q_minus_grid = jnp.sqrt(self.k4_grid ** 2 + q ** 2 / 4 - k4qx_grid)
        # Angle between the unit vector z^\hat and k4_vector - q_vector / 2
        cos_alphap_k4 = ((self.k4_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k4q_minus_grid)

        # Legendre polynomials
        P_L3_grid = self.sf.plm(jnp.cos(self.theta_grid), L_3, m_L)
        P_Ld_grid = self.sf.plm(cos_alphap_k4, L_d, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k4q_grid = self.dwf.compute_wf(k4q_minus_grid, L_d)
        
        # \delta U_{J_1, L_2, L_3, S = 1, T_1}(k_2, k_4)
        delta_U_L2_L3_grid = self.deltaU(k_2, self.k4_grid, J_1, L_2, L_3, 1,
                                         T_1)

        # Integrand over k_2, \theta, and k_4
        integrand = P_L3_grid * P_Ld_grid * psi_k4q_grid * delta_U_L2_L3_grid
        
        # Return integrand over k_2 by integrating over \theta and k_4
        return jnp.sum(self.jacobian * integrand)
    
    
class F1:
    """Class for calculating the overlap matrix element:
    < \phi |
        T^\dagger(\lambda) G_0^\dagger \delta U J_0 \delta U^\dagger
    | \psi_i(\lambda) >
    """
    
    def __init__(self, kvnn, kmax, kmid, ntot, lamb, L_max, sf, cg):
        
        # Set L_max as instance attribute
        self.L_max = L_max
        
        # Callable for getting special functions (e.g., spherical harmonics)
        self.sf = sf
        
        # Callable for getting Clebsch-Gordan coefficients
        self.cg = cg
        
        # Set-up deuteron wave function
        self.dwf = DeuteronWaveFunction(kvnn, kmax, kmid, ntot, lamb)

        # Initialize T-matrix class
        self.tmatrix = TMatrix(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Initialize \delta U class
        self.deltaU = DeltaU(kvnn, kmax, kmid, ntot, lamb, L_max)
        
        # Repackage nested sum over quantum numbers into one JAX array
        self.quantum_numbers = self.get_quantum_numbers(L_max)

        # Set cos(\theta) integration mesh
        ntot_theta = 21
        theta_array, theta_weights = GaussQuadrature(ntot_theta)(0, jnp.pi)
        
        # Set integration mesh for momenta
        self.k_array, self.k_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Meshgrids and Jacobian for \theta, k_4, and k_6 integrations
        self.theta_grid, self.k4_grid, self.k6_grid = jnp.meshgrid(
            theta_array, self.k_array, self.k_array, indexing='ij'
        )
        dtheta_grid, dk4_grid, dk6_grid = jnp.meshgrid(
            theta_weights, self.k_weights, self.k_weights, indexing='ij'
        )
        self.jacobian = (dtheta_grid * jnp.sin(self.theta_grid) * dk4_grid
                         * self.k4_grid ** 2 * dk6_grid * self.k6_grid ** 2)
    
    def get_quantum_numbers(self, L_max):
        """Repackage all quantum numbers into one big JAX array."""
        
        # Arrays for sums over quantum numbers
        T_1_array = jnp.array([0, 1])
        L_array = jnp.arange(0, L_max + 1, 1)
        m_s_array = jnp.array([-1, 0, 1])
        L_d_array = jnp.array([0, 2])
        
        quantum_numbers = []
        for T_1 in T_1_array:
            for L_1 in L_array:
                J_1_array = jnp.arange(jnp.abs(L_1 - 1), L_1 + 2, 1)
                for J_1 in J_1_array:
                    for L_2 in L_array:
                        for L_3 in L_array:
                            for m_s in m_s_array:
                                for L_4 in L_array:
                                    for L_d in L_d_array:
                                        
                                        # Check if the partial wave channels
                                        # are physical
                                        tmatrix_bool = (
                                            self.tmatrix.channel_is_physical(
                                                J_1, L_2, L_1, 1, T_1
                                            )
                                        )
                                        deltaU_23_bool = (
                                            self.deltaU.channel_is_physical(
                                                J_1, L_2, L_3, 1, T_1
                                            )
                                        )
                                        deltaU_d4_bool = (
                                            self.deltaU.channel_is_physical(
                                                1, L_d, L_4, 1, 0
                                            )
                                        )
                                        channel_bool = (
                                            tmatrix_bool and deltaU_23_bool
                                            and deltaU_d4_bool
                                        )
                                        
                                        # Check if T_1 and L_1 factor is 0
                                        if 1 + (-1) ** T_1 * (-1) ** L_1 == 0:
                                            lt_bool = False
                                        else:
                                            lt_bool = True
                                        
                                        # Append combination to quantum numbers
                                        if channel_bool and lt_bool:
                                            quantum_numbers.append(
                                                [T_1, L_1, J_1, L_2, L_3, m_s,
                                                 L_4, L_d]
                                            )
                                        
        # Return quantum numbers as JAX array
        return jnp.asarray(quantum_numbers)
                                        
    @partial(jit, static_argnums=(0,))
    def __call__(self, pp, thetap, q, gep, gen, fL_quantum_numbers):
        """Overlap matrix element [fm^3/2] for particular values of S_f, m_S_f,
        and m_J_d.
        """
        
        # Unpack f_L quantum numbers (no dependence on S_f)
        m_S_f, m_J_d = fL_quantum_numbers
        
        # Vectorize overlap method over quantum numbers
        args = pp, thetap, q, gep, gen, m_S_f, m_J_d
        overlap_array = vmap(
            self.overlap_wrt_quantum_numbers, in_axes=(0, None)
        )(self.quantum_numbers, args)
        
        # Sum over all quantum numbers
        return jnp.sum(overlap_array)
    
    @partial(jit, static_argnums=(0,))
    def overlap_wrt_quantum_numbers(self, quantum_numbers, args):
        """Overlap matrix element for particular quantum numbers."""
        
        # Unpack quantum numbers
        T_1, L_1, J_1, L_2, L_3, m_s, L_4, L_d = quantum_numbers
        
        # Unpack other arguments
        pp, thetap, q, gep, gen, m_S_f, m_J_d = args
        
        # Form factor part
        form_factors = gep + (-1) ** T_1 * gen
        
        # T_1 and L_1 factor
        factor_T1_L1 = 1 + (-1) ** T_1 * (-1) ** L_1
        
        # Spherical harmonic w.r.t. \theta'
        Y_L_1 = self.sf.ylm(thetap, 0.0, L_1, m_J_d - m_S_f)
        
        # CG < L_1 m_J_d - m_S_f S = 1 m_S_f | J_1 m_J_d >
        cg_L1_mSf = self.cg.get_coefficient(L_1, m_J_d - m_S_f, 1, m_S_f, J_1,
                                            m_J_d)
        
        # CG < J_1 m_J_d | L_3 m_J_d - m_s S = 1 m_s >
        cg_L3_ms = self.cg.get_coefficient(L_3, m_J_d - m_s, 1, m_s, J_1, m_J_d)
        
        # CG < L_4 m_J_d - m_s S = 1 m_s | 1 m_J_d >
        cg_L4_ms = self.cg.get_coefficient(L_4, m_J_d - m_s, 1, m_s, 1, m_J_d)
        
        # Call T-matrix for given partial wave channel [fm]
        T_matrix = jnp.conj(self.tmatrix.compute(pp, J_1, L_2, L_1, 1, T_1))

        # Half off-shell T-matrix with shape (ntot_k, 1)
        ntot_k = self.tmatrix.ntot
        T_half_offshell_array = T_matrix[:ntot_k, ntot_k]
        
        # On-shell T-matrix (scalar)
        T_onshell = T_matrix[ntot_k, ntot_k]
        
        # k_2 momentum mesh
        k2_array, k2_weights = self.dwf.k_array, self.dwf.k_weights
        
        # Denominator of Green's function
        greens_func_array = 1 / (pp ** 2 - k2_array ** 2)
        
        # Repackage arguments with L_d
        args = q, m_J_d, T_1, J_1, L_2, L_3, m_s, L_4, L_d

        # Vectorize overlap method over k_2 and attach k_2 ** 2 factor
        f_k2_array = (
            k2_array ** 2 * greens_func_array * T_half_offshell_array
            * self.vmap_overlap_wrt_k2(k2_array, args)
        )
        # Integrate over k_2
        integral_f_k2 = jnp.sum(k2_weights * f_k2_array)
        
        # Evaluate part with k_2 = p'
        f_pp = pp ** 2 * T_onshell * self.overlap_wrt_k2(pp, args)
        Lamb = self.tmatrix.Lamb
        # \int_0^\Lamb dk_2 / (p'^2 - k_2^2)
        integral_k2 = jnp.sum(k2_weights * greens_func_array)
        
        # Overlap matrix element
        overlap = 8 / jnp.pi ** 2 * jnp.sqrt(2 / jnp.pi) * (
            form_factors * factor_T1_L1 * Y_L_1 * cg_L1_mSf * cg_L3_ms 
            * cg_L4_ms * (
                integral_f_k2 - f_pp * (
                    integral_k2 - 1 / (2 * pp) * (
                        jnp.log((Lamb + pp) / (Lamb - pp)) + 1j * jnp.pi
                    )
                )
            )
        )
        
        return overlap

    @partial(jit, static_argnums=(0,))
    def vmap_overlap_wrt_k2(self, k2_array, args):
        """vmap of the method below w.r.t. k_2."""
        
        return vmap(self.overlap_wrt_k2, in_axes=(0, None))(k2_array, args)

    @partial(jit, static_argnums=(0,))
    def overlap_wrt_k2(self, k_2, args):
        """Overlap matrix element for particular values k_2."""
        
        # Unpack other arguments
        q, m_J_d, T_1, J_1, L_2, L_3, m_s, L_4, L_d = args
        
        # Fix orbital angular momentum projection
        m_L = m_J_d - m_s
        
        # Dot product k4_vector \dot q_vector
        k4qx_grid = self.k4_grid * q * jnp.cos(self.theta_grid)
        # Magnitude of momenta |k4_vector - q_vector / 2|
        k4q_minus_grid = jnp.sqrt(self.k4_grid ** 2 + q ** 2 / 4 - k4qx_grid)
        # Angle between the unit vector z^\hat and k4_vector - q_vector / 2
        cos_alphap_k4 = ((self.k4_grid * jnp.cos(self.theta_grid) - q / 2)
                         / k4q_minus_grid)

        # Legendre polynomials
        P_L3_grid = self.sf.plm(jnp.cos(self.theta_grid), L_3, m_L)
        P_L4_grid = self.sf.plm(cos_alphap_k4, L_4, m_L)
        
        # Deuteron wave function [fm^3/2]
        psi_k6_grid = self.dwf.compute_wf(self.k6_grid, L_d)
        
        # \delta U_{J_1, L_2, L_3, S = 1, T_1}(k_2, k_4)
        delta_U_L2_L3_grid = self.deltaU(k_2, self.k4_grid, J_1, L_2, L_3, 1,
                                         T_1)
        
        # \delta U_{J = 1, L_d, L_4, S = 1, T = 0}(k_6, |k_4 - q / 2|)
        delta_U_Ld_L4_grid = self.deltaU(self.k6_grid, k4q_minus_grid, 1, L_d,
                                         L_4, 1, 0)

        # Integrand over k_2, \theta, k_4, and k_6
        integrand = (P_L3_grid * P_L4_grid * psi_k6_grid * delta_U_L2_L3_grid
                     * delta_U_Ld_L4_grid)
        
        # Return integrand over k_2 by integrating over \theta, k_4, and k_6
        return jnp.sum(self.jacobian * integrand)