#!/usr/bin/env python3

"""
File: test_momentum_distribution_script.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. This particular version implements
several JAX speed-ups.

Last update: May 23, 2023

"""

# Python imports
from jax import config, jit, vmap
import jax.numpy as jnp
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.special import sph_harm
import time
import vegas

# Imports from scripts
from scripts.integration import momentum_mesh, unattach_weights_from_matrix
from scripts.potentials import Potential
from scripts.srg import load_srg_transformation
from scripts.tools import replace_periods

from clebsch_gordan import (
    compute_clebsch_gordan_table, compute_clebsch_gordan_array,
    clebsch_gordan_coefficient, clebsch_gordan_coefficient_vmap
)
from momentum_distribution_quantum_numbers import (
    PartialWaveChannel, get_delta_U_quantum_numbers,
    get_delta_U2_quantum_numbers, quantum_number_array
)
from single_particle_states import WoodsSaxon


# Enable double-precision with JAX arrays
config.update("jax_enable_x64", True)
    

@vegas.batchintegrand
class DeltaUIntegrand:
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""
    
    
    def __init__(self, q, tau, cg_array, N_j, delta_U_quantum_numbers,
                 phi_functions, delta_U_functions, delta_U_dagger_functions):
        
        # Set instance attributes
        self.q = q
        self.tau = tau
        self.cg_array = cg_array
        self.N_j = N_j
        self.delta_U_quantum_numbers = delta_U_quantum_numbers
        self.delU_Ntot = len(delta_U_quantum_numbers)
        self.phi_functions = phi_functions
        self.delta_U_functions = delta_U_functions
        self.delta_U_dagger_functions = delta_U_dagger_functions
        
        # Vectorize functions
        self.vectorized_phi = np.vectorize(get_phi_function,
                                           excluded=['phi_functions'])
        self.vectorized_get_delta_U = np.vectorize(self.get_delta_U)
        self.vectorized_get_delta_U_dag = np.vectorize(self.get_delta_U_dag)
        

    def __call__(self, x_array):
        """Evaluate the integrand at several points simultaneously."""
        
        # Relative momenta k
        k = x_array[:,0]
        theta_k = x_array[:,1]
        phi_k = x_array[:,2]
        k_vector = build_vector(k, theta_k, phi_k)
            
        # C.o.M. momenta K
        K = x_array[:,3]
        theta_K = x_array[:,4]
        phi_K = x_array[:,5]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Choose z-axis to be along q_vector
        q_vector = np.zeros_like(k_vector)
        q_vector[-1, :] = self.q
        q = np.repeat(self.q, k.size)
        theta_q = np.zeros_like(q)
        phi_q = np.zeros_like(q)
        
        # Calculate vector q-K/2
        qK_vector = q_vector - K_vector/2
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
        
        # Calculate vector K/2+k
        k1_vector = K_vector/2 + k_vector
        k1, theta_k1, phi_k1 = get_vector_components(k1_vector)
        
        # Calculate vector K/2-k
        k2_vector = K_vector/2 - k_vector
        k2, theta_k2, phi_k2 = get_vector_components(k2_vector)
        
        # Calculate vector K-q
        Kq_vector = K_vector-q_vector
        Kq, theta_Kq, phi_Kq = get_vector_components(Kq_vector)
            
        # Calculate the Jacobian determinant
        jacobian = k**2 * np.sin(theta_k) * K**2 * np.sin(theta_K)
        
        # Samples of quantum number sets
        quantum_number_array = self.get_quantum_numbers(x_array[:,6])

        # Unpack quantum numbers
        sigma_1 = quantum_number_array[:, 0]
        sigma_2 = quantum_number_array[:, 1]
        sigma = quantum_number_array[:, 2]
        sigmap = quantum_number_array[:, 3]
        tau_1 = quantum_number_array[:, 4]
        tau_2 = quantum_number_array[:, 5]
        taup = quantum_number_array[:, 6]
        
        s = np.repeat(1/2, sigma_1.size)
        t = s
        tau_array = np.repeat(self.tau, taup.size)
        
        n_alpha = quantum_number_array[:, 7]
        l_alpha = quantum_number_array[:, 8]
        j_alpha = quantum_number_array[:, 9]
        m_j_alpha = quantum_number_array[:, 10]
        m_t_alpha = quantum_number_array[:, 11]
        
        n_beta = quantum_number_array[:, 12]
        l_beta = quantum_number_array[:, 13]
        j_beta = quantum_number_array[:, 14]
        m_j_beta = quantum_number_array[:, 15]
        m_t_beta = quantum_number_array[:, 16]
        
        S = quantum_number_array[:, 17]
        M_S = quantum_number_array[:, 18]
        M_Sp = quantum_number_array[:, 19]                                  
        L = quantum_number_array[:, 20]
        M_L = quantum_number_array[:, 21]
        Lp = quantum_number_array[:, 22]
        M_Lp = quantum_number_array[:, 23]
        J = quantum_number_array[:, 24]
        M_J = quantum_number_array[:, 25]                             
        T = quantum_number_array[:, 26]
        M_T = quantum_number_array[:, 27]
        
        # < \sigma_1 \sigma_2 | S M_S >
        spin_12_cg = clebsch_gordan_coefficient_vmap(
            s, sigma_1, s, sigma_2, S, M_S, self.cg_array, self.N_j)
            
        # < S M_S' | \sigma \sigma' >
        spin_ssp_cg = clebsch_gordan_coefficient_vmap(
            s, sigma, s, sigmap, S, M_Sp, self.cg_array, self.N_j)
            
        # < \tau_1 \tau_2 | T M_T >
        isospin_12_cg = clebsch_gordan_coefficient_vmap(
            t, tau_1, t, tau_2, T, M_T, self.cg_array, self.N_j)
            
        # < T M_T | \tau \tau' >
        isospin_ttp_cg = clebsch_gordan_coefficient_vmap(
            t, tau_array, t, taup, T, M_T, self.cg_array, self.N_j)
            
        # < L M_L S M_S | J M_J >
        lsj_cg = clebsch_gordan_coefficient_vmap(L, M_L, S, M_S, J, M_J,
                                                 self.cg_array, self.N_j)
            
        # < J M_J | L' M_L' S M_S' >
        lpsj_cg = clebsch_gordan_coefficient_vmap(Lp, M_Lp, S, M_Sp, J, M_J,
                                                  self.cg_array, self.N_j)
            
        # 1 - (-1)^(L+S+T) factor
        # lst_factor = 1 - (-1) ** (L+S+T)
        lst_factor = 2
            
        # 1 - (-1)^(L'+S+T) factor
        # lpst_factor = 1 - (-1) ** (Lp+S+T)
        lpst_factor = 2
            
        # Spherical harmonics
        Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
        Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
            
        # < k (L S) J T | \delta U | |q-K/2| (L' S) J T >
        delta_U_partial_wave = self.vectorized_get_delta_U(k, qK, L, Lp, J, S,
                                                           T)
            
        # < |q-K/2| (L' S) J T | \delta U^\dagger | k (L S) J T >
        delta_U_dag_partial_wave = self.vectorized_get_delta_U_dag(qK, k, Lp, L,
                                                                   J, S, T)

        # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
        psi_alpha_1 = psi(
            n_alpha, l_alpha, s, j_alpha, m_j_alpha, m_t_alpha,
            k1, theta_k1, phi_k1, sigma_1, tau_1,
            self.cg_array, self.N_j, self.vectorized_phi, self.phi_functions
        )
            
        # \psi_\beta(K/2-k; \sigma_2, \tau_2)
        psi_beta_2 = psi(
            n_beta, l_beta, s, j_beta, m_j_beta, m_t_beta,
            k2, theta_k2, phi_k2, sigma_2, tau_2,
            self.cg_array, self.N_j, self.vectorized_phi, self.phi_functions
        )
            
        # \psi_\alpha(q; \sigma, \tau)
        psi_alpha_q = psi(
            n_alpha, l_alpha, s, j_alpha, m_j_alpha, m_t_alpha,
            q, theta_q, phi_q, sigma, tau_array,
            self.cg_array, self.N_j, self.vectorized_phi, self.phi_functions
        )
            
        # \psi_\beta(K-q; \sigma', \tau')
        psi_beta_Kq = psi(
            n_beta, l_beta, s, j_beta, m_j_beta, m_t_beta,
            Kq, theta_Kq, phi_Kq, sigmap, taup,
            self.cg_array, self.N_j, self.vectorized_phi, self.phi_functions
        )
            
        # \psi_\beta(q; \sigma, \tau)
        psi_beta_q = psi(
            n_beta, l_beta, s, j_beta, m_j_beta, m_t_beta,
            q, theta_q, phi_q, sigma, tau_array,
            self.cg_array, self.N_j, self.vectorized_phi, self.phi_functions
        )
            
        # \psi_\alpha(K-q; \sigma', \tau')
        psi_alpha_Kq = psi(
            n_alpha, l_alpha, s, j_alpha, m_j_alpha, m_t_alpha,
            Kq, theta_Kq, phi_Kq, sigmap, taup,
            self.cg_array, self.N_j, self.vectorized_phi, self.phi_functions
        )
            
        # \delta U term
        delta_U_term = (
            spin_12_cg * spin_ssp_cg * isospin_12_cg * isospin_ttp_cg
            * lsj_cg * lpsj_cg * lst_factor * lpst_factor
            * Y_L_k * np.conj(Y_Lp_qK) * delta_U_partial_wave
            * np.conj(psi_alpha_1) * np.conj(psi_beta_2) * (
                psi_beta_Kq * psi_alpha_q - psi_alpha_Kq * psi_beta_q
            )
        )
            
        # \delta U^\dagger term
        delta_U_dag_term = (
            spin_ssp_cg * spin_12_cg * isospin_ttp_cg * isospin_12_cg
            * lpsj_cg * lsj_cg * lst_factor * lpst_factor
            * Y_Lp_qK * np.conj(Y_L_k) * delta_U_dag_partial_wave
            * psi_alpha_1 * psi_beta_2 * (
                np.conj(psi_beta_Kq) * np.conj(psi_alpha_q)
                - np.conj(psi_alpha_Kq) * np.conj(psi_beta_q)
            )
        )
            
        # Add together for full integrand
        integrand = 1/2 * 1/2 * 2/np.pi * jacobian * self.delU_Ntot * (
            delta_U_term + delta_U_dag_term
        )

        return integrand.real
    
    
    def get_quantum_numbers(self, x):
        
        index = np.floor(x * self.delU_Ntot).astype(int)
        return self.delta_U_quantum_numbers[index]
    
    
    def get_delta_U(self, k, kp, L, Lp, J, S, T):
        """Gets the partial wave matrix element of \delta U."""
        
        key = (L, Lp, J, S, T)
        return self.delta_U_functions[key].ev(k, kp)
    
    
    def get_delta_U_dag(self, k, kp, L, Lp, J, S, T):
        """Gets the partial wave matrix element of \delta U^\dagger."""
        
        key = (L, Lp, J, S, T)
        return self.delta_U_dagger_functions[key].ev(k, kp)
    
    
# STILL NEED \DELTA U^2 STUFF
    
    
def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return jnp.array(x == y, dtype=int)


@jit
def kronecker_delta_vmap(x, y):
    return vmap(kronecker_delta)(x, y)


@jit   
def build_vector(k, theta, phi):
    """
    Build a vector from input spherical coordinates.

    Parameters
    ----------
    k : array_like
        Magnitude of the vector.
    theta : array_like
        Polar angle of the vector in the range [0, \pi].
    phi : array_like
        Azimuthal angle of the vector in the range [0, 2\pi].

    Returns
    -------
    k_vector : array_like
        Output vector with shape (3,1) or (3,N) where N is the size of each
        input.

    """

    k_vector = jnp.array([k * jnp.sin(theta) * jnp.cos(phi),
                          k * jnp.sin(theta) * jnp.sin(phi),
                          k * jnp.cos(theta)], dtype=jnp.float64)

    return k_vector


@jit
def get_vector_components(k_vector):
    """
    Get the spherical coordinates from an input vector.

    Parameters
    ----------
    k_vector : array_like
        Output vector with shape (3,1) or (3,N) where N is the size of each
        input.

    Returns
    -------
    k : array_like
        Magnitude of the vector.
    theta : array_like
        Polar angle of the vector in the range [0, \pi].
    phi : array_like
        Azimuthal angle of the vector in the range [0, 2\pi].

    """

    k = jnp.linalg.norm(k_vector, axis=0)
    theta = jnp.arccos(k_vector[2]/k)
    phi = jnp.arctan2(k_vector[1], k_vector[0])

    return k, theta, phi


def get_sp_wave_functions(woods_saxon, kmax, kmid, ntot):
    """Set interpolating functions for s.p. wave functions \phi."""
    
    occ_states = woods_saxon.occ_states

    phi_functions = {}
    for sp_state in occ_states: 
        key = (sp_state.n, sp_state.l, sp_state.j, sp_state.m_t)
        phi_functions[key] = woods_saxon.get_wf_kspace(sp_state, kmax, kmid,
                                                       ntot, interpolate=True)
            
    return phi_functions


def get_phi_function(n, l, j, m_t, k, phi_functions):
    
    key = (n, l, j, m_t)

    return phi_functions[key](k)


def psi(n, l, s, j, m_j, m_t, k, theta, phi, sigma, tau, cg_array, N_j,
        phi_vect, phi_functions):
    """Single-particle wave function."""
    
    # Calculate \phi_\alpha(q)
    phi_sp_wf = phi_vect(n, l, j, m_t, k, phi_functions)
    
    # Calculate spinor spherical harmonic
    Y_jml = spinor_sph_harm(theta, phi, l, s, j, m_j, sigma, cg_array, N_j)
    
    # Isospinor indexed by \tau \chi_{m_t}(\tau)
    chi_tau = kronecker_delta_vmap(tau, m_t).block_until_ready()

    return phi_sp_wf * Y_jml * chi_tau


def spinor_sph_harm(theta, phi, l, s, j, m_j, sigma, cg_array, N_j):
    
    m_s = sigma
    
    m_l = m_j - m_s
    
    # Calls clebsch_gordan_coefficient_vmap
    cg = clebsch_gordan_coefficient_vmap(l, m_l, s, m_s, j, m_j, cg_array,
                                         N_j).block_until_ready()
    
    # Calls sph_harm (which is already vectorized)
    Y_lm = np.where(np.abs(m_l) <= l, sph_harm(m_l, l, phi, theta), 0+0j)

    return cg * Y_lm


def interpolate_delta_U(channel, potential, generator, lamb,
                        hermitian_conjugate=False):
    """Interpolate \delta U(k, k') for the given channel."""

    # Get momentum mesh
    kmax, kmid, ntot = potential.kmax, potential.kmid, potential.ntot
    k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
    
    U_matrix_weights = load_srg_transformation(potential, generator, lamb)

    # Calculate \delta U = U - I
    I_matrix_weights = np.eye(len(U_matrix_weights))
    if hermitian_conjugate:
        delU_matrix_weights = (U_matrix_weights - I_matrix_weights).T
    else:
        delU_matrix_weights = U_matrix_weights - I_matrix_weights

    # Get specific sub-block if coupled-channel
    if channel in ['3S1-3D1', '3P2-3F2', '3D3-3G3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[:ntot,ntot:]
        )
    elif channel in ['3D1-3S1', '3F2-3P2', '3G3-3D3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[ntot:,:ntot]
        )
    elif channel in ['3D1-3D1', '3F2-3F2', '3G3-3G3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[ntot:,ntot:]
        )
    else:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[:ntot,:ntot]
        )
        
    # Interpolate \delta U(k, k')
    delU_func = RectBivariateSpline(k_array, k_array, delU_matrix)

    return delU_func


def get_delta_U_functions(channels, kvnn, kmax, kmid, ntot, generator, lamb):
    """Get \delta U and \delta U^\dagger functions."""

    delta_U_functions = {}
    delta_U_dagger_functions = {}
    for channel in channels:
        
        # Set channel argument to be compatible with potential functions
        if channel[:3] == '3D1':
            channel_arg = '3S1'
        elif channel[:3] == '3F2':
            channel_arg = '3P2'
        elif channel[:3] == '3G3':
            channel_arg == '3D3'
        else:
            channel_arg = channel[:3]
            
        potential = Potential(kvnn, channel_arg, kmax, kmid, ntot)
        
        pwc = PartialWaveChannel(channel)
        key = (pwc.L, pwc.Lp, pwc.J, pwc.S, pwc.T)
        
        delta_U_functions[key] = interpolate_delta_U(channel, potential,
                                                     generator, lamb)
        delta_U_dagger_functions[key] = interpolate_delta_U(
            channel, potential, generator, lamb, hermitian_conjugate=True
        )
        
    return delta_U_functions, delta_U_dagger_functions

    
def compute_I_term(q_array, tau, occ_states, cg_array, N_j, phi_functions):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in occ_states:

            # Single-particle wave function with z-axis along q_vector
            key = (alpha.n, alpha.l, alpha.j, alpha.m_t)
            phi_sp_wf = phi_functions[key](q_array)
            
            # Calculate spinor spherical harmonic
            m_l = alpha.m_j - sigma
            if np.abs(m_l) <= alpha.l:
                cg = clebsch_gordan_coefficient(
                    alpha.l, m_l, 1/2, sigma, alpha.j, alpha.m_j, cg_array, N_j
                )
                Y_lm = sph_harm(m_l, alpha.l, 0, 0)
                Y_jml = cg * Y_lm
            else:
                Y_jml = 0+0j
            
            # Isospinor indexed by \tau \chi_{m_t}(\tau)
            chi_tau = kronecker_delta(tau, alpha.m_t)

            psi_alpha_array = phi_sp_wf * Y_jml * chi_tau

            I_array += np.abs(psi_alpha_array) ** 2
                    
    return I_array


def compute_delta_U_term(
        q_array, tau, occupied_states, cg_array, N_j, phi_functions,
        delta_U_functions, delta_U_dagger_functions):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    if tau == 1/2:
        nucleon = 'proton'
    elif tau == -1/2:
        nucleon = 'neutron'
    file_name = f"{nucleus_name}_{nucleon}_delta_U_quantum_numbers.txt"
    
    # Try loading the file first
    try:
        
        delta_U_quantum_numbers = np.loadtxt('quantum_numbers/' + file_name)
        # delta_U_quantum_numbers = jnp.array(
        #     np.loadtxt('quantum_numbers/' + file_name)
        # )
    
    # Find all possible combinations and save file
    except OSError:
        
        print("Starting \delta U quantum numbers...")
        cg_table = compute_clebsch_gordan_table(4)
        quantum_numbers = get_delta_U_quantum_numbers(tau, occupied_states,
                                                      channels, cg_table)
        delta_U_quantum_numbers = quantum_number_array(quantum_numbers,
                                                       file_name)
        print("Finished with \delta U quantum numbers.\n")
        
    # Relative momenta from 0 to maximum of momentum mesh
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    # Polar angle
    theta_limits = [0, np.pi]
    # Azimuthal angle
    phi_limits = [0, 2*np.pi]
    # Sum over quantum numbers using vegas integration trick
    sum_limits = [0, 1]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits,
                              sum_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = DeltaUIntegrand(
            q, tau, cg_array, N_j, delta_U_quantum_numbers, phi_functions,
            delta_U_functions, delta_U_dagger_functions
        )

        # Train the integrator
        integ(integrand, nitn=5, neval=5e4)
        # Final result
        result = integ(integrand, nitn=10, neval=5e4)

        delta_U_array[i] = result.mean
        delta_U_errors[i] = result.sdev

        t1 = time.time()

        percent = (i+1)/len(q_array)*100
        mins = (t1-t0)/60
        print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                    
    return delta_U_array, delta_U_errors


def compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, kmax, kmid, ntot, lamb,
        generator='Wegner', print_normalization=False, ipm_only=False, save=True
):
    """Compute the single-nucleon momentum distribution."""
    
    # Compute table of Clebsch-Gordan coefficients
    cg_array, N_j = compute_clebsch_gordan_array(4)
    print("Done calculating Clebsch-Gordan table.\n")
    
    # Set-up single-particle states
    woods_saxon = WoodsSaxon(nucleus_name, Z, N, run_woodsaxon=False)
    phi_functions = get_sp_wave_functions(woods_saxon, 10.0, 2.0, 120)
    
    # Set-up \delta U and \delta U^\dagger functions
    delta_U_functions, delta_U_dagger_functions = get_delta_U_functions(
        channels, kvnn, kmax, kmid, ntot, generator, lamb
    )
    
    # Set momentum mesh
    q_array, q_weights = momentum_mesh(10.0, 4.0, 100, nmod=70)
    
    # Compute the I term
    I_array = compute_I_term(q_array, tau, woods_saxon.occ_states, cg_array,
                             N_j, phi_functions)
    
    # Independent particle model has no \delta U contributions
    if ipm_only:
        
        delta_U_array = np.zeros_like(I_array)
        delta_U_errors = np.zeros_like(I_array)
        delta_U2_array = np.zeros_like(I_array)
        delta_U2_errors = np.zeros_like(I_array)
    
    # Include \delta U, \delta U^\dagger, and \delta U \delta U^\dagger
    else:

        t0 = time.time()

        # Compute \delta U + \delta U^\dagger term using vegas
        print("Beginning \delta U linear terms.\n")
        t1 = time.time()
        delta_U_array, delta_U_errors = compute_delta_U_term(
            q_array, tau, woods_saxon.occ_states, cg_array, N_j, phi_functions,
            delta_U_functions, delta_U_dagger_functions
        )
        t2 = time.time()
        print(f"Done after {(t2-t1)/60:.3f} minutes.\n")
        
        # TESTING
        delta_U2_array = np.zeros_like(I_array)
        delta_U2_errors = np.zeros_like(I_array)
        
        t5 = time.time()
        print(f"Total time elapsed: {(t5-t0)/3600:.3f} hours.\n")
    
    # Combine each term for the total momentum distribution [fm^3]
    n_array = I_array + delta_U_array + delta_U2_array
    n_errors = np.sqrt(delta_U_errors**2 + delta_U2_errors**2)

    if print_normalization:
        normalization = compute_normalization(q_array, q_weights, n_array)
        print(f"Normalization = {normalization:.5f}.")
        
    if save and not(ipm_only):  # Do not save IPM-only data
        save_momentum_distribution(
            nucleus_name, tau, kvnn, lamb, q_array, q_weights, n_array,
            n_errors, I_array, delta_U_array, delta_U_errors, delta_U2_array,
            delta_U2_errors
        )
    
    return q_array, q_weights, n_array, n_errors


def compute_normalization(q_array, q_weights, n_array):
    """Compute the normalization of the momentum distribution."""

    return 4*np.pi * np.sum(q_weights * q_array**2 * n_array)


def save_momentum_distribution(
        nucleus_name, tau, kvnn, lamb, q_array, q_weights, n_array, n_errors,
        I_array, delta_U_array, delta_U_errors, delta_U2_array, delta_U2_errors
):
    """Save the momentum distribution along with the isolated contributions."""
    
    data = np.vstack(
        (q_array, q_weights, n_array, n_errors, I_array, delta_U_array,
         delta_U_errors, delta_U2_array, delta_U2_errors)
    ).T
            
    if tau == 1/2:
        nucleon = 'proton'
    elif tau == -1/2:
        nucleon = 'neutron'
                
    hdr = ("q, q weight, n(q), n(q) error, I, \delta U + \delta U^\dagger,"
           " \delta U + \delta U^\dagger error, \delta U^2, \delta U^2 error\n")
    
    directory = 'momentum_distributions/'

    file_name = replace_periods(f"{nucleus_name}_{nucleon}_momentum"
                                f"_distribution_kvnn_{kvnn}_lamb_{lamb}")
    
    np.savetxt(directory + file_name + '.txt', data, header=hdr)

    
def load_momentum_distribution(nucleus_name, nucleon, kvnn, lamb):
    """Load and return the momentum distribution along with the isolated
    contributions.
    """
    
    directory = 'momentum_distributions/'

    file_name = replace_periods(f"{nucleus_name}_{nucleon}_momentum"
                                f"_distribution_kvnn_{kvnn}_lamb_{lamb}")
    
    data = np.loadtxt(directory + file_name + '.txt')
    
    q_array = data[:, 0]
    q_weights = data[:, 1]
    n_array = data[:, 2]
    n_errors = data[:, 3]
    I_array = data[:, 4]
    delta_U_array = data[:, 5]
    delta_U_errors = data[:, 6]
    delta_U2_array = data[:, 7]
    delta_U2_errors = data[:, 8]
    
    return (q_array, q_weights, n_array, n_errors, I_array, delta_U_array,
            delta_U_errors, delta_U2_array, delta_U2_errors)


if __name__ == '__main__':
    
    # Nucleus
    # nucleus_name, Z, N = 'He4', 2, 2
    nucleus_name, Z, N = 'O16', 8, 8
    # nucleus_name, Z, N = 'Ca40', 20, 20
    # nucleus_name, Z, N = 'Ca48', 20, 28
    
    # Nucleon
    tau = 1/2
    
    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    
    # NN potential and momentum mesh
    # kvnn, kmax, kmid, ntot = 6, 30.0, 4.0, 120  # AV18
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120
    # kvnn, kmax, kmid, ntot = 111, 15.0, 3.0, 120  # SMS N4LO 450 MeV
    
    # SRG \lambda value
    # lamb = 1.35
    lamb = 1.5
    # lamb = 2.0
    # lamb = 3.0
    # lamb = 6.0

    # Compute and save the momentum distribution
    q_array, q_weights, n_array, n_errors = compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, kmax, kmid, ntot, lamb,
        print_normalization=True, save=False
    )