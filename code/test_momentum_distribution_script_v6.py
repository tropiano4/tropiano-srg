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
batch mode in vegas and directly sums over all non-zero combinations of partial
wave channels and single-particle quantum numbers.

Last update: June 1, 2023

"""

# Python imports
import numpy as np
import numpy.linalg as la
from scipy.interpolate import RectBivariateSpline
from scipy.special import sph_harm
import time
import vegas

# Imports from scripts
from scripts.integration import momentum_mesh, unattach_weights_from_matrix
from scripts.potentials import Potential
from scripts.srg import load_srg_transformation
from scripts.tools import replace_periods

from clebsch_gordan import compute_clebsch_gordan_table
from clebsch_gordan import clebsch_gordan_coefficient_v1 as cg_func
from momentum_distribution_quantum_numbers import (
    PartialWaveChannel, get_delta_U_quantum_numbers,
    get_delta_U2_quantum_numbers, quantum_number_array
)
from single_particle_states import WoodsSaxon
        
        
def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return int(x == y)

    
def build_vector(k, theta, phi):
    """
    Build a vector from input spherical coordinates.

    Parameters
    ----------
    k : float
        Magnitude of the vector.
    theta : float
        Polar angle of the vector in the range [0, \pi].
    phi : float
        Azimuthal angle of the vector in the range [0, 2\pi].

    Returns
    -------
    k_vector : 1-D ndarray
        Output vector with shape (3,1).

    """

    k_vector = np.array([k * np.sin(theta) * np.cos(phi),
                         k * np.sin(theta) * np.sin(phi),
                         k * np.cos(theta)])

    return k_vector


def get_vector_components(k_vector):
    """
    Get the spherical coordinates from an input vector.

    Parameters
    ----------
    k_vector : 1-D ndarray
        Input vector with shape (3,1).

    Returns
    -------
    k : float
        Magnitude of the vector.
    theta : float
        Polar angle of the vector in the range [0, \pi].
    phi : float
        Azimuthal angle of the vector in the range [0, 2\pi].

    """

    k = la.norm(k_vector, axis=0)
    theta = np.arccos(k_vector[2]/k)
    phi = np.arctan2(k_vector[1], k_vector[0])

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


def psi(n, l, j, m_j, m_t, k, theta, phi, sigma, tau, cg_table, phi_functions):
    """Single-particle wave function."""

    # Calculate \phi_\alpha(k)
    if l % 2 == 0:  # Even l
        phi_sp_wf = phi_functions[(n, l, j, m_t)](k)
    else:  # Odd l needs factor of i^-l
        phi_sp_wf = 1j ** (-l) * phi_functions[(n, l, j, m_t)](k)
    
    # Calculate spinor spherical harmonic
    Y_jml = spinor_spherical_harmonic(l, j, m_j, theta, phi, sigma, cg_table)
    
    # Isospinor indexed by \tau \chi_{m_t}(\tau)
    chi_tau = kronecker_delta(tau, m_t)

    return phi_sp_wf * Y_jml * chi_tau


def spinor_spherical_harmonic(l, j, m_j, theta, phi, sigma, cg_table):
    """Spinor spherical harmonic for a s.p. state described by the quantum
    numbers j, m_j, l, and s=1/2.
    """
    
    Y_jml = 0+0j

    # Spinor indexed by \sigma \eta_{m_s}^(\sigma) = \delta_{m_s, \sigma}
    m_s = sigma
    
    # m_l must be fixed since m_j and m_s are determined
    m_l = m_j - m_s
        
    # Check that |m_l| <= l
    if np.abs(m_l) <= l:
        
        # Clebsch-Gordan coefficient for l-s coupling
        cg = cg_table[(l, m_l, 1/2, m_s, j, m_j)]
        
        # Spherical harmonic
        Y_lm = sph_harm(m_l, l, phi, theta)
            
        Y_jml += cg * Y_lm
            
    return Y_jml


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


def compute_I_term(q_array, tau, occ_states, cg_table, phi_functions):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in occ_states:

            # Loop over q
            psi_alpha_array = np.zeros_like(q_array, dtype='complex')
            for i, q in enumerate(q_array):

                # Single-particle wave function with z-axis along q_vector
                psi_alpha_array[i] = psi(
                    alpha.n, alpha.l, alpha.j, alpha.m_j, alpha.m_t, q, 0, 0,
                    sigma, tau, cg_table, phi_functions
                )
                    
            I_array += abs(psi_alpha_array)**2
                    
    return I_array


@vegas.batchintegrand
class delta_U_term_integrand:
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""
    
    def __init__(
            self, q, tau, delta_U_quantum_numbers, cg_table, phi_functions,
            delta_U_functions, delta_U_dagger_functions
    ):
        self.q = q
        self.tau = tau
        self.delta_U_quantum_numbers = delta_U_quantum_numbers
        self.delU_Ntot = delta_U_quantum_numbers.shape[0]
        self.cg_table = cg_table
        self.phi_functions = phi_functions
        self.delta_U_functions = delta_U_functions
        self.delta_U_dagger_functions = delta_U_dagger_functions

    def __call__(self, momenta_array):
        """Evaluate the integrand at several points simultaneously."""
        
        k = momenta_array[:,0]
        theta_k = momenta_array[:,1]
        phi_k = momenta_array[:,2]
        k_vector = build_vector(k, theta_k, phi_k)
            
        # C.o.M. momenta K
        K = momenta_array[:,3]
        theta_K = momenta_array[:,4]
        phi_K = momenta_array[:,5]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Choose z-axis to be along q_vector
        q_vector = np.zeros_like(k_vector)
        q_vector[-1, :] = self.q
        
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
        
        integrand = np.zeros_like(jacobian, dtype='complex')
        for i in range(self.delU_Ntot):
            
            # Unpack quantum numbers
            sigma_1 = self.delta_U_quantum_numbers[i, 0]
            sigma_2 = self.delta_U_quantum_numbers[i, 1]
            sigma = self.delta_U_quantum_numbers[i, 2]
            sigmap = self.delta_U_quantum_numbers[i, 3]
            tau_1 = self.delta_U_quantum_numbers[i, 4]
            tau_2 = self.delta_U_quantum_numbers[i, 5]
            taup = self.delta_U_quantum_numbers[i, 6]
            
            n_alpha = self.delta_U_quantum_numbers[i, 7]
            l_alpha = self.delta_U_quantum_numbers[i, 8]
            j_alpha = self.delta_U_quantum_numbers[i, 9]
            m_j_alpha = self.delta_U_quantum_numbers[i, 10]
            m_t_alpha = self.delta_U_quantum_numbers[i, 11]
            
            n_beta = self.delta_U_quantum_numbers[i, 12]
            l_beta = self.delta_U_quantum_numbers[i, 13]
            j_beta = self.delta_U_quantum_numbers[i, 14]
            m_j_beta = self.delta_U_quantum_numbers[i, 15]
            m_t_beta = self.delta_U_quantum_numbers[i, 16]
            
            S = self.delta_U_quantum_numbers[i, 17]
            M_S = self.delta_U_quantum_numbers[i, 18]
            M_Sp = self.delta_U_quantum_numbers[i, 19]                                  
            L = self.delta_U_quantum_numbers[i, 20]
            M_L = self.delta_U_quantum_numbers[i, 21]
            Lp = self.delta_U_quantum_numbers[i, 22]
            M_Lp = self.delta_U_quantum_numbers[i, 23]
            J = self.delta_U_quantum_numbers[i, 24]
            M_J = self.delta_U_quantum_numbers[i, 25]                             
            T = self.delta_U_quantum_numbers[i, 26]
            M_T = self.delta_U_quantum_numbers[i, 27]
            
            # < \sigma_1 \sigma_2 | S M_S >
            spin_12_cg = cg_func(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                  self.cg_table)
            
            # < S M_S' | \sigma \sigma' >
            spin_ssp_cg = cg_func(1/2, sigma, 1/2, sigmap, S, M_Sp,
                                  self.cg_table)
            
            # < \tau_1 \tau_2 | T M_T >
            isospin_12_cg = cg_func(1/2, tau_1, 1/2, tau_2, T, M_T,
                                    self.cg_table)
            
            # < T M_T | \tau \tau' >
            isospin_ttp_cg = cg_func(1/2, self.tau, 1/2, taup, T, M_T,
                                      self.cg_table)
            
            # < L M_L S M_S | J M_J >
            lsj_cg = cg_func(L, M_L, S, M_S, J, M_J, self.cg_table)
            
            # < J M_J | L' M_L' S M_S' >
            lpsj_cg = cg_func(Lp, M_Lp, S, M_Sp, J, M_J, self.cg_table)
            
            # 1 - (-1)^(L+S+T) factor
            lst_factor = 1 - (-1) ** (L+S+T)
            
            # 1 - (-1)^(L'+S+T) factor
            lpst_factor = 1 - (-1) ** (Lp+S+T)
            
            # Spherical harmonics
            Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
            Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
            Y_L_qK = sph_harm(M_L, L, phi_qK, theta_qK)
            Y_Lp_k = sph_harm(M_Lp, Lp, phi_k, theta_k)
            
            # < k (L S) J T | \delta U | |q-K/2| (L' S) J T >
            delta_U_partial_wave = (
                self.delta_U_functions[(L, Lp, J, S, T)].ev(k, qK))
    
            # < |q-K/2| (L' S) J T | \delta U^\dagger | k (L S) J T >
            delta_U_dag_partial_wave = (
                self.delta_U_dagger_functions[(Lp, L, J, S, T)].ev(qK, k))
            
            # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
            psi_alpha_1 = psi(
                n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
                k1, theta_k1, phi_k1, sigma_1, tau_1,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\beta(K/2-k; \sigma_2, \tau_2)
            psi_beta_2 = psi(
                n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
                k2, theta_k2, phi_k2, sigma_2, tau_2,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\alpha(q; \sigma, \tau)
            psi_alpha_q = psi(
                n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
                self.q, 0, 0, sigma, self.tau,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\beta(K-q; \sigma', \tau')
            psi_beta_Kq = psi(
                n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
                Kq, theta_Kq, phi_Kq, sigmap, taup,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\beta(q; \sigma, \tau)
            psi_beta_q = psi(
                n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
                self.q, 0, 0, sigma, self.tau,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\alpha(K-q; \sigma', \tau')
            psi_alpha_Kq = psi(
                n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
                Kq, theta_Kq, phi_Kq, sigmap, taup,
                self.cg_table, self.phi_functions
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
                * Y_L_qK * np.conj(Y_Lp_k) * delta_U_dag_partial_wave
                * psi_alpha_1 * psi_beta_2 * (
                    np.conj(psi_beta_Kq) * np.conj(psi_alpha_q)
                    - np.conj(psi_alpha_Kq) * np.conj(psi_beta_q)
                )
            )

            # Add together for full integrand
            integrand += 1/2 * 1/2 * 2/np.pi * jacobian * (
                delta_U_term + delta_U_dag_term
            )

        return integrand.real
    

def compute_delta_U_term(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions
):
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
    
    # Find all possible combinations and save file
    except OSError:
        
        print("Starting \delta U quantum numbers...")
        cg_table = compute_clebsch_gordan_table(4)
        quantum_numbers = get_delta_U_quantum_numbers(tau, occ_states, channels,
                                                      cg_table)
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

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = delta_U_term_integrand(
            q, tau, delta_U_quantum_numbers, cg_table, phi_functions,
            delta_U_functions, delta_U_dagger_functions
        )

        # Train the integrator
        integ(integrand, nitn=5, neval=1e3)
        # Final result
        result = integ(integrand, nitn=10, neval=1e3)

        delta_U_array[i] = result.mean
        delta_U_errors[i] = result.sdev

        t1 = time.time()

        percent = (i+1)/len(q_array)*100
        mins = (t1-t0)/60
        print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                    
    return delta_U_array, delta_U_errors


@vegas.batchintegrand
class delta_U2_term_integrand:
    """Evaluate the integrand of the \delta U \delta U^\dagger term."""
    
    def __init__(
            self, q, tau, delta_U2_quantum_numbers, cg_table, phi_functions,
            delta_U_functions, delta_U_dagger_functions
    ):
        self.q = q
        self.tau = tau
        self.delta_U2_quantum_numbers = delta_U2_quantum_numbers
        self.delU2_Ntot = delta_U2_quantum_numbers.shape[0]
        self.cg_table = cg_table
        self.phi_functions = phi_functions
        self.delta_U_functions = delta_U_functions
        self.delta_U_dagger_functions = delta_U_dagger_functions
        
    def __call__(self, momenta_array):
        """Evaluate the integrand at several points simultaneously."""
        
        # Relative momenta k
        k = momenta_array[:,0]
        theta_k = momenta_array[:,1]
        phi_k = momenta_array[:,2]
        k_vector = build_vector(k, theta_k, phi_k)
            
        # Relative momenta k'
        kp = momenta_array[:,3]
        theta_kp = momenta_array[:,4]
        phi_kp = momenta_array[:,5]
        kp_vector = build_vector(kp, theta_kp, phi_kp)
            
        # C.o.M. momenta K
        K = momenta_array[:,6]
        theta_K = momenta_array[:,7]
        phi_K = momenta_array[:,8]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Choose z-axis to be along q_vector
        q_vector = np.zeros_like(k_vector)
        q_vector[-1, :] = self.q
        
        # Calculate vector q-K/2
        qK_vector = q_vector - K_vector/2
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
        
        # Calculate vector K/2+k
        k1_vector = K_vector/2 + k_vector
        k1, theta_k1, phi_k1 = get_vector_components(k1_vector)
        
        # Calculate vector K/2-k
        k2_vector = K_vector/2 - k_vector
        k2, theta_k2, phi_k2 = get_vector_components(k2_vector)
        
        # Calculate vector K/2+k'
        k3_vector = K_vector/2 + kp_vector
        k3, theta_k3, phi_k3 = get_vector_components(k3_vector)
        
        # Calculate vector K/2-k'
        k4_vector = K_vector/2 - kp_vector
        k4, theta_k4, phi_k4 = get_vector_components(k4_vector)
            
        # Calculate the Jacobian determinant
        jacobian = (k**2 * np.sin(theta_k) * kp**2 * np.sin(theta_kp) * K**2
                    * np.sin(theta_K))
        
        integrand = np.zeros_like(jacobian, dtype='complex')
        for i in range(self.delU2_Ntot):
            
            # Unpack dictionary
            sigma_1 = self.delta_U2_quantum_numbers[i, 0]
            sigma_2 = self.delta_U2_quantum_numbers[i, 1]
            sigma_3 = self.delta_U2_quantum_numbers[i, 2]
            sigma_4 = self.delta_U2_quantum_numbers[i, 3]
            tau_1 = self.delta_U2_quantum_numbers[i, 4]
            tau_2 = self.delta_U2_quantum_numbers[i, 5]
            tau_3 = self.delta_U2_quantum_numbers[i, 6]
            tau_4 = self.delta_U2_quantum_numbers[i, 7]
            taup = self.delta_U2_quantum_numbers[i, 8]
            
            n_alpha = self.delta_U2_quantum_numbers[i, 9]
            l_alpha = self.delta_U2_quantum_numbers[i, 10]
            j_alpha = self.delta_U2_quantum_numbers[i, 11]
            m_j_alpha = self.delta_U2_quantum_numbers[i, 12]
            m_t_alpha = self.delta_U2_quantum_numbers[i, 13]
            
            n_beta = self.delta_U2_quantum_numbers[i, 14]
            l_beta = self.delta_U2_quantum_numbers[i, 15]
            j_beta = self.delta_U2_quantum_numbers[i, 16]
            m_j_beta = self.delta_U2_quantum_numbers[i, 17]
            m_t_beta = self.delta_U2_quantum_numbers[i, 18]
            
            S = self.delta_U2_quantum_numbers[i, 19]
            M_S = self.delta_U2_quantum_numbers[i, 20]
            M_Sp = self.delta_U2_quantum_numbers[i, 21]                                  
            L = self.delta_U2_quantum_numbers[i, 22]
            M_L = self.delta_U2_quantum_numbers[i, 23]
            Lp = self.delta_U2_quantum_numbers[i, 24]
            M_Lp = self.delta_U2_quantum_numbers[i, 25]
            J = self.delta_U2_quantum_numbers[i, 26]
            M_J = self.delta_U2_quantum_numbers[i, 27]                             
            T = self.delta_U2_quantum_numbers[i, 28]
            M_T = self.delta_U2_quantum_numbers[i, 29]
            
            M_Sppp = self.delta_U2_quantum_numbers[i, 30]
            Lpp = self.delta_U2_quantum_numbers[i, 31]
            M_Lpp = self.delta_U2_quantum_numbers[i, 32]
            Lppp = self.delta_U2_quantum_numbers[i, 33]
            M_Lppp = self.delta_U2_quantum_numbers[i, 34]
            Jp = self.delta_U2_quantum_numbers[i, 35]
            M_Jp = self.delta_U2_quantum_numbers[i, 36]
            Tp = self.delta_U2_quantum_numbers[i, 37]
            M_Tp = self.delta_U2_quantum_numbers[i, 38]
            
            # < \sigma_1 \sigma_2 | S M_S >
            spin_12_cg = cg_func(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                 self.cg_table)
            
            # < S M_S' | \sigma \sigma' >
            spin_34_cg = cg_func(1/2, sigma_3, 1/2, sigma_4, S, M_Sppp,
                                 self.cg_table)
            
            # < \tau_1 \tau_2 | T M_T >
            isospin_12_cg = cg_func(1/2, tau_1, 1/2, tau_2, T, M_T,
                                    self.cg_table)
            
            # < T M_T | \tau \tau' >
            isospin_ttp_cg = cg_func(1/2, self.tau, 1/2, taup, T, M_T,
                                     self.cg_table)
            
            # < \tau \tau' | T' M_T' >
            isospin_ttpp_cg = cg_func(1/2, self.tau, 1/2, taup, Tp, M_Tp,
                                      self.cg_table)
            
            # < T' M_T' | \tau__3 \tau_4 >
            isospin_34_cg = cg_func(1/2, tau_3, 1/2, tau_4, Tp, M_Tp,
                                    self.cg_table)
            
            # < L M_L S M_S | J M_J >
            lsj_cg = cg_func(L, M_L, S, M_S, J, M_J, self.cg_table)
            
            # < J M_J | L' M_L' S M_S' >
            lpsj_cg = cg_func(Lp, M_Lp, S, M_Sp, J, M_J, self.cg_table)
            
            # < L'' M_L'' S M_S' | J' M_J' >
            lppsjp_cg = cg_func(Lpp, M_Lpp, S, M_Sp, Jp, M_Jp, self.cg_table)
            
            # < J' M_J' | L''' M_L''' S M_S''' >
            lpppsjp_cg = cg_func(Lppp, M_Lppp, S, M_Sppp, Jp, M_Jp,
                                 self.cg_table)
            
            # 1 - (-1)^(L+S+T) factor
            lst_factor = 1 - (-1) ** (L+S+T)
            
            # 1 - (-1)^(L'+S+T) factor
            lpst_factor = 1 - (-1) ** (Lp+S+T)
            
            # 1 - (-1)^(L''+S+T') factor
            lppstp_factor = 1 - (-1) ** (Lpp+S+Tp)
            
            # 1 - (-1)^(L'''+S+T') factor
            lpppstp_factor = 1 - (-1) ** (Lppp+S+Tp)
            
            # Spherical harmonics
            Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
            Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
            Y_Lpp_qK = sph_harm(M_Lpp, Lpp, phi_qK, theta_qK)
            Y_Lppp_kp = sph_harm(M_Lppp, Lppp, phi_kp, theta_kp)
            
            # < k (L S) J T | \delta U | |q-K/2| (L' S) J T >
            delta_U_partial_wave = (
                self.delta_U_functions[(L, Lp, J, S, T)].ev(k, qK))
            
            # < |q-K/2| (L'' S) J' T' | \delta U^\dagger | k' (L''' S) J' T' >
            delta_U_dag_partial_wave = (
                self.delta_U_dagger_functions[(Lpp, Lppp, Jp, S, Tp)].ev(qK, kp)
            )
            
            # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
            psi_alpha_1 = psi(
                n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
                k1, theta_k1, phi_k1, sigma_1, tau_1,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\beta(K/2-k; \sigma_2, \tau_2)
            psi_beta_2 = psi(
                n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
                k2, theta_k2, phi_k2, sigma_2, tau_2,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\alpha(K/2+k'; \sigma_3, \tau_3)
            psi_alpha_3 = psi(
                n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
                k3, theta_k3, phi_k3, sigma_3, tau_3,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\beta(K/2-k'; \sigma_4, \tau_4)
            psi_beta_4 = psi(
                n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
                k4, theta_k4, phi_k4, sigma_4, tau_4,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\beta(K/2+k'; \sigma_3, \tau_3)
            psi_beta_3 = psi(
                n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
                k3, theta_k3, phi_k3, sigma_3, tau_3,
                self.cg_table, self.phi_functions
            )
            
            # \psi_\alpha(K/2-k'; \sigma_4, \tau_4)
            psi_alpha_4 = psi(
                n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
                k4, theta_k4, phi_k4, sigma_4, tau_4,
                self.cg_table, self.phi_functions
            )
        
            integrand += 1/4 * (1/2 * 2/np.pi) ** 2 * jacobian * (
                spin_12_cg * spin_34_cg
                * isospin_12_cg * isospin_ttp_cg
                * isospin_ttpp_cg * isospin_34_cg
                * lsj_cg * lpsj_cg * lppsjp_cg * lpppsjp_cg
                * lst_factor * lpst_factor * lppstp_factor * lpppstp_factor
                * Y_L_k * np.conj(Y_Lp_qK) * Y_Lpp_qK * np.conj(Y_Lppp_kp)
                * delta_U_partial_wave * delta_U_dag_partial_wave
                * np.conj(psi_alpha_1) * np.conj(psi_beta_2) * (
                    psi_alpha_3 * psi_beta_4 - psi_beta_3 * psi_alpha_4
                )
            )

        return integrand.real


def compute_delta_U2_term(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions
):
    """Compute the \delta U * n(q) * \delta U^\dagger term."""
        
    delta_U2_array = np.zeros_like(q_array)
    delta_U2_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    if tau == 1/2:
        nucleon = 'proton'
    elif tau == -1/2:
        nucleon = 'neutron'
    file_name = f"{nucleus_name}_{nucleon}_delta_U2_quantum_numbers.txt"
    
    # Try loading the file first
    try:
        
        delta_U2_quantum_numbers = np.loadtxt('quantum_numbers/' + file_name)
    
    # Find all possible combinations and save file
    except OSError:
        
        print("Starting \delta U \delta U^\dagger quantum numbers...")
        cg_table = compute_clebsch_gordan_table(4)
        quantum_numbers = get_delta_U2_quantum_numbers(tau, occ_states,
                                                       channels, cg_table)
        delta_U2_quantum_numbers = quantum_number_array(quantum_numbers,
                                                        file_name)
        print("Finished with \delta U \delta U^\dagger quantum numbers.\n")
        
    # Relative momenta from 0 to maximum of momentum mesh
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    
    # Polar angle
    theta_limits = [0, np.pi]
    # Azimuthal angle
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = delta_U2_term_integrand(
            q, tau, delta_U2_quantum_numbers, cg_table, phi_functions,
            delta_U_functions, delta_U_dagger_functions
        )
        
        # Train the integrator
        integ(integrand, nitn=5, neval=1e3)
        # Final result
        result = integ(integrand, nitn=10, neval=1e3)

        delta_U2_array[i] = result.mean
        delta_U2_errors[i] = result.sdev

        t1 = time.time()

        percent = (i+1)/len(q_array)*100
        mins = (t1-t0)/60
        print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                    
    return delta_U2_array, delta_U2_errors


def compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, kmax, kmid, ntot, lamb,
        generator='Wegner', print_normalization=False, ipm_only=False, save=True
):
    """Compute the single-nucleon momentum distribution."""
    
    # Compute table of Clebsch-Gordan coefficients
    cg_table = compute_clebsch_gordan_table(4)
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
    I_array = compute_I_term(q_array, tau, woods_saxon.occ_states, cg_table,
                             phi_functions)
    
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
            q_array, tau, woods_saxon.occ_states, cg_table, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
        t2 = time.time()
        print(f"Done after {(t2-t1)/3600:.3f} hours.\n")
        
        # # TESTING
        # delta_U_array = np.zeros_like(I_array)
        # delta_U_errors = np.zeros_like(I_array)

        # Compute \delta U \delta U^\dagger term using vegas
        print("Beginning \delta U \delta U^\dagger term.\n")
        t3 = time.time()
        delta_U2_array, delta_U2_errors = compute_delta_U2_term(
            q_array, tau, woods_saxon.occ_states, cg_table, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
        t4 = time.time()
        print(f"Done after {(t4-t3)/3600:.3f} hours.\n")
        
        # # # TESTING
        # delta_U2_array = np.zeros_like(I_array)
        # delta_U2_errors = np.zeros_like(I_array)
        
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
    # kvnn, kmax, kmid, ntot = 6, 30.0, 4.0, 120  # AV18 kmax=30
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120  # AV18 kmax=15
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
        print_normalization=True, save=True
    )