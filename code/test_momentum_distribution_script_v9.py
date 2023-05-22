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

Last update: May 22, 2023

"""

# Python imports
from functools import partial
from jax import config, jit, vmap
import jax.numpy as jnp
import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.special import spherical_jn, sph_harm
import shutil
from sympy.physics.quantum.cg import CG
import time
import vegas

# Imports from scripts
from scripts.integration import momentum_mesh, unattach_weights_from_matrix
from scripts.potentials import Potential
from scripts.srg import load_srg_transformation
from scripts.tools import convert_l_to_string, coupled_channel, replace_periods
from scripts.woodsaxon import ws


# Enable double-precision with JAX arrays
config.update("jax_enable_x64", True)


class SingleParticleState:
    """
    Single-particle state class. Packs together the following single-particle
    quantum numbers into one object.
    
    Parameters
    ----------
    n : int
        Principal quantum number n = 1, 2, ...
    l : int
        Orbital angular momentum l = 0, 1, ...
    j : float
        Total angular momentum j = 1/2, 3/2, ...
    m_j : float
        Total angular momentum projection m_j = -j, -j+1, ..., j.
    m_t : float
        Isospin projection m_t = 1/2 or -1/2.
    
    """
    
    
    def __init__(self, n, l, j, m_j, m_t):
        
        # Check if m_j is valid
        if abs(m_j) > j:
            raise RuntimeError("m_j is not valid.")
            
        # Check that |m_t| = 1/2
        if abs(m_t) != 1/2:
            raise RuntimeError("m_t is not valid.")
            
        self.n = n
        self.l = l
        self.j = j
        self.m_j = m_j
        self.m_t = m_t
        
        if m_t == 1/2:
            self.nucleon = 'proton'
        elif m_t == -1/2:
            self.nucleon = 'neutron'
        
        
    def __eq__(self, sp_state):

        if (
            self.n == sp_state.n and self.l == sp_state.l
            and self.j == sp_state.j and self.m_j == sp_state.m_j
            and self.m_t == sp_state.m_t
        ):
            
            return True
        
        else:
            
            return False
        
        
    def __str__(self):
        
        # Spectroscopic notation of orbital angular momentum
        l_str = convert_l_to_string(self.l)  # E.g., 's', 'p', 'd', ...
        
        # Display j subscript as a fraction
        numerator = 2*int(self.j) + 1
        denominator = 2

        return fr"${self.n}{l_str}_{{{numerator}/{denominator}}}$"
    
    
class PartialWaveChannel:
    """
    Partial wave channel class. Packs together the quantum numbers of a partial
    wave channel into one object.
    
    Parameters
    ----------
    channel : str
        Name of the channel (e.g., '1S0').
    
    """

    
    def __init__(self, channel):
    
        # Set instance attributes
        self.channel = channel
        self.L, self.Lp = self.get_orbital_angular_momentum()
        self.J = self.get_angular_momentum()
        self.S = self.get_spin()
        self.T = self.get_isospin(self.L, self.S)

        
    def __eq__(self, channel):

        if (
            self.L == channel.L and self.Lp == channel.Lp
            and self.J == channel.J and self.S == channel.S
            and self.T == channel.T
        ):
            
            return True
        
        else:
            
            return False
    
    
    def get_orbital_angular_momentum(self):
        """Gets the total orbital angular momentum L and L'."""
        
        if self.channel[1] == 'S':
            L = 0
        elif self.channel[1] == 'P':
            L = 1
        elif self.channel[1] == 'D':
            L = 2
        elif self.channel[1] == 'F':
            L = 3
        elif self.channel[1] == 'G':
            L = 4
        elif self.channel[1] == 'H':
            L = 5
        else:
            raise RuntimeError("Channel L exceeds the range of the function.")
            
        # Coupled-channel
        if coupled_channel(self.channel[:3]):
            
            if self.channel[5] == 'S':
                Lp = 0
            elif self.channel[5] == 'P':
                Lp = 1
            elif self.channel[5] == 'D':
                Lp = 2
            elif self.channel[5] == 'F':
                Lp = 3
            elif self.channel[5] == 'G':
                Lp = 4
            elif self.channel[5] == 'H':
                Lp = 5
            else:
                raise RuntimeError(
                    "Channel L' exceeds the range of the function.")
        
        # L = L' if the channel is not coupled
        else:
            
            Lp = L
            
        return L, Lp
    
    
    def get_angular_momentum(self):
        """Total angular momentum J = 0, 1, 2, ..."""
        
        J = int(self.channel[2])
        
        return J
    
    
    def get_spin(self):
        """Total spin S = 0 or 1."""
        
        S = int((int(self.channel[0])-1)/2)
        
        return S
    
    
    def get_isospin(self, L, S):
        """Total isospin according to antisymmetry."""
    
        # Make sure [1-(-1)^(L+S+T)] is not zero
        if (1 - (-1) ** (L+S)) == 0:
            T = 1
        else:
            T = 0
        
        return T
    
    
class WoodsSaxon:
    """
    Woods-Saxon orbitals class. Handles the wave functions associated with the
    Woods-Saxon potential from the subroutine in woodsaxon.f90. Outputs wave
    functions in coordinate and momentum space.
    
    Parameters
    ----------
    nucleus_name : str
        Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
    Z : int
        Proton number of the nucleus.
    N : int
        Neutron number of the nucleus.
    run_woodsaxon : bool, optional
        Option to run the Woods-Saxon subroutine to generate orbital files.
    n_max : int, optional
        Maximum principal quantum number where n = 1, 2, ..., n_max.
    l_max : int, optional
        Maximum orbital angular momentum where l = 0, 1, ..., l_max.
    rmax : float, optional
        Maximum r for orbital tables.
    ntab : int, optional
        Number of points for orbital tables.
        
    """
    
    
    def __init__(
        self, nucleus_name, Z, N, run_woodsaxon=True, n_max=0, l_max=0, rmax=40,
        ntab=2000
    ):
        
        # Set instance attributes
        self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}/"
        self.dr = rmax / ntab
        self.r_array = np.arange(self.dr, rmax + self.dr, self.dr)

        # Generate orbitals?
        if run_woodsaxon:
            
            self.run_woods_saxon_code(nucleus_name, Z, N, n_max, l_max, rmax,
                                      ntab)

            # Move output files to relevant directory
            shutil.move("ws_log", self.woods_saxon_directory + "ws_log")
            shutil.move("ws_pot", self.woods_saxon_directory + "ws_pot")
            shutil.move("ws_rho", self.woods_saxon_directory + "ws_rho")
                
        # Order single-particle states with lowest energy first
        self.order_sp_states(Z, N)
        
        # Organize wave functions in dictionary with the file name as the key
        self.sp_wfs = {}
        for sp_state in self.sp_states:
            # Wave functions are independent of m_j, so fix m_j=j
            if sp_state.m_j == sp_state.j:
                file_name = get_orbital_file_name(sp_state.n, sp_state.l,
                                                  sp_state.j, sp_state.m_t)
                if run_woodsaxon:
                    shutil.move(file_name,
                                self.woods_saxon_directory + file_name)
                data = np.loadtxt(self.woods_saxon_directory + file_name)
                # Use file name as the key
                self.sp_wfs[file_name] = data[:, 1]

            
    def run_woods_saxon_code(
            self, nucleus_name, Z, N, n_max, l_max, rmax, ntab
    ):
        """Run Woods-Saxon code to generate data."""
        
        # Total number of nucleons
        A = Z + N
        
        # Type of orbitals: 1 - nucleons with no Coulomb
        #                   2 - distinguish protons and neutrons
        ntau = 2
        
        # Orbitals to consider (note, we track 2*j not j)
        norb, lorb, jorb = [], [], []
        for n in range(1, n_max+1):
            for l in range(0, l_max+1):
                norb.append(n)
                lorb.append(l)
                jorb.append(int(2*(l+1/2)))
                if int(2*(l-1/2)) > 0:  # Don't append negative j
                    norb.append(n)
                    lorb.append(l)
                    jorb.append(int(2*(l-1/2)))
        nrad = len(jorb)
        orbws = np.zeros(shape=(2,nrad,ntab), order='F')
    
        # Divide orbital by r? -> get R(r); false: get u(r)=r R(r)
        rdiv = False
        dens = True
    
        # Set parameters of the Woods-Saxon potential
        prm = np.zeros(shape=(2,9), order='F')
    
        # Starting with vws (p & n)
        if nucleus_name == 'He4':
            prm[:,0] = 76.8412
        elif nucleus_name == 'O16':
            prm[:,0] = 58.0611
        elif nucleus_name == 'Ca40':
            prm[:,0] = 54.3051
        elif nucleus_name == 'Ca48':
            prm[0,0] = 59.4522
            prm[1,0] = 46.9322
    
        # Not sure about these (better way to load these parameters?)
        prm[:,1] = 1.275
        prm[:,2] = 0.7
        prm[:,3] = 0.
        prm[:,4] = 1.
        prm[:,5] = 36
        prm[:,6] = 1.32
        prm[:,7] = 0.7
        prm[:,8] = 1.275
        
        # Print summary, potentials, and densities
        prnt = True
        prntorb = True

        # Run Fortran subroutine
        ws(ntau, A, Z, rmax, orbws, norb, lorb, jorb, prm, rdiv, prnt, prntorb,
           dens)
        
        
    def order_sp_states(self, Z, N):
        """Keep track of all s.p. states and occupied s.p. states"""

        self.sp_states = []
        self.occ_states = []
        proton_count = 0
        neutron_count = 0
        
        # File with details of the orbitals
        ws_file = self.woods_saxon_directory + "ws_log"
    
        # Order single-particle states using the ws_log file
        with open(ws_file, 'r') as f:
            for line in f:
                unit = line.strip().split()
                
                # Protons
                if len(unit) > 0 and unit[0] == '1':

                    j = int(unit[3])/2
                    for m_j in np.arange(-j, j+1, 1):
                        sp_state = SingleParticleState(
                            int(unit[1])+1, int(unit[2]), j, m_j, 1/2
                        )  # n, l, j, m_j, m_t
                    
                        self.sp_states.append(sp_state)
                    
                        if proton_count < Z:
                            self.occ_states.append(sp_state)
                            # Add up filled proton states
                            proton_count += 1
                    
                
                # Neutrons
                elif len(unit) > 0 and unit[0] == '2':

                    j = int(unit[3])/2
                    for m_j in np.arange(-j, j+1, 1):
                        sp_state = SingleParticleState(
                            int(unit[1])+1, int(unit[2]), j, m_j, -1/2
                        )  # n, l, j, m_j, m_t
                    
                        self.sp_states.append(sp_state)
                    
                        if neutron_count < N:
                            self.occ_states.append(sp_state)
                            # Add up filled neutron states
                            neutron_count += 1
                        
                        
    def get_wf_rspace(self, sp_state, print_normalization=False):
        """Single-particle wave function in coordinate space."""
        
        # Orbital file name is the key
        u_array = self.sp_wfs[get_orbital_file_name(sp_state.n, sp_state.l,
                                                    sp_state.j, sp_state.m_t)]

        # Normalization: \int dr |u(r)|^2 = 1
        if print_normalization:
            normalization = np.sum(self.dr*u_array**2)
            print(f"Coordinate space normalization = {normalization}.")

        return self.r_array, u_array
    
    
    def fourier_transformation(self, l, k_array):
        """Fourier transformation matrix for given orbital angular momentum."""
        
        # r_array column vectors and k_array row vectors where both grids are
        # n x m matrices
        r_cols, k_rows = np.meshgrid(self.r_array, k_array)
        
        # Transformation matrix with shape n x m, where m is the length of
        # r_array and n is the length of the k_array
        M = 1j**(-l) * np.sqrt(2/np.pi) * self.dr * r_cols * spherical_jn(
            l, k_rows*r_cols
        )
        
        return M
    
    
    def get_wf_kspace(
            self, sp_state, kmax, kmid, ntot, print_normalization=False,
            interpolate=False,
            
    ):
        """Single-particle wave function in momentum space."""
    
        # Set momentum mesh with more points at low momentum
        k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
    
        # Get coordinate-space s.p. wave function
        _, u_array = self.get_wf_rspace(sp_state)

        # Fourier-transform the wave function to momentum space
        phi_array = self.fourier_transformation(sp_state.l, k_array) @ u_array
    
        # Normalization: \int dk k^2 |\phi(k)|^2 = 1
        if print_normalization:
            normalization = np.sum(k_weights*k_array**2*abs(phi_array)**2)
            print(f"Momentum space normalization = {normalization}.")
            
        # Interpolate and return function?
        if interpolate:
            phi_func = interp1d(k_array, phi_array, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
            return phi_func
        
        # Otherwise return momentum, weights, and \phi(k)
        else:
            return k_array, k_weights, phi_array


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


def compute_cg_array(j_max):
    
    j_array = np.arange(0, j_max+1/2, 1/2)
    N_j = j_array.size
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
    
    return jnp.array(cg_array, dtype=jnp.float64), N_j


def cg_mapping(j, m, N_j):
    """Return the indices of the input angular momentum and projection for the
    array of Clebsch-Gordan coefficients.
    """
    
    j_index = jnp.array(j / 0.5, dtype=int)
    m_index = jnp.array(jnp.abs(m/0.5) + jnp.heaviside(-m, 0) * (N_j-1),
                        dtype=int)

    return j_index, m_index


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


def get_orbital_file_name(n, l, j, m_t):
    """Returns the file name of the orbital."""
        
    # Proton
    if m_t == 1/2:
        file_name = f"p.n{int(n-1)}.l{int(l)}.j{int(2*j)}.orb"
    # Neutron
    elif m_t == -1/2:
        file_name = f"n.n{int(n-1)}.l{int(l)}.j{int(2*j)}.orb"
        
    return file_name


def get_sp_wave_functions(sp_basis, kmax, kmid, ntot):
    """Set interpolating functions for s.p. wave functions \phi."""
    
    occ_states = sp_basis.occ_states

    phi_functions = {}
    for sp_state in occ_states: 
        file_name = get_orbital_file_name(sp_state.n, sp_state.l, sp_state.j,
                                          sp_state.m_t)
        phi_functions[file_name] = sp_basis.get_wf_kspace(
            sp_state, kmax, kmid, ntot, interpolate=True)
            
    return phi_functions


def get_phi_function(n, l, j, m_t, k, phi_functions):
        
    return phi_functions[get_orbital_file_name(n, l, j, m_t)](k)


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
    # Y_jml = jnp.where(jnp.abs(m_l) <= l, cg * sph_harm(m_l, l, phi, theta),
    #                   0+0j)
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

    
def compute_I_term(q_array, tau, occ_states, cg_array, N_j):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in occ_states:

            # Single-particle wave function with z-axis along q_vector
            key = get_orbital_file_name(alpha.n, alpha.l, alpha.j, alpha.m_t)
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
        q_array, tau, occ_states, cg_array, N_j, delta_U_quantum_numbers,
        phi_functions, delta_U_functions, delta_U_dagger_functions):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
        
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


def compute_normalization(q_array, q_weights, n_array):
    """Compute the normalization of the momentum distribution."""

    return 4*np.pi * np.sum(q_weights * q_array**2 * n_array)


def save_momentum_distribution(
        nucleus_name, tau, lamb, q_array, q_weights, n_array, n_errors, I_array,
        delta_U_array, delta_U_errors, delta_U2_array, delta_U2_errors
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

    file_name = replace_periods(
        f"{nucleus_name}_{nucleon}_momentum_distribution_{lamb}"
    )
    
    np.savetxt(file_name + '.txt', data, header=hdr)


def load_momentum_distribution(nucleus_name, nucleon, lamb):
    """Load and return the momentum distribution along with the isolated
    contributions.
    """
    
    file_name = replace_periods(
        f"{nucleus_name}_{nucleon}_momentum_distribution_{lamb}"
    )
    
    data = np.loadtxt(file_name + '.txt')
    
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
    # nucleus_name, Z, N, = 'Ca48', 20, 28

    # Nucleon
    tau = 1/2

    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')

    # NN potential and momentum mesh
    # kvnn, kmax, kmid, ntot = 6, 30.0, 4.0, 120  # AV18
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120  # AV18
    # kvnn, kmax, kmid, ntot = 111, 15.0, 3.0, 120  # SMS N4LO 450 MeV

    # SRG \lambda value
    # lamb = 1.35
    lamb = 1.5
    # lamb = 2.0
    # lamb = 3.0
    # lamb = 6.0

    # Compute table of Clebsch-Gordan coefficients
    cg_array, N_j = compute_cg_array(4)
    print("Done calculating Clebsch-Gordan table.\n")

    # Set-up single-particle states
    woods_saxon = WoodsSaxon(nucleus_name, Z, N, run_woodsaxon=False)
    phi_functions = get_sp_wave_functions(woods_saxon, 10.0, 2.0, 120)
    
    # Option to do just the independent particle model (IPM)
    # ipm_only = True
    ipm_only = False
    
    # Option to save the data
    # save = True
    save = False
    
    # Set momentum mesh
    q_array, q_weights = momentum_mesh(10.0, 4.0, 100, nmod=70)
    
    # Compute the I term
    I_array = compute_I_term(q_array, tau, woods_saxon.occ_states, cg_array,
                             N_j)
    
    # Independent particle model has no \delta U contributions
    if ipm_only:
        
        delta_U_array = np.zeros_like(I_array)
        delta_U_errors = np.zeros_like(I_array)
        delta_U2_array = np.zeros_like(I_array)
        delta_U2_errors = np.zeros_like(I_array)
    
    # Include \delta U, \delta U^\dagger, and \delta U \delta U^\dagger
    else:

        t0 = time.time()
        
        # Get an organized list of quantum numbers to sum over
        # delta_U_quantum_numbers = jnp.array(
        #     np.loadtxt('O16_quantum_numbers.txt')
        # )
        delta_U_quantum_numbers = np.loadtxt('O16_quantum_numbers.txt')

        # Set-up \delta U and \delta U^\dagger functions
        delta_U_functions, delta_U_dagger_functions = get_delta_U_functions(
            channels, kvnn, kmax, kmid, ntot, 'Wegner', lamb
        )

        # Compute \delta U + \delta U^\dagger term using vegas
        print("Beginning \delta U linear terms.\n")
        t1 = time.time()
        delta_U_array, delta_U_errors = compute_delta_U_term(
            q_array, tau, woods_saxon.occ_states, cg_array, N_j,
            delta_U_quantum_numbers, phi_functions, delta_U_functions,
            delta_U_dagger_functions
        )
        t2 = time.time()
        print(f"Done after {(t2-t1)/60:.3f} minutes.\n")
        
        # # TESTING
        # delta_U_array = np.zeros_like(I_array)
        # delta_U_errors = np.zeros_like(I_array)

        # # Compute \delta U \delta U^\dagger term using vegas
        # print("Beginning \delta U \delta U^\dagger term.\n")
        # t3 = time.time()
        # delta_U2_array, delta_U2_errors = compute_delta_U2_term(
        #     q_array, tau, woods_saxon.occ_states, cg_array, N_j, pwc_list,
        #     phi_functions, delta_U_functions, delta_U_dagger_functions
        # )
        # t4 = time.time()
        # print(f"Done after {(t4-t3)/60:.3f} minutes.\n")
        
        # TESTING
        delta_U2_array = np.zeros_like(I_array)
        delta_U2_errors = np.zeros_like(I_array)
        
        t5 = time.time()
        print(f"Total time elapsed: {(t5-t0)/60:.3f} minutes.\n")
    
    # Combine each term for the total momentum distribution [fm^3]
    n_array = I_array + delta_U_array + delta_U2_array
    n_errors = np.sqrt(delta_U_errors**2 + delta_U2_errors**2)

    normalization = compute_normalization(q_array, q_weights, n_array)
    print(f"Normalization = {normalization:.5f}.")
    
    if save and not(ipm_only):  # Do not save IPM-only data
        save_momentum_distribution(
            nucleus_name, tau, lamb, q_array, q_weights, n_array, n_errors,
            I_array, delta_U_array, delta_U_errors, delta_U2_array,
            delta_U2_errors
        )