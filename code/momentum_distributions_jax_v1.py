#!/usr/bin/env python3

"""
File: momentum_distributions_jax.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: November 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. This particular version implements
previous codes but in JAX.

Last update: November 7, 2023

"""

# Python imports
from functools import partial
import numpy as np
from scipy.special import spherical_jn
from sympy.physics.quantum.cg import CG
import time
import vegas

# JAX imports
from jax import jit, vmap
from jax.config import config
import jax.numpy as jnp
from jax.numpy.linalg import norm
from jax.scipy.special import sph_harm
from jaxinterp2d import interp2d

# Imports from scripts
from scripts.integration import momentum_mesh, unattach_weights_from_matrix
from scripts.potentials import Potential
from scripts.srg import load_srg_transformation
from scripts.tools import coupled_channel, replace_periods


@jit
def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return jnp.asarray(x == y, dtype=int)


@jit
def build_vector(k, theta, phi):
    """Build a vector from input spherical coordinates."""

    k_vector = jnp.array([k * jnp.sin(theta) * jnp.cos(phi),
                          k * jnp.sin(theta) * jnp.sin(phi),
                          k * jnp.cos(theta)])

    return k_vector


@jit
def get_vector_components(k_vector):
    """Get the spherical coordinates from an input vector."""

    k = norm(k_vector)
    theta = jnp.arccos(k_vector[2]/k)
    phi = jnp.arctan2(k_vector[1], k_vector[0])

    return k, theta, phi


class ClebschGordan:
    
    def __init__(self):

        # Load Clebsch-Gordan coefficients array
        try:
            
            self.cg_array = jnp.load('cg_array.npy')
            self.N_j = self.cg_array.shape[0]
        
        # Compute and save Clebsch-Gordan coefficients array as .npz file
        except FileNotFoundError:
            
            print("Need to compute CG array...")
            # jmax = 6.5  # This should work for nuclei <= Pb208
            # TESTING
            jmax = 4.5
            self.cg_array, self.N_j = self.compute_clebsch_gordan_array(jmax)
            jnp.save('cg_array', self.cg_array)
            print("Done computing CG array.")
        
    def compute_clebsch_gordan_array(self, j_max):
        """Calculate Clebsch-Gordan coefficients for combinations of j and m_j
        up to j_max.
        """

        j_array = np.arange(0, j_max + 1/2, 1/2)
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
    
        return jnp.array(cg_array, dtype=jnp.float64), N_j

    @partial(jit, static_argnums=(0,))
    def m_index(self, m):
        """Returns the unique index associated with m."""
        
        return jnp.asarray(
            jnp.abs(m / 0.5) + jnp.heaviside(-m, 0) * (self.N_j - 1), dtype=int
        )

    @partial(jit, static_argnums=(0,))
    def get_cg(self, j1, m1, j2, m2, j3, m3):
        """Clebsch-Gordan coefficient < j1 m1 j2 m2 | j3 m3 >."""

        i_j1 = jnp.asarray(j1 / 0.5, dtype=int)
        i_m1 = self.m_index(m1)
        i_j2 = jnp.asarray(j2 / 0.5, dtype=int)
        i_m2 = self.m_index(m2)
        i_j3 = jnp.asarray(j3 / 0.5, dtype=int)
        i_m3 = self.m_index(m3)
        
        return self.cg_array[i_j1, i_j2, i_j3, i_m1, i_m2, i_m3]
    
    
class WoodsSaxon:
    """Woods-Saxon single-particle states in coordinate- and momentum-space."""
    
    def __init__(
        self, nucleus_name: str, Z: int, N: int, clebsch_gordan: callable,
        kmax=10.0, kstep=0.1
    ):
        
        # Set instance attributes
        self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}/"
        self.clebsch_gordan = clebsch_gordan
        
        # Order single-particle states with lowest energy first
        occupied_states = self.order_sp_states(Z, N)
        
        # Set-up momentum points for single-particle wave functions
        self.k_array = np.arange(0.0, kmax + kstep, kstep)
        N_k = self.k_array.size

        # Set JAX array of occupied momentum-space wave functions \phi_\alpha(k)
        self.wave_functions = jnp.zeros((Z + N, N_k), dtype=complex)
        for i, occupied_state in enumerate(occupied_states):
            n, l, j, m_j, m_t = occupied_state  # Degenerate in m_j
            phi_k_array = self.compute_wf_kspace(n, l, j, m_t)
            self.wave_functions = self.wave_functions.at[i].set(phi_k_array)
            
        # Set-up mapping for [n, l, j, m_j, m_t] -> single integer index
        self.set_mapping(occupied_states)

    def order_sp_states(self, Z, N):
        """Keep track of all s.p. states and occupied s.p. states"""

        # Occupied single-particle states < E_F
        occupied_states = []
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
                        # (n, l, j, m_j, m_t) tuple
                        sp_state = (
                            int(unit[1])+1, int(unit[2]), j, m_j, 1/2
                        )
                        if proton_count < Z:
                            occupied_states.append(sp_state)
                            # Add up filled proton states
                            proton_count += 1
                
                # Neutrons
                elif len(unit) > 0 and unit[0] == '2':
                    j = int(unit[3])/2
                    for m_j in np.arange(-j, j+1, 1):
                        # (n, l, j, m_j, m_t) tuple
                        sp_state = (
                            int(unit[1])+1, int(unit[2]), j, m_j, -1/2
                        )
                        if neutron_count < N:
                            occupied_states.append(sp_state)
                            # Add up filled neutron states
                            neutron_count += 1
                            
        return occupied_states
    
    def get_wf_rspace(self, n, l, j, m_t, print_norm=False):
        """Single-particle wave function in coordinate space."""
        
        if m_t == 1/2:
            file_name = f"p.n{int(n-1)}.l{int(l)}.j{int(2*j)}.orb"
        elif m_t == -1/2:
            file_name = f"n.n{int(n-1)}.l{int(l)}.j{int(2*j)}.orb"
            
        data = np.loadtxt(self.woods_saxon_directory + file_name)
        
        r_array = data[:, 0]  # fm
        u_array = data[:, 1]  # fm^-1/2
        
        # Normalization: \int dr |u(r)|^2 = 1
        if print_norm:
            normalization = np.trapz(u_array ** 2)
            print(f"Coordinate space normalization = {normalization}.")

        return r_array, u_array

    def compute_wf_kspace(self, n, l, j, m_t):
        """Single-particle wave function in momentum space."""
        
        # Load Woods-Saxon orbital
        r_array, u_array = self.get_wf_rspace(n, l, j, m_t)
        dr = r_array[2] - r_array[1]  # Step-size in r
    
        # Get meshgrids of k and r for integration over r
        k_grid, r_grid = np.meshgrid(self.k_array, r_array, indexing='ij')
        _, u_grid = np.meshgrid(self.k_array, u_array, indexing='ij')

        integrand = r_grid * u_grid * spherical_jn(l, k_grid * r_grid) * dr
        
        return 1j ** (-l) * np.sqrt(2/np.pi) * np.sum(integrand, axis=-1)

    def set_mapping(self, occupied_states):
        """Maps s.p. state to corresponding index in wave functions array."""

        # Get maximum n, l, and j
        self.occupied_states_array = np.array(occupied_states)
        n_max = np.amax(self.occupied_states_array[:, 0])
        l_max = np.amax(self.occupied_states_array[:, 1])
        self.j_max = np.amax(self.occupied_states_array[:, 2])
        
        # Size of each component of the mapping array
        n_size = int(n_max + 1)
        l_size = int(l_max + 1)
        j_size = int(self.j_max + 1/2)
        m_j_size = int(2*j_size)

        # Initialize mapping array
        self.mapping_array = jnp.zeros((n_size, l_size, j_size, m_j_size, 2),
                                       dtype=int)
            
        # Loop over occupied states and fill in index i for state
        for i, occupied_state in enumerate(occupied_states):
            
            n, l, j, m_j, m_t = occupied_state
            n_index = int(n - 1)
            l_index = int(l)
            j_index = int(j - 1/2)
            m_j_index = int(self.j_max + j)
            if m_t == 1/2:
                m_t_index = 0
            elif m_t == -1/2:
                m_t_index = 1

            self.mapping_array = self.mapping_array.at[
                n_index, l_index, j_index, m_j_index, m_t_index
            ].set(i)
    
    @partial(jit, static_argnums=(0,))
    def mapping(self, n, l, j, m_j, m_t):
        """Returns the unique index i of the occupied state given orbital."""
        
        n_index = jnp.asarray(n-1, dtype=int)
        l_index = jnp.asarray(l, dtype=int)
        j_index = jnp.asarray(j - 1/2, dtype=int)
        m_j_index = jnp.asarray(self.j_max + j, dtype=int)
        m_t_index = jnp.where(m_t == 1/2, 0, 1)
        
        return self.mapping_array[n_index, l_index, j_index, m_j_index,
                                  m_t_index]

    @partial(jit, static_argnums=(0,))
    def spinor_spherical_harmonic(self, theta, phi, sigma, l, j, m_j):
        """Spinor spherical harmonic for a s.p. state described by the quantum
        numbers j, m_j, l, and s=1/2.
        """
        
        # Convert \theta and \phi to arrays
        theta = jnp.array([theta])
        phi = jnp.array([phi])
        
        # Spinor indexed by \sigma \eta_{m_s}^(\sigma) = \delta_{m_s, \sigma}
        m_s = sigma
    
        # m_l must be fixed since m_j and m_s are determined
        m_l = m_j - m_s
        
        # Make sure that |m_l| <= l
        condition = jnp.abs(m_l) <= l

        return jnp.where(
            condition,
            (
                self.clebsch_gordan.get_cg(l, m_l, 1/2, m_s, j, m_j)
                * sph_harm(
                    jnp.array([m_l], dtype=int), jnp.array([l], dtype=int), phi,
                    theta, n_max=10
                )
            ),
            jnp.zeros_like(theta, dtype=complex)
        )
    
    @partial(jit, static_argnums=(0,))
    def psi(self, k, theta, phi, sigma, tau, n, l, j, m_j, m_t):
        """Single-particle wave function \psi_\alpha(k_vector; \sigma, \tau)."""

        # Get \phi_\alpha(k)
        index = self.mapping(n, l, j, m_j, m_t)
        phi_k_array = self.wave_functions[index]
        
        # Interpolate to point k
        phi_k = jnp.interp(k, self.k_array, phi_k_array)
        
        # Calculate spinor spherical harmonic
        Y_jml = self.spinor_spherical_harmonic(theta, phi, sigma, l, j, m_j)
    
        # Isospinor indexed by \tau \chi_{m_t}(\tau)
        chi_tau = kronecker_delta(tau, m_t)

        return phi_k * Y_jml * chi_tau
    
    
class PartialWaveChannel:
    """Organizes the quantum numbers of a partial wave channel."""

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
    

class DeltaUMatrixElement:
    """Plane-wave matrix element of \delta U = U - I."""
    
    def __init__(self, kvnn: int, kmax: float, kmid: float, ntot: int,
                 generator: str, lamb: float, clebsch_gordan: callable):
        
        # Set instance attributes
        self.clebsch_gordan = clebsch_gordan
        k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
        self.k_array = jnp.asarray(k_array)
        self.k_weights = jnp.asarray(k_weights)

        # Arrays of partial wave channel quantum numbers
        self.channel_strings = ['1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1',
                                '3D1-3D1']
        self.channels = self.set_partial_wave_channels()
        N_channels = self.channels.shape[0]
        
        # Set JAX array of partial-wave matrix elements of \delta U
        self.deltaU_pw_me = jnp.zeros((N_channels, ntot, ntot))
        self.deltaUdag_pw_me = jnp.zeros((N_channels, ntot, ntot))
        for i, channel_str in enumerate(self.channel_strings):
            delta_U = self.get_delta_U(kvnn, channel_str, kmax, kmid, ntot,
                                       generator, lamb, hc=False)
            delta_U_dag = self.get_delta_U(kvnn, channel_str, kmax, kmid, ntot,
                                           generator, lamb, hc=True)
            self.deltaU_pw_me = self.deltaU_pw_me.at[i].set(delta_U)
            self.deltaUdag_pw_me = self.deltaUdag_pw_me.at[i].set(delta_U_dag)
            
        # Set-up mapping for [L, Lp, S, J, T] -> single integer index
        self.set_mapping(self.channels)
        
    def set_partial_wave_channels(self):
        """Sets all S-wave partial-wave channels."""
        
        channels = []
        for channel_str in self.channel_strings:
            
            # Get quantum numbers of channel
            pwc = PartialWaveChannel(channel_str)
            L, Lp, S, J, T = pwc.L, pwc.Lp, pwc.S, pwc.J, pwc.T
            
            # Loop over possible M_J values
            M_J_array = np.arange(-J, J+1)
            for M_J in M_J_array:
            
                channels.append([L, Lp, S, J, M_J, T])
            
        return jnp.array(channels, dtype=int)
    
    def get_delta_U(self, kvnn, channel_str, kmax, kmid, ntot, generator, lamb,
                    hc=False):
        """Get partial-wave matrix elements of \delta U."""
        
        if channel_str[:3] in ['3S1', '3D1']:
            potential = Potential(kvnn, '3S1', kmax, kmid, ntot)
        else:
            potential = Potential(kvnn, '1S0', kmax, kmid, ntot)

        # Get SRG transformation with integration weights [unitless]
        U_matrix_weights = load_srg_transformation(potential, generator, lamb)

        # Calculate \delta U = U - I
        I_matrix_weights = jnp.eye(len(U_matrix_weights))
        if hc:  # Hermitian conjugate
            delU_matrix_weights = (U_matrix_weights - I_matrix_weights).T
        else:
            delU_matrix_weights = U_matrix_weights - I_matrix_weights

        # Get specific sub-block if coupled-channel and unattach weights
        if channel_str == '3S1-3D1':
            delU_matrix = unattach_weights_from_matrix(
                self.k_array, self.k_weights, delU_matrix_weights[:ntot,ntot:]
            )
        elif channel_str == '3D1-3S1':
            delU_matrix = unattach_weights_from_matrix(
                self.k_array, self.k_weights, delU_matrix_weights[ntot:,:ntot]
            )
        elif channel_str == '3D1-3D1':
            delU_matrix = unattach_weights_from_matrix(
                self.k_array, self.k_weights, delU_matrix_weights[ntot:,ntot:]
            )
        else:
            delU_matrix = unattach_weights_from_matrix(
                self.k_array, self.k_weights, delU_matrix_weights[:ntot,:ntot]
            )
            
        return delU_matrix
        
    def set_mapping(self, channels):
        """Maps partial-wave channel to corresponding index in \delta U partial-
        wave matrix elements array.
        """

        # Get maximum L, L', S, J, and T
        L_max = jnp.amax(self.channels[:, 0])
        Lp_max = jnp.amax(self.channels[:, 1])
        S_max = jnp.amax(self.channels[:, 2])
        self.J_max = jnp.amax(self.channels[:, 3])
        T_max = jnp.amax(self.channels[:, 4])
        
        # Size of each component of the mapping array
        L_size = int(L_max + 1)
        Lp_size = int(Lp_max + 1)
        S_size = int(S_max + 1)
        J_size = int(self.J_max + 1)
        M_J_size = int(2 * J_size)
        T_size = int(T_max + 1)

        # Initialize mapping array
        self.mapping_array = jnp.zeros((L_size, Lp_size, S_size, J_size,
                                        M_J_size, T_size), dtype=int)
            
        # Loop over channels and fill in index i for state
        for i, channel in enumerate(channels):
            
            L, Lp, S, J, M_J, T = channel
            L_index = int(L)
            Lp_index = int(Lp)
            S_index = int(S)
            J_index = int(J)
            M_J_index = int(self.J_max + J)
            T_index = int(T)

            self.mapping_array = self.mapping_array.at[
                L_index, Lp_index, S_index, J_index, M_J_index, T_index
            ].set(i)
    
    @partial(jit, static_argnums=(0,))
    def mapping(self, L, Lp, S, J, T):
        """Returns the unique index i of the partial-wave channel."""
        
        L_index = jnp.asarray(L, dtype=int)
        Lp_index = jnp.asarray(Lp, dtype=int)
        S_index = jnp.asarray(S, dtype=int)
        J_index = jnp.asarray(J, dtype=int)
        M_J_index = jnp.asarray(self.J_max + J, dtype=int)
        T_index = jnp.asarray(T, dtype=int)
        
        return self.mapping_array[L_index, Lp_index, S_index, J_index,
                                  M_J_index, T_index]
    
    @partial(jit, static_argnums=(0,))
    def delta_U_partial_wave(self, k, kp, hc, L, Lp, S, J, T):
        """Evaluate partial-wave matrix element of \delta U at k and k'."""
        
        # Get \phi_\alpha(k)
        index = self.mapping(L, Lp, S, J, T)
        
        # Partial-wave matrix element (with or without \dagger)
        delta_U = jnp.where(hc, self.deltaUdag_pw_me[index],
                            self.deltaU_pw_me[index])

        # Interpolate to point (k, k')
        return interp2d(k, kp, self.k_array, self.k_array, delta_U)
    
    @partial(jit, static_argnums=(0,))
    def __call__(
            self, k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2,
            sigma_3, sigma_4, tau_1, tau_2, tau_3, tau_4, hc
    ):
        """Sum over partial wave channels."""

        return jnp.sum(self.vmap_channels(
            k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2, sigma_3,
            sigma_4, tau_1, tau_2, tau_3, tau_4, hc, self.channels)
        )
    
    @partial(jit, static_argnums=(0,))
    def vmap_channels(
            self, k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2,
            sigma_3, sigma_4, tau_1, tau_2, tau_3, tau_4, hc, channel
    ):
        """vmap of check_selection_rules method with shape (N_channels,)."""
        
        return vmap(
            self.check_selection_rules,
            in_axes=(None, None, None, None, None, None, None, None, None, None, 
                     None, None, None, None, None, 0)
        )(k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2, sigma_3,
          sigma_4, tau_1, tau_2, tau_3, tau_4, hc, channel)

    @partial(jit, static_argnums=(0,))
    def check_selection_rules(
            self, k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2,
            sigma_3, sigma_4, tau_1, tau_2, tau_3, tau_4, hc, channel
    ):
        """Check that quantum number selection rules are satisfied and evaluate
        the matrix element.
        """
        
        # Unpack partial-wave channel quantum numbers
        L, Lp, S, J, M_J, T = channel
        
        # Set total spin and isospin projections
        M_T = tau_1 + tau_2
        M_S = sigma_1 + sigma_2
        M_Sp = sigma_3 + sigma_4
        
        # Set total orbital angular momentum projection
        M_L = M_J - M_S
        M_Lp = M_J - M_Sp
        
        # Make sure quantum numbers are consistent
        isospin_conservation = tau_1 + tau_2 == tau_3 + tau_4
        selection_rules = jnp.logical_and(
            isospin_conservation, jnp.logical_and(
                jnp.abs(M_T) <= T, jnp.logical_and(
                    jnp.abs(M_S) <= S, jnp.logical_and(
                        jnp.abs(M_S) <= S, jnp.logical_and(
                            jnp.abs(M_Sp) <= S, jnp.logical_and(
                                jnp.abs(M_J) <= J, jnp.logical_and(
                                    jnp.abs(M_L) <= L, jnp.abs(M_Lp) <= Lp
                                )
                            )
                        )
                    )
                )
            )
        )
        
        return jnp.where(
            selection_rules,
            self.evaluate(
                k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2,
                sigma_3, sigma_4, tau_1, tau_2, tau_3, tau_4, hc, channel, M_T,
                M_S, M_Sp, M_L, M_Lp
            ),
            jnp.zeros_like(k)
        )
        
    @partial(jit, static_argnums=(0,))
    def evaluate(
            self, k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2,
            sigma_3, sigma_4, tau_1, tau_2, tau_3, tau_4, hc, channel, M_T, M_S,
            M_Sp, M_L, M_Lp
    ):
        """Evaluate the matrix element of \delta U given the input quantum
        numbers.
        """
        
        # Unpack partial-wave channel quantum numbers
        L, Lp, S, J, M_J, T = channel

        # Get spin Clebsch-Gordan coefficients
        spin_12_cg = self.clebsch_gordan.get_cg(1/2, sigma_1, 1/2, sigma_2, S,
                                                M_S)
        spin_34_cg = self.clebsch_gordan.get_cg(1/2, sigma_3, 1/2, sigma_4, S,
                                                M_Sp)
    
        # Get isospin Clebsch-Gordan coefficients
        isospin_12_cg = self.clebsch_gordan.get_cg(1/2, tau_1, 1/2, tau_2, T,
                                                   M_T)
        isospin_34_cg = self.clebsch_gordan.get_cg(1/2, tau_3, 1/2, tau_4, T,
                                                   M_T)
                
        # \delta U or \delta U^\dagger in partial-wave basis
        delU_partial_wave = self.delta_U_partial_wave(k, kp, hc, L, Lp, S, J, T)
                    
        # Get L-S coupling Clebsch-Gordan coefficients
        lsj_cg = self.clebsch_gordan.get_cg(L, M_L, S, M_S, J, M_J)
        lpsj_cg = self.clebsch_gordan.get_cg(Lp, M_Lp, S, M_Sp, J, M_J)
                        
        # Calculate spherical harmonics
        Y_k = sph_harm(
            jnp.array([M_L], dtype=int), jnp.array([L], dtype=int),
            jnp.array([phi_k]), jnp.array([theta_k]), n_max=10
        )
        Y_kp = sph_harm(
            jnp.array([M_Lp], dtype=int), jnp.array([Lp], dtype=int),
            jnp.array([phi_kp]), jnp.array([theta_kp]), n_max=10
        )

        # Factors of [1-(-1)^(L+S+T)] [1-(-1)^(L'+S+T)] included
        matrix_element = (
            1.0/2.0 * 2.0/jnp.pi * 4.0 * spin_12_cg * spin_34_cg
            * isospin_12_cg * isospin_34_cg * delU_partial_wave * lsj_cg
            * lpsj_cg * Y_k * jnp.conj(Y_kp)
        )
                    
        return matrix_element
    

class DeltaUIntegrand:
    """Integrand of \delta U + \delta U^\dagger terms."""
    
    def __init__(self, tau: float, woods_saxon: callable,
                 delta_U_matrix_element: callable):
        
        # Set instance attributes
        self.tau = tau
        self.woods_saxon = woods_saxon
        self.delta_U_matrix_element = delta_U_matrix_element
        
        # Arrays of isospins [(\tau', \tau_1, \tau_2), ...]
        self.isospin_configurations = self.set_isospin_configurations()
        
        # Arrays of spins [(\sigma_1, \sigma_2, \sigma, \sigma', ...)]
        self.spin_configurations = self.set_spin_configurations()
        
        # Arrays of occupied Woods-Saxon orbitals [(n, l, j, m_j, m_t), ...]
        self.occupied_states = woods_saxon.occupied_states_array
        
    def set_spin_configurations(self):
        """Spin projection configurations for \delta U + \delta U^\dagger term.
        """
    
        spin_projections = [1/2, -1/2]
        spin_configurations = []
        for spin_1 in spin_projections:
            for spin_2 in spin_projections:
                for spin_3 in spin_projections:
                    for spin_4 in spin_projections:
                        spin_configurations.append((spin_1, spin_2, spin_3,
                                                    spin_4))
                        
        return jnp.array(spin_configurations)
    
    def set_isospin_configurations(self):
        """Isospin projection configurations for \delta U + \delta U^\dagger
        term.
        """
    
        isospin_projections = [1/2, -1/2]
        isospin_configurations = []
        for taup in isospin_projections:
            for tau_1 in isospin_projections:
                for tau_2 in isospin_projections:
                    # Check that M_T is conserved
                    if self.tau + taup == tau_1 + tau_2:
                        isospin_configurations.append((taup, tau_1, tau_2))
    
        return jnp.array(isospin_configurations)
    
    @partial(jit, static_argnums=(0,))
    def __call__(self, q, x_array):
        """Sum over \alpha."""

        return jnp.sum(self.vmap_alpha(q, x_array, self.occupied_states))
    
    @partial(jit, static_argnums=(0,))
    def vmap_alpha(self, q, x_array, alpha):
        """vmap of sum_over_beta method with shape (N_occupied,)."""

        return vmap(self.sum_over_beta, in_axes=(None, None, 0))(q, x_array,
                                                                 alpha)
    
    @partial(jit, static_argnums=(0,))
    def sum_over_beta(self, q, x_array, alpha):
        """Sum over \beta."""

        return jnp.sum(self.vmap_beta(q, x_array, alpha, self.occupied_states))
    
    @partial(jit, static_argnums=(0,))
    def vmap_beta(self, q, x_array, alpha, beta):
        """vmap of sum_over_spin method with shape (N_occupied,)."""

        return vmap(self.sum_over_spin, in_axes=(None, None, None, 0))(
            q, x_array, alpha, beta
        )
    
    @partial(jit, static_argnums=(0,))
    def sum_over_spin(self, q, x_array, alpha, beta):
        """Sum over \sigma_1, \sigma_2, \sigma, \sigma'."""
        
        return jnp.sum(self.vmap_spin(q, x_array, alpha, beta,
                                      self.spin_configurations))
    
    @partial(jit, static_argnums=(0,))
    def vmap_spin(self, q, x_array, alpha, beta, spin):
        """vmap of sum_over_isospins method with shape (N_spin,)."""

        return vmap(self.sum_over_isospin, in_axes=(None, None, None, None, 0))(
            q, x_array, alpha, beta, spin
        )
    
    @partial(jit, static_argnums=(0,))
    def sum_over_isospin(self, q, x_array, alpha, beta, spin):
        """Sum over \tau', \tau_1, and \tau_2."""

        return jnp.sum(self.vmap_isospin(q, x_array, alpha, beta, spin,
                                         self.isospin_configurations))
    
    @partial(jit, static_argnums=(0,))
    def vmap_isospin(self, q, x_array, alpha, beta, spin, isospin):
        """vmap of evaluate method with shape (N_isospin,)."""

        return vmap(self.evaluate, in_axes=(None, None, None, None, None, 0))(
            q, x_array, alpha, beta, spin, isospin
        )
    
    @partial(jit, static_argnums=(0,))
    def evaluate(self, q, x_array, alpha, beta, spin, isospin):
        """Evaluate integrand for particular point in momenta and configuration
        of quantum numbers. This method assumes sampling over the following:
        
            k, theta_k, phi_k -> vector momenta k,
            K, theta_K, phi_K -> vector momenta K,
            
        where x_array contains each variable meaning x_array.shape = (6,).
        """

        # Unpack quantum numbers
        n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha = alpha
        n_beta, l_beta, j_beta, m_j_beta, m_t_beta = beta
        sigma_1, sigma_2, sigma, sigmap = spin
        taup, tau_1, tau_2 = isospin

        # Choose z-axis to be along q_vector
        # q_vector = jnp.array([0.0, 0.0, self.q])
        # TESTING
        q_vector = jnp.array([0.0, 0.0, q])
        
        # Relative momenta k
        k, theta_k, phi_k = x_array[0:3]
        k_vector = build_vector(k, theta_k, phi_k)
        
        # C.o.M. momenta K
        K, theta_K, phi_K = x_array[3:6]
        K_vector = build_vector(K, theta_K, phi_K)

        # Calculate vector q - K/2
        qK_vector = q_vector - K_vector / 2.0
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
        
        # Calculate vector k_1 = K/2 + k
        k1_vector = K_vector / 2.0 + k_vector
        k1, theta_k1, phi_k1 = get_vector_components(k1_vector)
        
        # Calculate vector k_2 = K/2 - k
        k2_vector = K_vector / 2.0 - k_vector
        k2, theta_k2, phi_k2 = get_vector_components(k2_vector)
        
        # Calculate vector K - q
        Kq_vector = K_vector - q_vector
        Kq, theta_Kq, phi_Kq = get_vector_components(Kq_vector)
        
        # Calculate the Jacobian determinant
        jacobian = k ** 2 * jnp.sin(theta_k) * K ** 2 * jnp.sin(theta_K)
  
        # Plane-wave matrix elements of \delta U and \delta U^\dagger
        delta_U_plane_wave = self.delta_U_matrix_element(
            k, theta_k, phi_k, qK, theta_qK, phi_qK, sigma_1, sigma_2, sigma,
            sigmap, tau_1, tau_2, self.tau, taup, False
        )
        delta_U_dag_plane_wave = self.delta_U_matrix_element(
            qK, theta_qK, phi_qK, k, theta_k, phi_k, sigma, sigmap, sigma_1,
            sigma_2, self.tau, taup, tau_1, tau_2, True
        )
        
        # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
        psi_alpha_1 = self.woods_saxon.psi(
            k1, theta_k1, phi_k1, sigma_1, tau_1, n_alpha, l_alpha, j_alpha,
            m_j_alpha, m_t_alpha
        )
        # \psi_\beta(K/2-k; \sigma_2, \tau_2)
        psi_beta_2 = self.woods_saxon.psi(
            k2, theta_k2, phi_k2, sigma_2, tau_2, n_beta, l_beta, j_beta,
            m_j_beta, m_t_beta
        )
        
        # \psi_\alpha(q; \sigma, \tau)
        # psi_alpha_q = self.woods_saxon.psi(
        #     self.q, 0.0, 0.0, sigma, tau, n_alpha, l_alpha, j_alpha, m_j_alpha,
        #     m_t_alpha
        # )
        psi_alpha_q = self.woods_saxon.psi(
            q, 0.0, 0.0, sigma, self.tau, n_alpha, l_alpha, j_alpha, m_j_alpha,
            m_t_alpha
        )
        # \psi_\beta(K-q; \sigma', \tau')
        psi_beta_Kq = self.woods_saxon.psi(
            Kq, theta_Kq, phi_Kq, sigmap, taup, n_beta, l_beta, j_beta,
            m_j_beta, m_t_beta
        )
        
        # \psi_\alpha(K-q; \sigma', \tau')
        psi_alpha_Kq = self.woods_saxon.psi(
            Kq, theta_Kq, phi_Kq, sigmap, taup, n_alpha, l_alpha, j_alpha,
            m_j_alpha, m_t_alpha
        )
        # \psi_\beta(q; \sigma, \tau)
        # psi_beta_q = self.woods_saxon.psi(
        #     self.q, 0.0, 0.0, sigma, tau, n_beta, l_beta, j_beta, m_j_beta,
        #     m_t_beta
        # )
        psi_beta_q = self.woods_saxon.psi(
            q, 0.0, 0.0, sigma, self.tau, n_beta, l_beta, j_beta, m_j_beta,
            m_t_beta
        )
        
        # \delta U term
        delta_U_term = (
            delta_U_plane_wave * jnp.conj(psi_alpha_1) * jnp.conj(psi_beta_2)
            * (psi_beta_Kq * psi_alpha_q - psi_alpha_Kq * psi_beta_q)
        )
        
        # \delta U^\dagger term
        delta_U_dag_term = (
            delta_U_dag_plane_wave * psi_alpha_1 * psi_beta_2
            * (
                jnp.conj(psi_beta_Kq) * jnp.conj(psi_alpha_q)
                - jnp.conj(psi_alpha_Kq) * jnp.conj(psi_beta_q)
            )
        )

        # Add together for full integrand
        integrand = jacobian * (delta_U_term + delta_U_dag_term) / 2.0
        
        return integrand.real
    
    
class DeltaU2Integrand:
    
    def __init__(self):
        pass
    
    def __call__(self):
        pass
    
    
class MomentumDistribution:
    """Single-nucleon momentum distribution using SRG transformations and
    a Woods-Saxon mean-field treatment.
    """
    
    def __init__(
        self, nucleus_name, Z, N, kvnn, lamb, kmax=15.0, kmid=3.0, ntot=120,
        generator='Wegner'
    ):
        
        # Instance attributes
        self.nucleus_name = nucleus_name
        self.kvnn = kvnn
        self.lamb = lamb
        
        # Set table of Clebsch-Gordan coefficients
        self.clebsch_gordan = ClebschGordan()
        
        # Set single-particle basis
        self.woods_saxon = WoodsSaxon(nucleus_name, Z, N, self.clebsch_gordan)
        
        # Set \delta U matrix element function
        self.delta_U_matrix_element = DeltaUMatrixElement(
            kvnn, kmax, kmid, ntot, generator, lamb, self.clebsch_gordan)
        
    @partial(jit, static_argnums=(0,))
    def psi_q(self, q, sigma, tau, alpha):
        """Compute \psi_\alpha(q;\sigma,\tau)."""
        
        # Unpack quantum numbers
        n, l, j, m_j, m_t = alpha
        
        # Single-particle wave function with z-axis along q_vector
        theta = jnp.zeros_like(q)
        phi = jnp.zeros_like(q)
        psi_alpha = self.woods_saxon.psi(q, theta, phi, sigma, tau, n, l, j,
                                         m_j, m_t)

        return psi_alpha
    
    @partial(jit, static_argnums=(0,))
    def vmap_ipm_alpha(self, q, sigma, tau, alpha):
        """vmap of ipm_term method with shape (N_alpha,)."""
        
        return vmap(self.psi_q, in_axes=(None, None, None, 0))(q, sigma, tau, 
                                                               alpha)
    
    @partial(jit, static_argnums=(0,))
    def sum_over_alpha(self, tau, q):
        """Sum over \alpha and \sigma in the IPM term."""
        
        # Evaluate IPM on array of occupied \alpha
        alpha_array = self.woods_saxon.occupied_states_array
        
        # Sum over \sigma = +/- 1/2
        ipm_array = (
            self.vmap_ipm_alpha(q, jnp.asarray(1/2), tau, alpha_array)
            + self.vmap_ipm_alpha(q, jnp.asarray(-1/2), tau, alpha_array)
        )
                
        # Sum over \alpha
        return jnp.sum(ipm_array)
    
    @partial(jit, static_argnums=(0,))
    def vmap_q_ipm(self, tau, q):
        """vmap of sum_over_alpha method with shape (N_q,)."""
        
        return vmap(self.sum_over_alpha, in_axes=(None, 0))(tau, q)
    
    @partial(jit, static_argnums=(0,))
    def ipm_term(self, tau, q_array):
        """Compute the I * n(q) * I term."""
        
        return jnp.abs(self.vmap_q_ipm(tau, q_array)) ** 2
    
    def delU_term(self, tau, q_array, neval):
        """Compute \delta U * n(q) * I + I * n(q) * \delta U^\dagger terms."""

        # Relative momenta from 0 to 10 fm^-1
        k_limits = [0, 10]
        # C.o.M. momenta up to 3 fm^-1
        K_limits = [0, 3]
        # Polar angle from 0 to \pi
        theta_limits = [0, np.pi]
        # Azimuthal angle from 0 to 2\pi
        phi_limits = [0, 2*np.pi]

        # Set-up integrator with multiple processors
        integ = vegas.Integrator([
            k_limits, theta_limits, phi_limits,
            K_limits, theta_limits, phi_limits
        ], nproc=8)
    
        print("Starting \delta U + \delta U^\dagger term...\n")
        
        # TESTING
        delta_U_integrand = DeltaUIntegrand(tau, self.woods_saxon,
                                            self.delta_U_matrix_element)

        # Evaluate the \delta U + \delta U^\dagger term for each q
        delta_U_array = np.zeros_like(q_array)
        delta_U_errors = np.zeros_like(q_array)
        for i, q in enumerate(q_array):

            # TESTING
            t0 = time.time()
        
            # integrand = DeltaUIntegrand(tau, q, self.woods_saxon,
            #                             self.delta_U_matrix_element)
            # TESTING
            integrand = partial(delta_U_integrand, q)
        
            # Train the integrator
            integ(integrand, nitn=5, neval=neval)
    
            # Final result
            result = integ(integrand, nitn=10, neval=neval)
        
            delta_U_array[i] = result.mean
            delta_U_errors[i] = result.sdev
        
            # TESTING
            t1 = time.time()
            mins = (t1-t0)/60
            percent = (i+1) / q_array.size * 100
            print(f"{percent:.1f} percent done after {mins:.3f} minutes.")
        
        print("\nDone with \delta U + \delta U^\dagger term.")

        return delta_U_array, delta_U_errors
    
    def delU2_term(self, tau, q_array, neval):
        """Compute the \delta U * n(q) * \delta U^\dagger term."""

        # Relative momenta from 0 to 10 fm^-1
        k_limits = [0, 10]
        # C.o.M. momenta up to 3 fm^-1
        K_limits = [0, 3]
        # Polar angle from 0 to \pi
        theta_limits = [0, np.pi]
        # Azimuthal angle from 0 to 2\pi
        phi_limits = [0, 2*np.pi]

        # Set-up integrator with multiple processors
        integ = vegas.Integrator([
            k_limits, theta_limits, phi_limits,
            k_limits, theta_limits, phi_limits,
            K_limits, theta_limits, phi_limits
        ], nproc=8)
    
        print("Starting \delta U \delta U^\dagger term...\n")

        # Evaluate the \delta U \delta U^\dagger term for each q
        delta_U2_array = np.zeros_like(q_array)
        delta_U2_errors = np.zeros_like(q_array)
        for i, q in enumerate(q_array):

            # TESTING
            t0 = time.time()
        
            integrand = DeltaU2Integrand(tau, q, self.woods_saxon,
                                         self.delta_U_matrix_element)
        
            # Train the integrator
            integ(integrand, nitn=5, neval=neval)
    
            # Final result
            result = integ(integrand, nitn=10, neval=neval)
        
            delta_U2_array[i] = result.mean
            delta_U2_errors[i] = result.sdev
        
            # TESTING
            t1 = time.time()
            mins = (t1-t0)/60
            percent = (i+1) / q_array.size * 100
            print(f"{percent:.1f} percent done after {mins:.3f} minutes.")
        
        print("\nDone with \delta U \delta U^\dagger term.")

        return delta_U2_array, delta_U2_errors
        
    def compute_normalization(self, q_array, q_weights, n_array):
        """Normalization of single-nucleon momentum distribution."""
        
        return 4 * np.pi * np.sum(q_weights * q_array ** 2 * n_array)
    
    def evaluate(
            self, tau, delU_neval, delU2_neval, qmax=10.0, qmid=3.0, ntotq=100,
            nmodq=60, print_normalization=False, save=False
    ):
        """Compute the single-nucleon momentum distribution."""
    
        # Set points in q
        q_array, q_weights = momentum_mesh(qmax, qmid, ntotq, nmod=nmodq)
    
        # Compute the I term
        I_array = self.ipm_term(tau, q_array).block_until_ready()
        print("Done with IPM term.")
        
        # Compute the \delta U + \delta U^\dagger terms
        delta_U_array, delta_U_errors = self.delU_term(tau, q_array, delU_neval)
        
        # Compute the \delta U \delta U^\dagger term
        # delta_U2_array, delta_U2_errors = self.delU2_term(tau, q_array, delU2_neval)
        # TESTING
        delta_U2_array = np.zeros_like(I_array)
        delta_U2_errors = np.zeros_like(I_array)
        
        n_array = I_array + delta_U_array + delta_U2_array
        n_errors = np.sqrt(delta_U_errors ** 2 + delta_U2_errors ** 2)
        
        # Option to print normalization of the total momentum distribution
        if print_normalization:
            normalization = self.compute_normalization(q_array, q_weights,
                                                       n_array)
            if tau == 1/2:
                print(f"Z = {normalization:.5f}.")
            elif tau == -1/2:
                print(f"N = {normalization:.5f}.")
    
        # Option to save the momentum distribution as a .txt file
        if save:
            self.save_momentum_distribution(
                tau, q_array, q_weights, n_array, n_errors, I_array,
                delta_U_array, delta_U_errors, delta_U2_array, delta_U2_errors
            )
    
        return q_array, q_weights, n_array, n_errors
    
    def save_momentum_distribution(
            self, tau, q_array, q_weights, n_array, n_errors, I_array,
            delta_U_array, delta_U_errors, delta_U2_array, delta_U2_errors
    ):
        """Save the momentum distribution along with the isolated
        contributions.
        """
    
        data = np.vstack(
            (q_array, q_weights, n_array, n_errors, I_array, delta_U_array,
             delta_U_errors, delta_U2_array, delta_U2_errors)
        ).T
            
        if tau == 1/2:
            nucleon = 'proton'
        elif tau == -1/2:
            nucleon = 'neutron'
                
        hdr = (
            "q, q weight, n(q), n(q) error, I, \delta U + \delta U^\dagger,"
            " \delta U + \delta U^\dagger error, \delta U^2, \delta U^2 error\n"
        )
    
        directory = 'momentum_distributions/'

        file_name = replace_periods(
            f"{self.nucleus_name}_{nucleon}_momentum_distribution"
            f"_kvnn_{self.kvnn}_lamb_{self.lamb}"
        )
    
        np.savetxt(directory + file_name + '.txt', data, header=hdr)
        
    def load_momentum_distribution(self, nucleon):
        """Load and return the momentum distribution along with the isolated
        contributions.
        """
    
        directory = 'momentum_distributions/'

        file_name = replace_periods(
            f"{self.nucleus_name}_{nucleon}_momentum_distribution"
            f"_kvnn_{self.kvnn}_lamb_{self.lamb}"
        )
    
        data = np.loadtxt(directory + file_name + '.txt')

        return data


if __name__ == '__main__':
    
    # Enable double precision
    config.update("jax_enable_x64", True)
    
    # Nucleus
    nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'C12', 6, 6
    # nucleus_name, Z, N = 'O16', 8, 8
    # nucleus_name, Z, N = 'Ca40', 20, 20
    # nucleus_name, Z, N = 'Ca48', 20, 28
    # nucleus_name, Z, N = 'Fe54', 26, 28
    # nucleus_name, Z, N = 'Sn118', 50, 68
    # nucleus_name, Z, N, = 'Pb208', 82, 126
    
    # Nucleon
    tau = 1/2
    # tau = -1/2
    
    # NN potential and momentum mesh
    # kvnn, kmax, kmid, ntot = 5, 15.0, 3.0, 120  # Nijmegen II
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120  # AV18
    # kvnn, kmax, kmid, ntot = 7, 15.0, 3.0, 120  # CD-Bonn
    # kvnn, kmax, kmid, ntot = 79, 10.0, 2.0, 120  # EMN N4LO 500 MeV
    # kvnn, kmax, kmid, ntot = 111, 10.0, 2.0, 120  # SMS N4LO 450 MeV
    # kvnn, kmax, kmid, ntot = 222, 10.0, 2.0, 120  # GT+ N2LO 1 fm
    
    
    # SRG \lambda value
    # lamb = 1.35
    lamb = 1.5
    # lamb = 1.7
    # lamb = 2.0
    # lamb = 3.0
    # lamb = 6.0
    
    # Max evaluations of the integrand
    delU_neval, delU2_neval = 1e3, 5e3
    # delU_neval, delU2_neval = 1e4, 5e4
    # delU_neval, delU2_neval = 5e4, 1e5
    # delU_neval, delU2_neval = 1e5, 5e5
    # delU_neval, delU2_neval = 2.5e5, 7.5e5
    # delU_neval, delU2_neval = 3e5, 1e6

    # Compute and save the momentum distribution
    md = MomentumDistribution(nucleus_name, Z, N, kvnn, lamb)
    q_array, q_weights, n_array, n_errors = md.evaluate(tau, delU_neval,
                                                        delU2_neval, save=False)