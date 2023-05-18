#!/usr/bin/env python3

"""
File: test_momentum_distribution_script.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. This particular version evaluates the
sum over all non-zero combinations of partial wave channels and single-particle
quantum numbers as an additional integral and also uses JAX to speed-up
various functions.

Last update: May 17, 2023

"""

# Python imports
from functools import partial
from jax import config, jit
import jax.numpy as jnp
from numba import njit
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
    Single-particle basis class. Handles the wave functions associated with the
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
        key = get_orbital_file_name(sp_state.n, sp_state.l, sp_state.j,
                                    sp_state.m_t)
        u_array = self.sp_wfs[key]

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
        key = get_orbital_file_name(sp_state.n, sp_state.l, sp_state.j,
                                    sp_state.m_t)
        phi_functions[key] = sp_basis.get_wf_kspace(
            sp_state, kmax, kmid, ntot, interpolate=True)
            
    return phi_functions


def interpolate_delta_U(
        kvnn, channel, kmax, kmid, ntot, generator, lamb,
        hermitian_conjugate=False
):
    """Interpolate \delta U(k, k') for the given channel."""
    
    # Set channel argument to be compatible with potential functions
    if channel[:3] in ['3S1', '3D1']:
        channel_arg = '3S1'
    elif channel[:3] in ['3P2', '3F2']:
        channel_arg = '3P2'
    elif channel[:3] in ['3D3', '3G3']:
        channel_arg = '3D3'
    else:
        channel_arg = channel[:3]

    # Get momentum mesh
    k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
    
    potential = Potential(kvnn, channel_arg, kmax, kmid, ntot)
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
        
        pwc = PartialWaveChannel(channel)
        key = (pwc.L, pwc.Lp, pwc.J, pwc.S, pwc.T)
        
        delta_U_functions[key] = interpolate_delta_U(kvnn, channel, kmax, kmid,
                                                     ntot, generator, lamb)
        delta_U_dagger_functions[key] = interpolate_delta_U(
            kvnn, channel, kmax, kmid, ntot, generator, lamb,
            hermitian_conjugate=True
        )
        
    return delta_U_functions, delta_U_dagger_functions


@jit
def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return jnp.array(x == y, dtype=int)


@jit
# @njit
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
    # k_vector = np.array([k * np.sin(theta) * np.cos(phi),
    #                      k * np.sin(theta) * np.sin(phi),
    #                      k * np.cos(theta)])

    return k_vector


@jit
# @njit
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

    k = jnp.linalg.norm(k_vector)
    theta = jnp.arccos(k_vector[2]/k)
    phi = jnp.arctan2(k_vector[1], k_vector[0])
    # k = np.linalg.norm(k_vector)
    # theta = np.arccos(k_vector[2]/k)
    # phi = np.arctan2(k_vector[1], k_vector[0])

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


def psi(n, l, j, m_j, m_t, k, theta, phi, sigma, tau, cg_array, N_j,
        phi_functions):
    """Single-particle wave function."""

    # Calculate \phi_\alpha(q)
    key = get_orbital_file_name(n, l, j, m_t)
    phi_sp_wf = phi_functions[key](k)
    
    # Calculate spinor spherical harmonic
    Y_jml = spinor_spherical_harmonic(l, j, m_j, theta, phi, sigma, cg_array,
                                      N_j)
    
    # Isospinor indexed by \tau \chi_{m_t}(\tau)
    chi_tau = kronecker_delta(tau, m_t)

    return phi_sp_wf * Y_jml * chi_tau


def spinor_spherical_harmonic(l, j, m_j, theta, phi, sigma, cg_array, N_j):
    """Spinor spherical harmonic for a s.p. state described by the quantum
    numbers j, m_j, l, and s=1/2.
    """

    # Spinor indexed by \sigma \eta_{m_s}^(\sigma) = \delta_{m_s, \sigma}
    m_s = sigma
    
    # m_l must be fixed since m_j and m_s are determined
    # m_l = m_j - m_s
    m_l = np.array(m_j - m_s).astype(int)
    l = np.array(l).astype(int)
        
    # Check that |m_l| <= l
    if np.abs(m_l) <= l:
        
        # Clebsch-Gordan coefficient for l-s coupling
        cg = clebsch_gordan_coefficient(l, m_l, 1/2, m_s, j, m_j, cg_array, N_j)
        
        # Spherical harmonic
        Y_lm = sph_harm(m_l, l, phi, theta)
            
        return cg * Y_lm
    
    else:
        
        return 0+0j
    

def compute_I_term(q_array, tau, occ_states, cg_array, N_j, phi_functions):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in occ_states:

            # Single-particle wave function with z-axis along q_vector
            psi_alpha_array = psi(
                alpha.n, alpha.l, alpha.j, alpha.m_j, alpha.m_t, q_array, 0, 0,
                sigma, tau, cg_array, N_j, phi_functions
            )

            I_array += np.abs(psi_alpha_array) ** 2
                    
    return I_array


@partial(jit, static_argnames=['n'])
def get_quantum_numbers(x, quantum_numbers, n):
    
    index = np.floor(x * n).astype(int)
    return quantum_numbers[index]


def delta_U_term_integrand(
        q, tau, delta_U_quantum_numbers, N_delU, cg_array, N_j, phi_functions,
        delta_U_functions, delta_U_dagger_functions, x_array
):
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""

    # Choose z-axis to be along q_vector
    q_vector = np.array([0, 0, q])

    # Relative momenta k
    k, theta_k, phi_k = x_array[:3]
    k_vector = build_vector(k, theta_k, phi_k)
        
    # C.o.M. momenta K
    K, theta_K, phi_K = x_array[3:6]
    K_vector = build_vector(K, theta_K, phi_K)
    
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
    
    # Sample of quantum numbers
    # quantum_numbers = get_quantum_numbers(x_array[6], delta_U_quantum_numbers,
    #                                       N_delU)
    index = np.floor(x_array[6] * N_delU).astype(int)
    quantum_numbers = delta_U_quantum_numbers[index]
    
    # Unpack quantum numbers
    sigma_1, sigma_2, sigma, sigmap = quantum_numbers[0:4]
    tau_1, tau_2, taup = quantum_numbers[4:7]
    n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha = quantum_numbers[7:12]
    n_beta, l_beta, j_beta, m_j_beta, m_t_beta = quantum_numbers[12:17]
    S, M_S, M_Sp = quantum_numbers[17:20]
    L, M_L, Lp, M_Lp = quantum_numbers[20:24]               
    J, M_J = quantum_numbers[24:26]
    T, M_T = quantum_numbers[26:]                
        
    # < \sigma_1 \sigma_2 | S M_S >
    spin_12_cg = clebsch_gordan_coefficient(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                            cg_array, N_j)
    
    # < S M_S' | \sigma \sigma' >
    spin_ssp_cg = clebsch_gordan_coefficient(1/2, sigma, 1/2, sigmap, S, M_Sp,
                                             cg_array, N_j)
    
    # < \tau_1 \tau_2 | T M_T >
    isospin_12_cg = clebsch_gordan_coefficient(1/2, tau_1, 1/2, tau_2, T, M_T,
                                               cg_array, N_j)
    
    # < T M_T | \tau \tau' >
    isospin_ttp_cg = clebsch_gordan_coefficient(1/2, tau, 1/2, taup, T, M_T,
                                                cg_array, N_j)
        
    # < L M_L S M_S | J M_J >
    lsj_cg = clebsch_gordan_coefficient(L, M_L, S, M_S, J, M_J, cg_array, N_j)
    
    # < J M_J | L' M_L' S M_S' >
    lpsj_cg = clebsch_gordan_coefficient(Lp, M_Lp, S, M_Sp, J, M_J, cg_array,
                                         N_j)
    
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
    key_1 = (L, Lp, J, S, T)
    delta_U_partial_wave = delta_U_functions[key_1].ev(k, qK)
    
    # < |q-K/2| (L' S) J T | \delta U^\dagger | k (L S) J T >
    key_2 = (Lp, L, J, S, T)
    delta_U_dag_partial_wave = delta_U_dagger_functions[key_2].ev(qK, k)
    
    # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
    psi_alpha_1 = psi(
        n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
        k1, theta_k1, phi_k1, sigma_1, tau_1,
        cg_array, N_j, phi_functions
    )
    
    # \psi_\beta(K/2-k; \sigma_2, \tau_2)
    psi_beta_2 = psi(
        n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
        k2, theta_k2, phi_k2, sigma_2, tau_2,
        cg_array, N_j, phi_functions
    )
    
    # \psi_\alpha(q; \sigma, \tau)
    psi_alpha_q = psi(
        n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
        q, 0, 0, sigma, tau,
        cg_array, N_j, phi_functions
    )
    
    # \psi_\beta(K-q; \sigma', \tau')
    psi_beta_Kq = psi(
        n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
        Kq, theta_Kq, phi_Kq, sigmap, taup,
        cg_array, N_j, phi_functions
    )
    
    # \psi_\beta(q; \sigma, \tau)
    psi_beta_q = psi(
        n_beta, l_beta, j_beta, m_j_beta, m_t_beta,
        q, 0, 0, sigma, tau,
        cg_array, N_j, phi_functions
    )
    
    # \psi_\alpha(K-q; \sigma', \tau')
    psi_alpha_Kq = psi(
        n_alpha, l_alpha, j_alpha, m_j_alpha, m_t_alpha,
        Kq, theta_Kq, phi_Kq, sigmap, taup,
        cg_array, N_j, phi_functions
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
    integrand = 1/2 * 1/2 * 2/np.pi * jacobian * N_delU * (
        delta_U_term + delta_U_dag_term
    )

    return integrand.real


def compute_delta_U_term(
        q_array, tau, occ_states, cg_array, N_j, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions,
):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    delta_U_quantum_numbers = np.loadtxt('O16_quantum_numbers.txt')
    N_delU = len(delta_U_quantum_numbers)
        
    # Relative momenta from 0 to maximum of momentum mesh
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    # Polar angle
    theta_limits = [0, np.pi]
    # Azimuthal angle
    phi_limits = [0, 2*np.pi]
    # Sum over quantum numbers using vegas integration trick
    qn_limits = [0, 1]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits,
                              qn_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = partial(
            delta_U_term_integrand, q, tau, delta_U_quantum_numbers, N_delU,
            cg_array, N_j, phi_functions, delta_U_functions,
            delta_U_dagger_functions,
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
    cg_array, N_j = compute_cg_array(4)
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
            q_array, tau, woods_saxon.occ_states, cg_array, N_j, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
        t2 = time.time()
        print(f"Done after {(t2-t1)/3600:.3f} hours.\n")
        
        # # TESTING
        # delta_U_array = np.zeros_like(I_array)
        # delta_U_errors = np.zeros_like(I_array)

        # # Compute \delta U \delta U^\dagger term using vegas
        # print("Beginning \delta U \delta U^\dagger term.\n")
        # t3 = time.time()
        # delta_U2_array, delta_U2_errors = compute_delta_U2_term(
        #     q_array, tau, woods_saxon.occ_states, cg_table, channels,
        #     phi_functions, delta_U_functions, delta_U_dagger_functions
        # )
        # t4 = time.time()
        # print(f"Done after {(t4-t3)/3600:.3f} hours.\n")
        
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
            nucleus_name, tau, lamb, q_array, q_weights, n_array, n_errors,
            I_array, delta_U_array, delta_U_errors, delta_U2_array,
            delta_U2_errors
        )
    
    return q_array, q_weights, n_array, n_errors


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
    
    # Nucleon
    tau = 1/2
    
    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    
    # NN potential and momentum mesh
    kvnn, kmax, kmid, ntot = 6, 30.0, 4.0, 120  # AV18
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