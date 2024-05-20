#!/usr/bin/env python3

"""
File: test_spectroscopic_overlaps.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 14, 2024

This script serves as a testbed for calculating spectroscopic overlaps using
mean field approximations for initial and final states and applying SRG
transformations to the operator.

Last update: May 20, 2024

"""

# Python imports
from numba import njit
import numpy as np
import numpy.linalg as la
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from scipy.special import sph_harm, spherical_jn
from sympy.physics.quantum.cg import CG
import time
import vegas

# Imports from scripts
from scripts.integration import momentum_mesh, unattach_weights_from_matrix
from scripts.potentials import Potential
from scripts.srg import compute_srg_transformation, load_srg_transformation
from scripts.tools import convert_l_to_string, coupled_channel, replace_periods


class SingleParticleState:
    """
    Single-particle state class. Packs together the following single-particle
    quantum numbers into one object.
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
   
    
class WoodsSaxon:
    """
    Woods-Saxon orbitals class. Handles the radial wave functions associated
    with the Woods-Saxon potential from the subroutine in woodsaxon.f90. Outputs
    radial functions in coordinate and momentum space.
    
    Parameters
    ----------
    nucleus_name : str
        Name of the nucleus (e.g., 'O16', 'Ca40', etc.)
    Z : int
        Proton number of the nucleus.
    N : int
        Neutron number of the nucleus.
    run_woods_saxon : bool, optional
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
    
    def __init__(self, nucleus_name, Z, N, cg_table, rmax=40, ntab=2000,
                 kmax=10.0, kmid=2.0, ntot=120, parametrization='seminole'):
        
        # Set instance attributes
        self.woods_saxon_directory = (
            f"../data/woods_saxon/{parametrization}/{nucleus_name}/"
        )
        self.cg_table = cg_table
        self.A = int(Z + N)

        # Order single-particle states with lowest energy first
        self.order_sp_states(Z, N)
        
        # Organize wave functions in dictionary
        self.sp_wfs = {}
        for sp_state in self.sp_states:
            
            # Wave functions are independent of m_j, so fix m_j=j
            if sp_state.m_j == sp_state.j:
                
                file_name = self.get_orbital_file_name(sp_state)
 
                data = np.loadtxt(self.woods_saxon_directory + file_name)
                
                # Use n, l, j, m_t as the key
                key = (sp_state.n, sp_state.l, sp_state.j, sp_state.m_t)
                self.sp_wfs[key] = data[:, 1]
                
        # r_array and dr are the same for every s.p. state
        self.r_array = data[:, 0]
        self.dr = max(self.r_array) / len(self.r_array)
        
        # Interpolate occupied s.p. wave functions w.r.t. momentum k
        self.interpolate_wavefunctions(kmax, kmid, ntot)
        
        
    def get_orbital_file_name(self, sp_state):
        """Returns the file name of the orbital."""
            
        # Proton
        if sp_state.m_t == 1/2:
            file_name = (f"p.n{int(sp_state.n-1)}.l{int(sp_state.l)}"
                          f".j{int(2*sp_state.j)}.orb")
        # Neutron
        elif sp_state.m_t == -1/2:
            file_name = (f"n.n{int(sp_state.n-1)}.l{int(sp_state.l)}"
                          f".j{int(2*sp_state.j)}.orb")

        return file_name
        
        
    def order_sp_states(self, Z, N):
        """Keep track of all s.p. states and occupied s.p. states"""

        self.sp_states = []  # All single-particle states
        self.occupied_states = []  # Occupied single-particle states < E_F
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
                            self.occupied_states.append(sp_state)
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
                            self.occupied_states.append(sp_state)
                            # Add up filled neutron states
                            neutron_count += 1


    def get_wf_rspace(self, sp_state, print_normalization=False):
        """Single-particle wave function in coordinate space."""
        
        # n, l, j, m_t
        key = (sp_state.n, sp_state.l, sp_state.j, sp_state.m_t) 
        u_array = self.sp_wfs[key]

        # Normalization: \int dr |u(r)|^2 = 1
        if print_normalization:
            normalization = np.sum(self.dr * u_array ** 2)
            print(f"Coordinate space normalization = {normalization}.")
        
        # Return r and u(r)
        return self.r_array, u_array
    
    
    def fourier_transformation(self, l, r_array, k_array):
        """Fourier transformation matrix for given orbital angular momentum."""
        
        # r_array column vectors and k_array row vectors where both grids are
        # n x m matrices
        r_cols, k_rows = np.meshgrid(r_array, k_array)
        
        # Transformation matrix with shape n x m, where m is the length of
        # r_array and n is the length of the k_array
        M = 1j ** (-l) * np.sqrt(2/np.pi) * self.dr * r_cols * spherical_jn(
            l, k_rows*r_cols
        )
        
        return M
    
    
    def get_wf_kspace(
            self, sp_state, kmax, kmid, ntot, print_normalization=False,
            interpolate=False
    ):
        """Single-particle wave function in momentum space."""
    
        # Set momentum mesh with more points at low momentum
        k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
    
        # Get coordinate-space s.p. wave function
        r_array, u_array = self.get_wf_rspace(sp_state)

        # Fourier-transform the wave function to momentum space
        phi_array = (self.fourier_transformation(sp_state.l, r_array, k_array)
                     @ u_array)
    
        # Normalization: \int dk k^2 |\phi(k)|^2 = 1
        if print_normalization:
            normalization = np.sum(k_weights * k_array ** 2
                                   * abs(phi_array) ** 2)
            print(f"Momentum space normalization = {normalization}.")
            
        # Interpolate and return function?
        if interpolate:
            
            if sp_state.l % 2 == 0:  # Even l
                
                phi_func = InterpolatedUnivariateSpline(k_array, phi_array.real)
                
            else:  # Odd l
            
                phi_func = InterpolatedUnivariateSpline(k_array, phi_array.imag)
                
            return phi_func
        
        # Otherwise return momentum, weights, and \phi(k)
        else:
            return k_array, k_weights, phi_array
        
        
    def interpolate_wavefunctions(self, kmax, kmid, ntot):
        """Create dictionary of \phi(k) interpolated functions where the key
        is the single-particle state.
        """
        
        # Organize interpolated wave functions in dictionary with s.p. quantum
        # numbers as the key
        self.phi_functions = {}
        
        for sp_state in self.occupied_states:
            
            key = (sp_state.n, sp_state.l, sp_state.j, sp_state.m_t)
            self.phi_functions[key] = self.get_wf_kspace(
                sp_state, kmax, kmid, ntot, interpolate=True
            )
            
    
    def psi(self, sp_state, k, theta, phi, sigma, tau):
        """Single-particle wave function \psi_\alpha(k_vector; \sigma, \tau)."""

        # Calculate \phi_\alpha(k)
        key = (sp_state.n, sp_state.l, sp_state.j, sp_state.m_t)
        if sp_state.l % 2 == 0:  # Even l
            phi_sp_wf = self.phi_functions[key](k)
        else:  # Odd l needs factor of i^-l
            phi_sp_wf = 1j ** (-sp_state.l) * self.phi_functions[key](k)
    
        # Calculate spinor spherical harmonic
        Y_jml = self.spinor_spherical_harmonic(sp_state.l, sp_state.j,
                                               sp_state.m_j, theta, phi, sigma)
    
        # Isospinor indexed by \tau \chi_{m_t}(\tau)
        chi_tau = kronecker_delta(tau, sp_state.m_t)

        return phi_sp_wf * Y_jml * chi_tau
    

    def spinor_spherical_harmonic(self, l, j, m_j, theta, phi, sigma):
        """Spinor spherical harmonic for a s.p. state described by the quantum
        numbers j, m_j, l, and s=1/2.
        """
        
        # Spinor indexed by \sigma \eta_{m_s}^(\sigma) = \delta_{m_s, \sigma}
        m_s = sigma
    
        # m_l must be fixed since m_j and m_s are determined
        m_l = m_j - m_s
        
        # Check that |m_l| <= l
        if np.abs(m_l) <= l:
        
            # Clebsch-Gordan coefficient for l-s coupling
            cg = self.cg_table[(l, m_l, 1/2, m_s, j, m_j)]
        
            # Spherical harmonic
            Y_lm = sph_harm(m_l, l, phi, theta)
            
            return cg * Y_lm
        
        else:
            
            return np.zeros_like(theta, dtype=complex)
        
    
class PartialWaveChannel:
    """
    Partial wave channel class. Packs together the quantum numbers of a partial
    wave channel into one object.
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
    
    
class DeltaUMatrixElement:
    """Computes the plane-wave matrix elements of \delta U or \delta U^\dagger.
    Allows for inverse-SRG transformations to match operators of different
    interactions.
    """
    
    def __init__(self, cg_table, kvnn, kmax, kmid, ntot, generator, lamb,
                 channels, kvnn_hard=None, lambda_m=None):
        
        # Set instance attributes
        self.cg_table = cg_table
        self.channels = channels
        
        # Set spin singlet and triplet channels
        self.spin_singlet_channels = self.get_singlet_channels()
        self.spin_triplet_channels = self.get_triplet_channels()
        
        # Set \delta U and \delta U^\dagger functions
        self.get_delta_U_functions(kvnn, kmax, kmid, ntot, generator, lamb,
                                   kvnn_hard, lambda_m)
    
    
    def __call__(
            self, k, theta_k, phi_k, kp, theta_kp, phi_kp, sigma_1, sigma_2,
            sigma_3, sigma_4, tau_1, tau_2, tau_3, tau_4, hc=False
    ):
        """Sum over partial waves to evaluate \delta U matrix element in
        plane-wave basis.
        """
        
        # matrix_element = 0+0j
        matrix_element = np.zeros_like(k, dtype=complex)
        
        # Check that M_T is conserved
        if tau_1 + tau_2 != tau_3 + tau_4:
            return matrix_element
        else:
            M_T = tau_1 + tau_2
            
        # Set total spin projections
        M_S = sigma_1 + sigma_2
        M_Sp = sigma_3 + sigma_4
        
        # Loop over only spin triplet channels
        if np.abs(M_S) == 1 or np.abs(M_Sp) == 1:
            channels = self.spin_triplet_channels
        # Loop over only spin singlet channels
        # NOTE: THIS WON'T WORK FOR P-WAVES CORRECTLY
        elif np.abs(M_T) == 1:
            channels = self.spin_singlet_channels
        # Loop over all channels
        else:
            channels = self.channels
            
        for channel in channels:

            # Get quantum numbers of channel
            pwc = PartialWaveChannel(channel)
            L, Lp, S, J, T = pwc.L, pwc.Lp, pwc.S, pwc.J, pwc.T
            
            # Check that |M_T| <= T and |M_S|, |M_S'| <= S
            if np.abs(M_T) <= T and np.abs(M_S) <= S and np.abs(M_Sp) <= S:
                
                # Get spin Clebsch-Gordan coefficients
                spin_12_cg = self.cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
                spin_34_cg = self.cg_table[(1/2, sigma_3, 1/2, sigma_4, S,
                                            M_Sp)]
    
                # Get isospin Clebsch-Gordan coefficients
                isospin_12_cg = self.cg_table[(1/2, tau_1, 1/2, tau_2, T, M_T)]
                isospin_34_cg = self.cg_table[(1/2, tau_3, 1/2, tau_4, T, M_T)]
                
                # \delta U^\dagger in partial-wave basis
                key = (L, Lp, J, S, T)
                if hc:
                    delta_U_partial_wave = self.delUdag_functions[key].ev(k, kp)
                # \delta U in partial-wave basis
                else:
                    delta_U_partial_wave = self.delU_functions[key].ev(k, kp)
            
                # Sum over M_J
                for M_J in np.arange(-J, J+1):
                    
                    # Set total orbital angular momentum projection
                    M_L = M_J - M_S
                    M_Lp = M_J - M_Sp
                    
                    # Check |M_L| <= L conditions
                    if np.abs(M_L) <= L and np.abs(M_Lp) <= Lp:
                        
                        # Get L-S coupling Clebsch-Gordan coefficients
                        lsj_cg = self.cg_table[(L, M_L, S, M_S, J, M_J)]
                        lpsj_cg = self.cg_table[(Lp, M_Lp, S, M_Sp, J, M_J)]
                        
                        # Calculate spherical harmonics
                        Y_k = sph_harm(M_L, L, phi_k, theta_k)
                        Y_kp = sph_harm(M_Lp, Lp, phi_kp, theta_kp)

                        # Factors of [1-(-1)^(L+S+T)] [1-(-1)^(L'+S+T)] included
                        matrix_element += (
                            1/2 * 2/np.pi * 4 * spin_12_cg * spin_34_cg
                            * isospin_12_cg * isospin_34_cg
                            * delta_U_partial_wave * lsj_cg * lpsj_cg
                            * Y_k * np.conj(Y_kp)
                        )
    
        return matrix_element
    
    
    def interpolate_delta_U(
            self, channel, potential, generator, lamb, hc=False,
            potential_hard=None, lambda_m=None
    ):
        """Interpolate \delta U(k, k') for the given channel."""

        # Get momentum mesh
        kmax, kmid, ntot = potential.kmax, potential.kmid, potential.ntot
        k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
        
        # Inverse-SRG-evolved Hamiltonian as starting point
        if lambda_m is not None:
            H_initial = potential.load_hamiltonian()
            H_hard_initial = potential_hard.load_hamiltonian()
            if generator == 'Block-diag':
                H_evolved = potential.load_hamiltonian('srg', generator, 1.0,
                                                       lambda_bd=lamb)
                H_hard_evolved = potential_hard.load_hamiltonian(
                    'srg', generator, 1.0, lambda_bd=lambda_m)
            else:
                H_evolved = potential.load_hamiltonian('srg', generator, lamb)
                H_hard_evolved = potential_hard.load_hamiltonian(
                    'srg', generator, lambda_m)
            # Get SRG transformation for inverse transformation
            U_hard = compute_srg_transformation(H_hard_initial, H_hard_evolved)
            # Do inverse transformation on initial Hamiltonian
            H_initial = U_hard.T @ H_initial @ U_hard
            # Get SRG transformation with integration weights [unitless]
            U_matrix_weights = compute_srg_transformation(H_initial, H_evolved)
            
        else:
            # Get SRG transformation with integration weights [unitless]
            U_matrix_weights = load_srg_transformation(potential, generator,
                                                       lamb)

        # Calculate \delta U = U - I
        I_matrix_weights = np.eye(len(U_matrix_weights))
        if hc:  # Hermitian conjugate
            delU_matrix_weights = (U_matrix_weights - I_matrix_weights).T
        else:
            delU_matrix_weights = U_matrix_weights - I_matrix_weights

        # Get specific sub-block if coupled-channel and unattach weights
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
        
        
    def get_delta_U_functions(
            self, kvnn, kmax, kmid, ntot, generator, lamb, kvnn_hard=None,
            lambda_m=None
    ):
        """Get \delta U and \delta U^\dagger functions."""
        
        self.delU_functions = {}
        self.delUdag_functions = {}
        for channel in self.channels:
        
            # Set channel argument to be compatible with potential functions
            if channel[:3] == '3D1':
                channel_arg = '3S1'
            elif channel[:3] == '3F2':
                channel_arg = '3P2'
            elif channel[:3] == '3G3':
                channel_arg == '3D3'
            else:
                channel_arg = channel[:3]
            
            # Potential in specified partial wave channel
            potential = Potential(kvnn, channel_arg, kmax, kmid, ntot)
            if kvnn_hard is not None:
                potential_hard = Potential(kvnn_hard, channel_arg, kmax, kmid,
                                           ntot)
            else:
                potential_hard = None
            
            # Use partial wave channel as key
            pwc = PartialWaveChannel(channel)
            key = (pwc.L, pwc.Lp, pwc.J, pwc.S, pwc.T)
        
            # Set \delta U and \delta U^\dagger functions
            self.delU_functions[key] = self.interpolate_delta_U(
                channel, potential, generator, lamb,
                potential_hard=potential_hard, lambda_m=lambda_m
            )
            self.delUdag_functions[key] = self.interpolate_delta_U(
                channel, potential, generator, lamb, hc=True,
                potential_hard=potential_hard, lambda_m=lambda_m
            )
            
            
    def get_singlet_channels(self):
        """Set all S=0 partial wave channels."""
        
        spin_singlet_channels = []
        for channel in self.channels:
            pwc = PartialWaveChannel(channel)
            if pwc.S == 0:
                spin_singlet_channels.append(channel)
        
        return spin_singlet_channels
                
        
    def get_triplet_channels(self):
        """Set all S=1 partial wave channels."""
        
        spin_triplet_channels = []
        for channel in self.channels:
            pwc = PartialWaveChannel(channel)
            if pwc.S == 1:
                spin_triplet_channels.append(channel)         
        
        return spin_triplet_channels
    

class DeltaUDaggerIntegrand:
    """Evaluate the integrand of the \delta U^\dagger term at several
    integration points.
    """
    
    def __init__(
            self, q, alpha, psi_states, overlap_matrix_inv, overlap_matrix_det,
            spin_configurations, isospin_configurations, delta_U_matrix_element,
            ws_psi, ws_phi
    ):
        """Set class attributes for evaluation of integrand."""
        
        self.q = q
        self.alpha = alpha
        self.psi_states = psi_states
        self.overlap_matrix_inv = overlap_matrix_inv
        self.overlap_matrix_det = overlap_matrix_det
        self.spin_configurations = spin_configurations
        self.isospin_configurations = isospin_configurations
        self.delta_U_matrix_element = delta_U_matrix_element
        self.ws_psi = ws_psi
        self.ws_phi = ws_phi
    
    def alpha_index(self, alpha):
        """Obtain the index of the s.p. state \alpha within the list containing
        occupied states of < \Psi |.
        """
        
        index = 0
        for j, sp_state in enumerate(self.psi_states):
            if sp_state == alpha:
                index = j
                break
                
        return index
    
    def evaluate(self, x_array):
        """Evaluates the integrand of the \delta U^\dagger term."""

        # Relative momenta k
        k = x_array[:, 0]
        theta_k = x_array[:, 1]
        phi_k = x_array[:, 2]
        k_vector = build_vector(k, theta_k, phi_k)
        
        # C.o.M. momenta K
        K = x_array[:, 3]
        theta_K = x_array[:, 4]
        phi_K = x_array[:, 5]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Choose z-axis to be along q_vector
        q_vector = np.zeros_like(k_vector)
        q_vector[-1, :] = self.q

        # Calculate vector q - K/2
        qK_vector = q_vector - K_vector/2
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
        
        # Calculate vector k_1 = K/2 + k
        k1_vector = K_vector/2 + k_vector
        k1, theta_k1, phi_k1 = get_vector_components(k1_vector)
        
        # Calculate vector k_2 = K/2 - k
        k2_vector = K_vector/2 - k_vector
        k2, theta_k2, phi_k2 = get_vector_components(k2_vector)
        
        # Calculate vector K - q
        Kq_vector = K_vector - q_vector
        Kq, theta_Kq, phi_Kq = get_vector_components(Kq_vector)
        
        # Calculate the Jacobian determinant
        jacobian = k ** 2 * np.sin(theta_k) * K ** 2 * np.sin(theta_K)
        
        # Index of \alpha in occupied states of < \Psi_0^A(\lambda) |
        jalpha = self.alpha_index(self.alpha)
        
        # Overlap matrix M_{\beta, \alpha} (will be used for \beta \alpha and
        # \beta' \alpha)
        M_ba = self.overlap_matrix_inv[jalpha, jalpha]
        
        # Sum over spin projections
        integrand = 0+0j
        for spin_projections in self.spin_configurations:
            
            sigma, sigmap, sigma_1, sigma_2 = spin_projections
            
            # Sum over isospin projections
            for isospin_projections in self.isospin_configurations:
                
                tau, taup, tau_1, tau_2 = isospin_projections
                
                # Plane-wave matrix elements of \delta U^\dagger
                delta_U_dag = self.delta_U_matrix_element(
                    qK, theta_qK, phi_qK, k, theta_k, phi_k, sigma, sigmap,
                    sigma_1, sigma_2, tau, taup, tau_1, tau_2, hc=True
                )
                
                # \psi_\beta(K/2+k; \sigma_1, \tau_1)
                psi_beta_1 = self.ws_phi.psi(self.alpha, k1, theta_k1, phi_k1,
                                             sigma_1, tau_1)
                
                # \psi_\beta(K/2-k; \sigma_2, \tau_2)
                psi_beta_2 = self.ws_phi.psi(self.alpha, k2, theta_k2, phi_k2,
                                             sigma_2, tau_2)
                
                # Sum over occupied states \alpha'
                for alphap in self.ws_psi.occupied_states:
                    
                    # Index of \alpha'
                    jalphap = self.alpha_index(alphap)
                    
                    # Overlap matrix M_{\beta, \alpha'} (will be used for
                    # \beta' \alpha' and \beta \alpha')
                    M_bap = self.overlap_matrix_inv[jalphap, jalphap]
                    
                    # \psi_\alpha'(K-q; \sigma', \tau')
                    psi_alphap = self.ws_psi.psi(alphap, Kq, theta_Kq, phi_Kq,
                                                 sigmap, taup)

                    # \psi_\beta'(K/2-k; \sigma_2, \tau_2)
                    psi_betap_2 = self.ws_phi.psi(alphap, k2, theta_k2, phi_k2,
                                                  sigma_2, tau_2)
                    
                    # \psi_\beta'(K/2+k; \sigma_1, \tau_1)
                    psi_betap_1 = self.ws_phi.psi(alphap, k1, theta_k1, phi_k1,
                                                  sigma_1, tau_1)
                    
                    # Add together for full integrand
                    integrand += 1/2 * (
                        self.overlap_matrix_det ** 2 * jacobian * delta_U_dag
                        * psi_alphap * M_ba * M_bap * (
                            np.conj(psi_beta_1) * np.conj(psi_betap_2)
                            - np.conj(psi_beta_2) * np.conj(psi_betap_1)
                        )
                    )

        return integrand


@vegas.batchintegrand
class Real(DeltaUDaggerIntegrand):
    """Real part of the \delta U^\dagger integrand."""
    
    def __init__(
            self, q, alpha, psi_states, overlap_matrix_inv, overlap_matrix_det,
            spin_configurations, isospin_configurations, delta_U_matrix_element,
            ws_psi, ws_phi
    ):
        """Inherit the attributes and methods of DeltaUDaggerIntegrand."""
        
        super().__init__(
            q, alpha, psi_states, overlap_matrix_inv, overlap_matrix_det,
            spin_configurations, isospin_configurations, delta_U_matrix_element,
            ws_psi, ws_phi
        )
        
    def __call__(self, x_array):
        """Take real part of the integrand."""
        
        return self.evaluate(x_array).real
    
    
@vegas.batchintegrand
class Imag(DeltaUDaggerIntegrand):
    """Imaginary part of the \delta U^\dagger integrand."""
    
    def __init__(
            self, q, alpha, psi_states, overlap_matrix_inv, overlap_matrix_det,
            spin_configurations, isospin_configurations, delta_U_matrix_element,
            ws_psi, ws_phi
    ):
        """Inherit the attributes and methods of DeltaUDaggerIntegrand."""
        
        super().__init__(
            q, alpha, psi_states, overlap_matrix_inv, overlap_matrix_det,
            spin_configurations, isospin_configurations, delta_U_matrix_element,
            ws_psi, ws_phi
        )
        
    def __call__(self, x_array):
        """Take imaginary part of the integrand."""
        
        return self.evaluate(x_array).imag
    

# TODO: Make option to do C11 instead of B11 based on input m_t
# TODO: Add real and complex capability to integrand
class SpectroscopicOverlap:
    """Compute the spectroscopic overlap between two nuclei."""
    
    def __init__(
            self, nucleus_name, Z, N, kvnn, lamb, channels, jmax=4,
            ws_prm='seminole', kmax=15.0, kmid=3.0, ntot=120,
            generator='Wegner', neval=1e3
    ):
        
        # Class attributes
        # ...
        self.neval = neval
        
        # Set table of Clebsch-Gordan coefficients
        cg_table = compute_clebsch_gordan_table(jmax)
        
        # Woods-Saxon for | \Phi_0^A(\lambda) >
        ws_phi = WoodsSaxon(nucleus_name, Z, N, cg_table,
                            parametrization=ws_prm)
        A = Z + N
        # Woods-Saxon for | \Psi_0^A(\lambda) >
        if A == 4:
            ws_psi = WoodsSaxon('H3', Z-1, N, cg_table, parametrization=ws_prm)
        elif A == 12:
            ws_psi = WoodsSaxon('B11', Z-1, N, cg_table, parametrization=ws_prm)
        elif A == 16:
            ws_psi = WoodsSaxon('N15', Z-1, N, cg_table, parametrization=ws_prm)
        elif A == 40:
            ws_psi = WoodsSaxon('Cl39', Z-1, N, cg_table,
                                parametrization=ws_prm)
        elif A == 56:
            ws_psi = WoodsSaxon('Co55', Z-1, N, cg_table,
                                parametrization=ws_prm)
        
        # Occupied states of | \Phi_0^A(\lambda) >
        self.phi_states = ws_phi.occupied_states
        
        # Occupied states of < \Psi_0^A(\lambda) |
        self.psi_states = []
        i, iZ, iN = 0, 0, 0
        while iZ + iN < A:
            sp_state = ws_psi.sp_states[i]
            if sp_state.m_t == 1/2 and iZ < Z:
                self.psi_states.append(sp_state)
                iZ += 1
            elif sp_state.m_t == -1/2 and iN < N:
                self.psi_states.append(sp_state)
                iN += 1
            i += 1

        # Set Woods-Saxon attributes
        self.dr = ws_phi.dr
        self.ws_phi = ws_phi
        self.ws_psi = ws_psi
        
        # Compute overlap matrix
        self.overlap_matrix = np.zeros((A, A))
        for i, phi in enumerate(self.phi_states):
            for j, psi in enumerate(self.psi_states):
                self.overlap_matrix[i, j] = self.M_ij(phi, psi)
                
        # Inverse of overlap matrix
        self.overlap_matrix_inv = la.inv(self.overlap_matrix)
        
        # Determinant of overlap matrix
        self.overlap_matrix_det = la.det(self.overlap_matrix)
        
        # Initialize \delta U matrix element class
        if lamb == np.inf:
            
            self.flag = False
            
        else:

            self.flag = True
            self.delta_U_matrix_element = DeltaUMatrixElement(
                cg_table, kvnn, kmax, kmid, ntot, generator, lamb, channels
            )
    
    def M_ij(self, phi, psi):
        """Compute the overlap between two s.p. states \phi and \psi."""
        
        # Orthogonality of isospinors and vector spherical harmonics
        cond = (phi.l == psi.l and phi.j == psi.j and phi.m_j == psi.m_j
                and phi.m_t == psi.m_t)
        
        if cond:

            _, u_phi = self.ws_phi.get_wf_rspace(phi)
            _, u_psi = self.ws_psi.get_wf_rspace(psi)
            
            # Overlap is \int_0^\infty dr u_\phi(r) * u_\psi(r)
            return self.dr * np.sum(u_phi * u_psi)
        
        else:
            
            return 0
        
    def alpha_index(self, alpha):
        """Obtain the index of the s.p. state \alpha within the list containing
        occupied states of < \Psi |.
        """
        
        index = 0
        for j, sp_state in enumerate(self.psi_states):
            if sp_state == alpha:
                index = j
                break
                
        return index
    
    def compute_I(self, q_array, alpha, sigma, tau):
        """Compute the I term in the overlap A(q;\sigma,\tau).
        """
        
        # Angles are set to 0 since q_vector = q * \hat{z}
        theta_array = np.zeros_like(q_array)
        phi_array = np.zeros_like(q_array)

        # Index of \alpha in occupied states of < \Psi_0^A(\lambda) |
        jalpha = self.alpha_index(alpha)

        # \beta must have same quantum numbers as \alpha
        ibeta, beta = jalpha, alpha
            
        # Single-particle wave function in momentum space
        phi_dagger_array = np.conj(
            self.ws_phi.psi(beta, q_array, theta_array, phi_array, sigma, tau)
        )
        
        # I term
        I_array = (self.overlap_matrix_det * phi_dagger_array
                   * self.overlap_matrix_inv[ibeta, jalpha])
        
        # Print information
        print("I term completed.")

        return I_array
    
    def compute_delta_U_dagger(self, q_array, alpha, sigma, tau, neval):
        """Compute the \delta U^\dagger term."""
    
        # Get sets of four spin projection configurations
        # (\sigma, \sigma', \sigma_1, \sigma_2)
        spin_configurations = set_spin_configurations(sigma)
    
        # Get sets of four isospin projection configurations
        # (tau, \tau', \tau_1, \tau_2)
        isospin_configurations = set_isospin_configurations(tau)

        # Relative momenta from 0 to 10 fm^-1
        k_limits = [0, 10]
        # C.o.M. momenta up to 3 fm^-1
        K_limits = [0, 3]
        # Polar angle from 0 to \pi
        theta_limits = [0, np.pi]
        # Azimuthal angle from 0 to 2\pi
        phi_limits = [0, 2*np.pi]

        # Evaluate the real part of \delta U^\dagger term for each q
        delta_U_real_array = np.zeros_like(q_array)
        delta_U_real_errors = np.zeros_like(q_array)

        # Set-up integrator with multiple processors
        integ_real = vegas.Integrator([
            k_limits, theta_limits, phi_limits,
            K_limits, theta_limits, phi_limits,
        ], nproc=8)

        t0 = time.time()
        for i, q in enumerate(q_array):
        
            # Initialize integrand given particular q and alpha
            integrand_real = Real(
                q, alpha, self.psi_states, self.overlap_matrix_inv,
                self.overlap_matrix_det, spin_configurations,
                isospin_configurations, self.delta_U_matrix_element,
                self.ws_psi, self.ws_phi
            )

            # Train the integrator
            integ_real(integrand_real, nitn=5, neval=neval)
            
            # Final result
            result = integ_real(integrand_real, nitn=10, neval=neval)
        
            delta_U_real_array[i] = result.mean
            delta_U_real_errors[i] = result.sdev
        
            # percent = (i+1)/(2*len(q_array))*100
            # print(f"{percent:.2f} percent done with \delta U^\dagger.")
            
        # Print information
        t1 = time.time()
        mins = (t1-t0)/60
        print("Real part of \delta U^\dagger term completed after "
              f"{mins:.2f} minutes.")

        # Evaluate the imaginary part of \delta U^\dagger term for each q
        delta_U_imag_array = np.zeros_like(q_array)
        delta_U_imag_errors = np.zeros_like(q_array)

        # Set-up integrator with multiple processors
        integ_imag = vegas.Integrator([
            k_limits, theta_limits, phi_limits,
            K_limits, theta_limits, phi_limits,
        ], nproc=8)

        t0 = time.time()
        for i, q in enumerate(q_array):
        
            # Initialize integrand given particular q and alpha
            integrand_imag = Imag(
                q, alpha, self.psi_states, self.overlap_matrix_inv,
                self.overlap_matrix_det, spin_configurations,
                isospin_configurations, self.delta_U_matrix_element,
                self.ws_psi, self.ws_phi
            )

            # Train the integrator
            integ_imag(integrand_imag, nitn=5, neval=neval)
            
            # Final result
            result = integ_imag(integrand_imag, nitn=10, neval=neval)
        
            delta_U_imag_array[i] = result.mean
            delta_U_imag_errors[i] = result.sdev
        
            # percent = (i+1)/(2*len(q_array))*100 + 50
            # print(f"{percent:.2f} percent done with \delta U^\dagger.")
        
        # Print information
        t1 = time.time()
        mins = (t1-t0)/60
        print("Imaginary part of \delta U^\dagger term completed after "
              f"{mins:.2f} minutes.")

        delta_U_array = delta_U_real_array + 1j * delta_U_imag_array
        delta_U_errors = np.sqrt(delta_U_real_errors ** 2
                                 + delta_U_imag_errors ** 2)

        return delta_U_array, delta_U_errors
    
    def compute_overlap(self, q_array, n, l, j, m_t, sigma, tau):
        """Compute the spectroscopic overlap fixing the z-axis along q_vector
        and averaging over m_j.
        """
        
        overlap_mj_array = np.zeros_like(q_array)
        error_mj_array = np.zeros_like(q_array)
        ipm_mj_array = np.zeros_like(q_array)
        
        # Average over m_j_\alpha
        m_j_array = np.arange(-j, j + 1, 1)
        for m_j in m_j_array:
            
            # Woods-Saxon orbital \alpha with m_j specified
            alpha = SingleParticleState(n, l, j, m_j, m_t)
        
            # Get contribution from I term
            I_term = self.compute_I(q_array, alpha, sigma, tau)
            
            # Evolve operator?
            if self.flag:
                
                # Print information
                print(f"\nStarting m_j = {m_j}.")
                
                # Get contribution from \delta U^\dagger term
                delta_U_dagger_term, errors = self.compute_delta_U_dagger(
                    q_array, alpha, sigma, tau, self.neval
                )
                
                # Print information
                print(f"Completed m_j = {m_j}.")
            
            # No operator evolution
            else:
                
                delta_U_dagger_term = np.zeros_like(I_term, dtype=complex)
                errors = np.zeros_like(q_array)
            
            # Sum up contributions for each m_j
            overlap_mj_array += np.abs(I_term + delta_U_dagger_term) ** 2
            error_mj_array += errors ** 2
            ipm_mj_array += np.abs(I_term) ** 2
            
        # Take square root and divide by 2j+1
        overlap_array = np.sqrt(overlap_mj_array / (2*j+1))
        ipm_array = np.sqrt(ipm_mj_array / (2*j+1))
            
        # Take square root of errors
        error_array = np.sqrt(error_mj_array / (2*j+1))
            
        return overlap_array, error_array, ipm_array

    def spectroscopic_factor(self, j, q_array, q_weights, overlap_array):
        """Compute the spectroscopic factor associated with the overlap."""
        
        # Factor of 2 for \sigma = +1/2 or -1/2
        factor = np.sqrt(2 * 4*np.pi * (2*j+1) * (2*np.pi) ** 3)
        
        # Give overlap the same normalization as VMC overlaps
        overlap_array *= factor

        return (np.sum(q_weights * q_array ** 2 * np.abs(overlap_array) ** 2)
                / (2*np.pi) ** 3 )
    
    def save(
            self, filename, j, q_array, q_weights, overlap_array, error_array,
            ipm_array
    ):
        """Convert A(q;\sigma,\tau) to have the same normalization as VMC
        overlaps and save.
        """
        
        # Factor of 2 for \sigma = +1/2 or -1/2
        factor = np.sqrt(2 * 4*np.pi * (2*j+1) * (2*np.pi) ** 3)
        
        # Give overlap the same normalization as VMC overlaps
        overlap_array *= factor
        error_array *= factor
        ipm_array *= factor
        
        data = np.vstack((q_array, q_weights, overlap_array, error_array,
                          ipm_array)).T
                    
        hdr = ("q, q weight, overlap, error, I term\n")
        
        np.savetxt(filename + '.txt', data, header=hdr)
        
def load(filename):
    """Load overlap and its error."""
        
    data = np.loadtxt(filename + '.txt')
        
    q_array = data[:, 0]
    q_weights = data[:, 1]
    overlap_array = data[:, 2]
    error_array = data[:, 3]
    ipm_array = data[:, 4]
        
    return q_array, q_weights, overlap_array, error_array, ipm_array
    

@njit
def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return int(x == y)


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
                            cg_table[(j_1, m_1, j_2, m_2, j_3, m_3)] = float(
                                CG(j_1, m_1, j_2, m_2, j_3, m_3).doit()
                            )
                            
    return cg_table


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

    # Must specify axis=0 for vegas batchmode
    k = la.norm(k_vector, axis=0)
    theta = np.arccos(k_vector[2]/k)
    phi = np.arctan2(k_vector[1], k_vector[0])

    return k, theta, phi


def set_spin_configurations(sigma):
    """Spin projection configurations for \delta U^\dagger term."""
    
    spin_projections = [1/2, -1/2]
    spin_configurations = []
    for sigmap in spin_projections:
        for sigma_1 in spin_projections:
            for sigma_2 in spin_projections:
                spin_configurations.append((sigma, sigmap, sigma_1, sigma_2))
                        
    return spin_configurations


def set_isospin_configurations(tau):
    """Isospin projection configurations for \delta U^\dagger term."""
    
    isospin_projections = [1/2, -1/2]
    isospin_configurations = []
    for taup in isospin_projections:
        for tau_1 in isospin_projections:
            for tau_2 in isospin_projections:
                # Check that M_T is conserved
                if tau + taup == tau_1 + tau_2:
                    isospin_configurations.append((tau, taup, tau_1, tau_2))
    
    return isospin_configurations


if __name__ == '__main__':
    
    # Nucleus
    # nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'C12', 6, 6
    # nucleus_name, Z, N = 'O16', 8, 8
    # nucleus_name, Z, N = 'Ar40', 18, 22
    nucleus_name, Z, N = 'Ni56', 28, 28

    # Quantum state
    # n, l, j, m_t = 1, 0, 1/2, 1/2  # 1s_{1/2}
    # n, l, j, m_t = 1, 1, 3/2, 1/2  # 1p_{3/2}
    # n, l, j, m_t = 1, 1, 1/2, 1/2  # 1p_{1/2}
    # n, l, j, m_t = 2, 0, 1/2, 1/2  # 2s_{1/2}
    n, l, j, m_t = 1, 3, 7/2, 1/2  # 1f_{7/2}

    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    
    # NN potential and momentum mesh
    # kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120  # AV18
    # kvnn, kmax, kmid, ntot = 7, 15.0, 3.0, 120  # CD-Bonn
    # kvnn, kmax, kmid, ntot = 79, 15.0, 3.0, 120  # EMN N4LO 500 MeV
    kvnn, kmax, kmid, ntot = 113, 15.0, 3.0, 120  # SMS N4LO 550 MeV
    
    # SRG \lambda value
    lamb = 1.5
    # lamb = 3.0
    # lamb = np.inf
    
    # Number of evaluations for vegas
    neval = 1e3
    # neval = 5e3
    
    # Initialize overlap class
    so = SpectroscopicOverlap(nucleus_name, Z, N, kvnn, lamb, channels,
                              neval=neval)

    # Compute overlap
    q_array, q_weights = momentum_mesh(10.0, 2.0, 120)
    overlap_array, error_array, ipm_array = so.compute_overlap(q_array, n, l, j,
                                                                m_t, 1/2, 1/2)
    
    # Save overlap
    if m_t == 1/2:
        nucleon = 'proton'
    elif m_t == -1/2:
        nucleon = 'neutron'
    filename = (replace_periods(f"{nucleus_name}_{nucleon}_n{n}_l{l}_j{2*j}"
                                f"_overlap_kvnn_{kvnn}_lamb_{lamb}") + '.txt')
    so.save(filename, j, q_array, q_weights, overlap_array, error_array,
            ipm_array)