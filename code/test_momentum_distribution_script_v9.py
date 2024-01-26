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
plane-wave matrix elements of \delta U, integrating over products of the matrix
elements and s.p. wavefunctions \psi. The MC package vegas is used to integrate
over vector momenta and s.p. quantum numbers \alpha, \beta, and spin
projections.

Last update: November 10, 2023

"""

# Python imports
from numba import njit
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from scipy.special import sph_harm, spherical_jn
from sympy.physics.quantum.cg import CG
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
    
    ### TESTING
    def __init__(self, nucleus_name, Z, N, cg_table, rmax=40, ntab=2000,
                 kmax=10.0, kmid=2.0, ntot=120, parametrization='Match'):
        
        # Set instance attributes
        # self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}/"
        ### TESTING
        self.woods_saxon_directory = (
            f"../data/woods_saxon/{parametrization}/{nucleus_name}/"
        )
        # if test:
        #     # self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}_test/"
        #     # self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}_seminole/"
        #     self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}_universal/"
        # else:
        #     self.woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}/"
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
            
            return 0+0j
        
        
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
    
    
class DeltaUMatrixElement:
    """DOCSTRING"""
    
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
        
        matrix_element = 0+0j
        
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
    
    
class DeltaUIntegrand:
    """Evaluates the integrand of the \delta U + \delta U^\dagger term."""
    
    
    def __init__(
            self, q, tau, woods_saxon, delta_U_matrix_element,
            spin_configurations, isospin_configurations, occupied_state_pairs
    ):
        
        # Set momenta and isospin as instance attributes
        self.q = q
        self.tau = tau
        self.isospin_configurations = isospin_configurations
        
        # Set Woods-Saxon and \delta U matrix element classes as instance attributes
        self.woods_saxon = woods_saxon
        self.delta_U_matrix_element = delta_U_matrix_element
        
        # Possible pairs of spin projections
        self.spin_configurations = spin_configurations
        self.N_spin = len(spin_configurations)
    
        # Possible pairs of occupied s.p. states
        self.occupied_state_pairs = occupied_state_pairs
        self.N_occupied_states = len(occupied_state_pairs)
    
    
    def __call__(self, x_array):
        """This method assumes sampling over the following variables:
        
            k, theta_k, phi_k -> vector momenta k,
            K, theta_K, phi_K -> vector momenta K,
            spin projections -> sigma_1, sigma_2, sigma, sigma',
            s.p. state pairs -> \alpha, \beta
            
        where x_array contains each variable meaning x_array.shape = (8,).
        """
        
        # Choose z-axis to be along q_vector
        q_vector = np.array([0, 0, self.q])
        
        # Relative momenta k
        k, theta_k, phi_k = x_array[0:3]
        k_vector = build_vector(k, theta_k, phi_k)
        
        # C.o.M. momenta K
        K, theta_K, phi_K = x_array[3:6]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Spin projections \sigma_1, \sigma_2, \sigma, \sigma'
        sigma_1, sigma_2, sigma, sigmap = self.get_spin_configuration(
            x_array[6])
        
        # S.p. state pairs \alpha, \beta < F
        alpha, beta = self.get_occupied_state_pair(x_array[7])
        
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
        
        # \psi_\alpha(q; \sigma, \tau)
        psi_alpha_q = self.woods_saxon.psi(alpha, self.q, 0, 0, sigma, self.tau)
        
        # \psi_\beta(q; \sigma, \tau)
        psi_beta_q = self.woods_saxon.psi(beta, self.q, 0, 0, sigma, self.tau)
        
        # Sum over isospin projections
        integrand = 0+0j
        for isospin_projections in self.isospin_configurations:
            
            taup, tau_1, tau_2 = isospin_projections
                
            # Plane-wave matrix elements of \delta U and \delta U^\dagger
            delta_U_plane_wave = self.delta_U_matrix_element(
                k, theta_k, phi_k, qK, theta_qK, phi_qK, sigma_1, sigma_2,
                sigma, sigmap, tau_1, tau_2, self.tau, taup, hc=False
            )
            delta_U_dag_plane_wave = self.delta_U_matrix_element(
                qK, theta_qK, phi_qK, k, theta_k, phi_k, sigma, sigmap, sigma_1,
                sigma_2, self.tau, taup, tau_1, tau_2, hc=True
            )
        
            # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
            psi_alpha_1 = self.woods_saxon.psi(alpha, k1, theta_k1, phi_k1,
                                               sigma_1, tau_1)
        
            # \psi_\beta(K/2-k; \sigma_2, \tau_2)
            psi_beta_2 = self.woods_saxon.psi(beta, k2, theta_k2, phi_k2,
                                              sigma_2, tau_2)

            # \psi_\beta(K-q; \sigma', \tau')
            psi_beta_Kq = self.woods_saxon.psi(beta, Kq, theta_Kq, phi_Kq,
                                               sigmap, taup)
        
            # \psi_\alpha(K-q; \sigma', \tau')
            psi_alpha_Kq = self.woods_saxon.psi(alpha, Kq, theta_Kq, phi_Kq,
                                                sigmap, taup)
        
            # \delta U term
            delta_U_term = (
                delta_U_plane_wave * np.conj(psi_alpha_1) * np.conj(psi_beta_2)
                * (psi_beta_Kq * psi_alpha_q - psi_alpha_Kq * psi_beta_q)
            )
        
            # \delta U^\dagger term
            delta_U_dag_term = (
                delta_U_dag_plane_wave * psi_alpha_1 * psi_beta_2
                * (
                    np.conj(psi_beta_Kq) * np.conj(psi_alpha_q)
                    - np.conj(psi_alpha_Kq) * np.conj(psi_beta_q)
                )
            )

            # Add together for full integrand
            integrand += (
                1/2 * jacobian * self.N_spin * self.N_occupied_states * (
                    delta_U_term + delta_U_dag_term
                )
            )
        
        return integrand.real    
    
    
    def get_spin_configuration(self, x):
        """Given a number between 0 and 1, return a set of four spin
        projections.
        """
        
        index = np.floor(x * self.N_spin).astype(int)
        return self.spin_configurations[index]
    
    
    def get_occupied_state_pair(self, x):
        """Given a number between 0 and 1, return a pair of occupied quantum
        states.
        """
        
        index = np.floor(x * self.N_occupied_states).astype(int)
        return self.occupied_state_pairs[index]
    
    
class DeltaU2Integrand:
    """Evaluates the integrand of the \delta U \delta U^\dagger term."""
    
    
    def __init__(
            self, q, tau, woods_saxon, delta_U_matrix_element,
            spin_configurations, isospin_configurations, occupied_state_pairs
    ):
        
        # Set momenta and isospin as instance attributes
        self.q = q
        self.tau = tau
        self.isospin_configurations = isospin_configurations
        
        # Set Woods-Saxon and \delta U matrix element classes as instance attributes
        self.woods_saxon = woods_saxon
        self.delta_U_matrix_element = delta_U_matrix_element
        
        # Possible pairs of spin projections
        self.spin_configurations = spin_configurations
        self.N_spin = len(spin_configurations)
    
        # Possible pairs of occupied s.p. states
        self.occupied_state_pairs = occupied_state_pairs
        self.N_occupied_states = len(occupied_state_pairs)
    
    
    def __call__(self, x_array):
        """This method assumes sampling over the following variables:
        
            k, theta_k, phi_k -> vector momenta k,
            kp, theta_kp, phi_kp -> vector momenta k',
            K, theta_K, phi_K -> vector momenta K,
            spin projections -> sigma_1, sigma_2, sigma, sigma', sigma_3,
                sigma_4
            s.p. state pairs -> \alpha, \beta
            
        where x_array contains each variable meaning x_array.shape = (11,).
        """
        
        # Choose z-axis to be along q_vector
        q_vector = np.array([0, 0, self.q])
        
        # Relative momenta k
        k, theta_k, phi_k = x_array[0:3]
        k_vector = build_vector(k, theta_k, phi_k)
        
        # Relative momenta k'
        kp, theta_kp, phi_kp = x_array[3:6]
        kp_vector = build_vector(kp, theta_kp, phi_kp)
        
        # C.o.M. momenta K
        K, theta_K, phi_K = x_array[6:9]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Spin projections \sigma_1, \sigma_2, \sigma, \sigma'
        sigma_1, sigma_2, sigma, sigmap, sigma_3, sigma_4 = (
            self.get_spin_configuration(x_array[9]))
        
        # S.p. state pairs \alpha, \beta < F
        alpha, beta = self.get_occupied_state_pair(x_array[10])
        
        # Calculate vector q - K/2
        qK_vector = q_vector - K_vector/2
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
        
        # Calculate vector k_1 = K/2 + k
        k1_vector = K_vector/2 + k_vector
        k1, theta_k1, phi_k1 = get_vector_components(k1_vector)
        
        # Calculate vector k_2 = K/2 - k
        k2_vector = K_vector/2 - k_vector
        k2, theta_k2, phi_k2 = get_vector_components(k2_vector)
        
        # Calculate vector k_3 = K/2 + k'
        k3_vector = K_vector/2 + kp_vector
        k3, theta_k3, phi_k3 = get_vector_components(k3_vector)
        
        # Calculate vector k_4 = K/2 - k'
        k4_vector = K_vector/2 - kp_vector
        k4, theta_k4, phi_k4 = get_vector_components(k4_vector)
        
        # Calculate the Jacobian determinant
        jacobian = (k ** 2 * np.sin(theta_k) * kp ** 2 * np.sin(theta_kp)
                    * K ** 2 * np.sin(theta_K))
        
        # Sum over isospin projections
        integrand = 0+0j
        for isospin_projections in self.isospin_configurations:
            
            taup, tau_1, tau_2, tau_3, tau_4 = isospin_projections
                
            # Plane-wave matrix elements of \delta U and \delta U^\dagger
            delta_U_plane_wave = self.delta_U_matrix_element(
                k, theta_k, phi_k, qK, theta_qK, phi_qK, sigma_1, sigma_2,
                sigma, sigmap, tau_1, tau_2, self.tau, taup, hc=False
            )
            delta_U_dag_plane_wave = self.delta_U_matrix_element(
                qK, theta_qK, phi_qK, kp, theta_kp, phi_kp, sigma, sigmap,
                sigma_3, sigma_4, self.tau, taup, tau_3, tau_4, hc=True
            )
        
            # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
            psi_alpha_1 = self.woods_saxon.psi(alpha, k1, theta_k1, phi_k1,
                                               sigma_1, tau_1)
        
            # \psi_\beta(K/2-k; \sigma_2, \tau_2)
            psi_beta_2 = self.woods_saxon.psi(beta, k2, theta_k2, phi_k2,
                                              sigma_2, tau_2)
            
            # \psi_\alpha(K/2+kp; \sigma_3, \tau_3)
            psi_alpha_3 = self.woods_saxon.psi(alpha, k3, theta_k3, phi_k3,
                                               sigma_3, tau_3)
        
            # \psi_\beta(K/2-kp; \sigma_4, \tau_4)
            psi_beta_4 = self.woods_saxon.psi(beta, k4, theta_k4, phi_k4,
                                              sigma_4, tau_4)
            
            # \psi_\alpha(K/2-kp; \sigma_4, \tau_4)
            psi_alpha_4 = self.woods_saxon.psi(alpha, k4, theta_k4, phi_k4,
                                               sigma_4, tau_4)

            # \psi_\beta(K/2+kp; \sigma_3, \tau_3)
            psi_beta_3 = self.woods_saxon.psi(beta, k3, theta_k3, phi_k3,
                                              sigma_3, tau_3)

            # Add together for full integrand
            integrand += (
                1/4 * jacobian * self.N_spin * self.N_occupied_states
                * delta_U_plane_wave * delta_U_dag_plane_wave
                * np.conj(psi_alpha_1) * np.conj(psi_beta_2) * (
                    psi_alpha_3 * psi_beta_4 - psi_alpha_4 * psi_beta_3
                )
            )
        
        return integrand.real    
    
    
    def get_spin_configuration(self, x):
        """Given a number between 0 and 1, return a set of six spin
        projections.
        """
        
        index = np.floor(x * self.N_spin).astype(int)
        return self.spin_configurations[index]
    
    
    def get_occupied_state_pair(self, x):
        """Given a number between 0 and 1, return a pair of occupied quantum
        states.
        """
        
        index = np.floor(x * self.N_occupied_states).astype(int)
        return self.occupied_state_pairs[index]
    
    
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
                            
    print(f"Done calculating Clebsch-Gordan table up to j_max = {j_max}.")
                                
    return cg_table


@njit
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


@njit
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

    k = norm(k_vector)
    theta = np.arccos(k_vector[2]/k)
    phi = np.arctan2(k_vector[1], k_vector[0])

    return k, theta, phi


def set_occupied_state_pairs(woods_saxon):
    """List of occupied s.p. state pairs."""
    
    occupied_state_pairs = []
    for alpha in woods_saxon.occupied_states:
        for beta in woods_saxon.occupied_states:
            occupied_state_pairs.append((alpha, beta))
    return occupied_state_pairs


def set_delU_spin_configurations():
    """Spin projection configurations for \delta U + \delta U^\dagger term."""
    
    spin_projections = [1/2, -1/2]
    spin_configurations = []
    for spin_1 in spin_projections:
        for spin_2 in spin_projections:
            for spin_3 in spin_projections:
                for spin_4 in spin_projections:
                    spin_configurations.append((spin_1, spin_2, spin_3, spin_4))
                        
    return spin_configurations


def set_delU2_spin_configurations():
    """Spin projection configurations for \delta U \delta U^\dagger term."""
    
    spin_projections = [1/2, -1/2]
    spin_configurations = []
    for spin_1 in spin_projections:
        for spin_2 in spin_projections:
            for spin_3 in spin_projections:
                for spin_4 in spin_projections:
                    for spin_5 in spin_projections:
                        for spin_6 in spin_projections:
                            spin_configurations.append((spin_1, spin_2, spin_3,
                                                        spin_4, spin_5, spin_6))
                        
    return spin_configurations


def set_delU_isospin_configurations(tau):
    """Isospin projection configurations for \delta U + \delta U^\dagger term.
    """
    
    isospin_projections = [1/2, -1/2]
    isospin_configurations = []
    for taup in isospin_projections:
        for tau_1 in isospin_projections:
            for tau_2 in isospin_projections:
                # Check that M_T is conserved
                if tau + taup == tau_1 + tau_2:
                    isospin_configurations.append((taup, tau_1, tau_2))
    
    return isospin_configurations


def set_delU2_isospin_configurations(tau):
    """Isospin projection configurations for \delta U + \delta U^\dagger term.
    """
    
    isospin_projections = [1/2, -1/2]
    isospin_configurations = []
    for taup in isospin_projections:
        for tau_1 in isospin_projections:
            for tau_2 in isospin_projections:
                for tau_3 in isospin_projections:
                    for tau_4 in isospin_projections:
                        # Check that M_T is conserved
                        if tau + taup == tau_1 + tau_2 == tau_3 + tau_4:
                            isospin_configurations.append(
                                (taup, tau_1, tau_2, tau_3, tau_4)
                            )
                            
    return isospin_configurations


def compute_I_term(q_array, tau, woods_saxon):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
    theta_array = np.zeros_like(q_array)
    phi_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in woods_saxon.occupied_states:

            # Single-particle wave function with z-axis along q_vector
            psi_alpha_array = woods_saxon.psi(alpha, q_array, theta_array,
                                              phi_array, sigma, tau)

            I_array += np.abs(psi_alpha_array) ** 2
            
    print("Done with I term.")
                    
    return I_array


def compute_delta_U_term(q_array, tau, woods_saxon, delta_U_matrix_element,
                         occupied_state_pairs, neval):
    """Compute the \delta U * n(q) * I + I * n(q) * \delta U^\dagger terms."""
    
    # Get sets of four spin projection configurations
    # (\sigma_1, \sigma_2, \sigma, \sigma')
    spin_configurations = set_delU_spin_configurations()
    
    # Get sets of isospin projection configurations
    # ( \tau', \tau_1, \tau_2) where \tau is fixed
    isospin_configurations = set_delU_isospin_configurations(tau)

    # Relative momenta from 0 to 10 fm^-1
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    # Polar angle from 0 to \pi
    theta_limits = [0, np.pi]
    # Azimuthal angle from 0 to 2\pi
    phi_limits = [0, 2*np.pi]
    # Limits for sums over \alpha and \beta pairs and spin projections
    sum_limits = [0, 1]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([
        k_limits, theta_limits, phi_limits,
        K_limits, theta_limits, phi_limits,
        sum_limits, sum_limits
    ], nproc=8)

    # Evaluate the \delta U + \delta U^\dagger term for each q
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
    for i, q in enumerate(q_array):

        integrand = DeltaUIntegrand(
            q, tau, woods_saxon, delta_U_matrix_element, spin_configurations,
            isospin_configurations, occupied_state_pairs
        )
        
        # Train the integrator
        integ(integrand, nitn=5, neval=neval)
    
        # Final result
        result = integ(integrand, nitn=10, neval=neval)
        
        delta_U_array[i] = result.mean
        delta_U_errors[i] = result.sdev
        
        percent = (i+1)/len(q_array)*100
        print(f"{percent:.2f} percent done with \delta U + \delta U^\dagger.")

    return delta_U_array, delta_U_errors


def compute_delta_U2_term(q_array, tau, woods_saxon, delta_U_matrix_element,
                          occupied_state_pairs, neval):
    """Compute the \delta U * n(q) \delta U^\dagger term."""
    
    # Get sets of six spin projection configurations
    # (\sigma_1, \sigma_2, \sigma, \sigma', \sigma_3, \sigma_4)
    spin_configurations = set_delU2_spin_configurations()
    
    # Get sets of isospin projection configurations
    # (\tau', \tau_1, \tau_2, \tau_3, \tau_4)
    isospin_configurations = set_delU2_isospin_configurations(tau)

    # Relative momenta from 0 to 10 fm^-1
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    # Polar angle from 0 to \pi
    theta_limits = [0, np.pi]
    # Azimuthal angle from 0 to 2\pi
    phi_limits = [0, 2*np.pi]
    # Limits for sums over \alpha and \beta pairs and spin projections
    sum_limits = [0, 1]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([
        k_limits, theta_limits, phi_limits,
        k_limits, theta_limits, phi_limits,
        K_limits, theta_limits, phi_limits,
        sum_limits, sum_limits
    ], nproc=8)

    # Evaluate the \delta U \delta U^\dagger term for each q
    delta_U2_array = np.zeros_like(q_array)
    delta_U2_errors = np.zeros_like(q_array)
    for i, q in enumerate(q_array):

        integrand = DeltaU2Integrand(
            q, tau, woods_saxon, delta_U_matrix_element, spin_configurations,
            isospin_configurations, occupied_state_pairs
        )
        
        # Train the integrator
        integ(integrand, nitn=5, neval=neval)
    
        # Final result
        result = integ(integrand, nitn=10, neval=neval)
        
        delta_U2_array[i] = result.mean
        delta_U2_errors[i] = result.sdev
        
        percent = (i+1)/len(q_array)*100
        print(f"{percent:.2f} percent done with \delta U \delta U^\dagger.")

    return delta_U2_array, delta_U2_errors


def compute_normalization(q_array, q_weights, n_array):
    """Compute the normalization of the momentum distribution."""

    return 4 * np.pi * np.sum(q_weights * q_array ** 2 * n_array)


def compute_momentum_distribution(
        nucleus_name, Z, N, tau, kvnn, lamb, channels, kmax=15.0, kmid=3.0,
        ntot=120, generator='Wegner', neval=5e4, print_normalization=False,
        kvnn_hard=None, lambda_m=None, parametrization='Match', ipm=False,
        save=False
):
    """Compute the single-nucleon momentum distribution."""
    
    # Set table of Clebsch-Gordan coefficients
    jmax = 4  # This should cover nuclei as heavy as Ca48
    cg_table = compute_clebsch_gordan_table(jmax)
    
    # Set single-particle basis
    # woods_saxon = WoodsSaxon(nucleus_name, Z, N, cg_table)
    ### TESTING
    woods_saxon = WoodsSaxon(nucleus_name, Z, N, cg_table,
                             parametrization=parametrization)

    # Set points in q
    q_array, q_weights = momentum_mesh(10.0, 2.0, 120)
    
    # Compute the I term
    I_array = compute_I_term(q_array, tau, woods_saxon)
    
    # IPM only?
    if ipm:
        
        delta_U_array = np.zeros_like(q_array)
        delta_U_errors = np.zeros_like(q_array)
        delta_U2_array = np.zeros_like(q_array)
        delta_U2_errors = np.zeros_like(q_array)
    
    # Full distribution
    else:
        
        # Get pairs of occupied states
        occupied_state_pairs = set_occupied_state_pairs(woods_saxon)
    
        # Initialize \delta U matrix element class
        delta_U_matrix_element = DeltaUMatrixElement(
            cg_table, kvnn, kmax, kmid, ntot, generator, lamb, channels,
            kvnn_hard, lambda_m
        )
    
        # Compute the \delta U + \delta U^\dagger term
        delta_U_array, delta_U_errors = compute_delta_U_term(
            q_array, tau, woods_saxon, delta_U_matrix_element,
            occupied_state_pairs, neval
        )
    
        # Compute the \delta U \delta U^\dagger term
        delta_U2_array, delta_U2_errors = compute_delta_U2_term(
            q_array, tau, woods_saxon, delta_U_matrix_element,
            occupied_state_pairs, neval
        )
    
    # Combine each term for the total momentum distribution [fm^3]
    n_array = I_array + delta_U_array + delta_U2_array
    n_errors = np.sqrt(delta_U_errors ** 2 + delta_U2_errors ** 2)

    # Option to print normalization of the total momentum distribution
    if print_normalization:
        normalization = compute_normalization(q_array, q_weights, n_array)
        print(f"Normalization = {normalization:.5f}.")
    
    # Option to save the momentum distribution as a .txt file
    if save:
        save_momentum_distribution(
            nucleus_name, tau, kvnn, lamb, q_array, q_weights, n_array,
            n_errors, I_array, delta_U_array, delta_U_errors, delta_U2_array,
            delta_U2_errors, kvnn_hard, lambda_m
        )
    
    return q_array, q_weights, n_array, n_errors


def save_momentum_distribution(
        nucleus_name, tau, kvnn, lamb, q_array, q_weights, n_array, n_errors,
        I_array, delta_U_array, delta_U_errors, delta_U2_array, delta_U2_errors,
        kvnn_hard=None, lambda_m=None, parametrization='Match'
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
    
    directory = f"../data/momentum_distributions/{nucleus_name}/"

    ### TESTING
    if kvnn_hard is not None:
        file_name = replace_periods(
            f"{nucleus_name}_{nucleon}_momentum_distribution_kvnn_{kvnn}_lamb"
            f"_{lamb}_kvnn_hard_{kvnn_hard}_lambda_m_{lambda_m}"
            f"_{parametrization}"
        )
    else:
        file_name = replace_periods(
            f"{nucleus_name}_{nucleon}_momentum_distribution_kvnn_{kvnn}_lamb"
            f"_{lamb}_{parametrization}"
        )
    
    np.savetxt(directory + file_name + '.txt', data, header=hdr)


def load_momentum_distribution(
        nucleus_name, nucleon, kvnn, lamb, kvnn_hard=None, lambda_m=None
):
    """Load and return the momentum distribution along with the isolated
    contributions.
    """
    
    directory = f"../data/momentum_distributions/{nucleus_name}/"

    if kvnn_hard is not None:
        file_name = replace_periods(
            f"{nucleus_name}_{nucleon}_momentum_distribution_kvnn_{kvnn}_lamb"
            f"_{lamb}_kvnn_hard_{kvnn_hard}_lambda_m_{lambda_m}"
        )
    else:
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
    nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'C12', 6, 6
    # nucleus_name, Z, N = 'O16', 8, 8
    # nucleus_name, Z, N = 'Ca40', 20, 20
    # nucleus_name, Z, N = 'Ca48', 20, 28
    # nucleus_name, Z, N = 'Pb208', 82, 126
    
    # Nucleon
    tau = 1/2
    
    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    
    # NN potential and momentum mesh
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120  # AV18
    # kvnn, kmax, kmid, ntot = 7, 15.0, 3.0, 120  # CD-Bonn
    # kvnn, kmax, kmid, ntot = 79, 15.0, 3.0, 120  # EMN N4LO 500 MeV
    # kvnn, kmax, kmid, ntot = 113, 15.0, 3.0, 120  # SMS N4LO 550 MeV
    
    # SRG \lambda value
    # lamb = 1.35
    lamb = 1.5
    # lamb = 2.0
    # lamb = 2.5
    
    neval = 5e4  # 4He
    # neval = 7.5e4  # 12C
    # neval = 1e5  # 16O
    # neval = 5e5  # 40Ca and 48Ca
    
    # Inverse-SRG evolution?
    kvnn_hard = None
    lambda_m = None
    # kvnn_hard = 6
    # lambda_m = 5.5
    # lambda_m = 5.0
    # lambda_m = 4.5
    # lambda_m = 4.0
    
    # Woods-Saxon parametrization TESTING
    # prm = 'Seminole'
    # prm = 'Universal'
    prm = 'Match'

    # Compute and save the momentum distribution TESTING
    q_array, q_weights, n_array, n_errors = compute_momentum_distribution(
        nucleus_name, Z, N, tau, kvnn, lamb, channels, neval=neval,
        kvnn_hard=kvnn_hard, lambda_m=lambda_m, parametrization=prm, save=True
    )
