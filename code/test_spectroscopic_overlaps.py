#!/usr/bin/env python3

"""
File: test_spectroscopic_overlaps.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 14, 2024

This script serves as a testbed for calculating spectroscopic overlaps using
mean field approximations for initial and final states and applying SRG
transformations to the operator.

Last update: February 23, 2024

"""

# Python imports
from numba import njit
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from scipy.special import sph_harm, spherical_jn
from sympy.physics.quantum.cg import CG

# Imports from scripts
from scripts.integration import (
    gaussian_quadrature_mesh, momentum_mesh, unattach_weights_from_matrix
)
from scripts.potentials import Potential
from scripts.srg import compute_srg_transformation, load_srg_transformation
from scripts.tools import convert_l_to_string, coupled_channel, replace_periods


#### TODO: Update sum over m_j
#### TODO: Use JAX instead of NumPy?


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
                 kmax=10.0, kmid=2.0, ntot=120, parametrization='match'):
        
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
            
            # return 0+0j
            ### TESTING
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


class DeltaUDagger:
    """Evaluates the \delta U^\dagger term."""
    
    
    def __init__(
            self, alpha, woods_saxon, delta_U_matrix_element,
            spin_configurations, isospin_configurations
    ):
        
        # Set momenta and isospin as instance attributes
        self.alpha = alpha
        self.spin_configurations = spin_configurations
        self.isospin_configurations = isospin_configurations
        
        # Set Woods-Saxon and \delta U matrix element classes as instance
        # attributes
        self.woods_saxon = woods_saxon
        self.delta_U_matrix_element = delta_U_matrix_element
        
        # Set momenta and angles
        kmax, ktot = 10.0, 100
        k_array, k_weights = gaussian_quadrature_mesh(kmax, ktot)
        Kmax, Ktot = 5.0, 50
        K_array, K_weights = gaussian_quadrature_mesh(Kmax, Ktot)
        theta_array, theta_weights = gaussian_quadrature_mesh(np.pi, 11)
        phi_array, phi_weights = gaussian_quadrature_mesh(2*np.pi, 15)
        
        # Create meshgrids
        (self.k_grid, self.thetak_grid, self.phik_grid, self.K_grid,
         self.thetaK_grid, self.phiK_grid) = np.meshgrid(
             k_array, theta_array, phi_array, K_array, theta_array, phi_array,
             indexing='ij'
        )
        (self.dk_grid, self.dthetak_grid, self.dphik_grid, self.dK_grid,
         self.dthetaK_grid, self.dphiK_grid) = np.meshgrid(
            k_weights, theta_weights, phi_weights, K_weights, theta_weights,
            phi_weights, indexing='ij'
        )
             
        # Jacobian determinant
        self.jacobian = (
            self.k_grid ** 2 * self.dk_grid * np.sin(self.thetak_grid)
            * self.dthetak_grid * self.dphik_grid * self.K_grid ** 2
            * self.dK_grid * np.sin(self.thetaK_grid) * self.dthetaK_grid
            * self.dphiK_grid
        )
             
        
    def __call__(self, q):
        """Integrate over K and k, and sum over quantum numbers."""
        
        # Calculate vector q - K/2
        qK, theta_qK, phi_qK = get_vector_components(
            q, 0.0, 0.0, -self.K_grid/2, self.thetaK_grid, self.phiK_grid
        )
        
        # Calculate vector K - q
        Kq, theta_Kq, phi_Kq = get_vector_components(
            self.K_grid, self.thetaK_grid, self.phiK_grid, -q, 0.0, 0.0
        )
        
        # Calculate vector k_1 = K/2 + k
        k1, theta_k1, phi_k1 = get_vector_components(
            self.K_grid/2, self.thetaK_grid, self.phiK_grid, self.k_grid,
            self.thetak_grid, self.phik_grid
        )
        
        # Calculate vector k_2 = K/2 - k
        k2, theta_k2, phi_k2 = get_vector_components(
            self.K_grid/2, self.thetaK_grid, self.phiK_grid, -self.k_grid,
            self.thetak_grid, self.phik_grid
        )
        
        # Sum over spin projections
        integrand = np.zeros_like(self.k_grid, dtype=complex)
        for spin_projections in self.spin_configurations:
                    
            sigma, sigmap, sigma_1, sigma_2 = spin_projections
                    
            # Sum over isospin projections
            for isospin_projections in self.isospin_configurations:
                        
                tau, taup, tau_1, tau_2 = isospin_projections
                
                # Plane-wave matrix elements of \delta U^\dagger
                delta_U_dag_plane_wave = self.delta_U_matrix_element(
                    qK, theta_qK, phi_qK, self.k_grid, self.thetak_grid,
                    self.phik_grid, sigma, sigmap, sigma_1, sigma_2, tau, taup,
                    tau_1, tau_2, hc=True
                )
                
                # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
                psi_alpha_1 = self.woods_saxon.psi(self.alpha, k1, theta_k1,
                                                   phi_k1, sigma_1, tau_1)
                
                # \psi_\alpha(K/2-k; \sigma_2, \tau_2)
                psi_alpha_2 = self.woods_saxon.psi(self.alpha, k2, theta_k2,
                                                   phi_k2, sigma_2, tau_2)
                        
                # Sum over occupied states
                for beta in self.woods_saxon.occupied_states:
                    
                    # \psi_\beta(K-q; \sigma', \tau')
                    psi_beta_Kq = self.woods_saxon.psi(beta, Kq, theta_Kq,
                                                       phi_Kq, sigmap, taup)

                    # \psi_\beta(K/2-k; \sigma_2, \tau_2)
                    psi_beta_2 = self.woods_saxon.psi(beta, k2, theta_k2,
                                                      phi_k2, sigma_2, tau_2)

                    # \psi_\beta(K/2+k; \sigma_1, \tau_1)
                    psi_beta_1 = self.woods_saxon.psi(beta, k1, theta_k1,
                                                      phi_k1, sigma_1, tau_1)

                    # Add together for full integrand
                    integrand += 1/2 * (
                        delta_U_dag_plane_wave * np.conj(psi_beta_Kq)
                        * (psi_beta_2 * psi_alpha_1 - psi_beta_1 * psi_alpha_2)
                    )
          
        # Integrate over K and k
        return np.sum(self.jacobian * integrand)


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
def get_vector_components(x, theta_x, phi_x, y, theta_y, phi_y):
    """Output the spherical components of a vector v = x + y."""
    
    # Components of x vector
    x_x = x * np.sin(theta_x) * np.cos(phi_x)
    x_y = x * np.sin(theta_x) * np.sin(phi_x)
    x_z = x * np.cos(theta_x)
    
    # Components of y vector
    y_x = y * np.sin(theta_y) * np.cos(phi_y)
    y_y = y * np.sin(theta_y) * np.sin(phi_y)
    y_z = y * np.cos(theta_y)
    
    # Dot product of x and y
    xy = x_x * y_x + x_y * y_y + x_z * y_z
    
    # Vector norm
    z = np.sqrt(x ** 2 + y ** 2 + 2 * xy)
    
    # Polar angle
    theta_z = np.arccos((x_z + y_z) / z)
    
    # Azimuthal angle
    phi_z = np.arctan2(x_y + y_y, x_x + y_x)

    return z, theta_z, phi_z


def set_spin_configurations():
    """Spin projection configurations for \delta U^\dagger term."""
    
    spin_projections = [1/2, -1/2]
    spin_configurations = []
    for spin_1 in spin_projections:
        for spin_2 in spin_projections:
            for spin_3 in spin_projections:
                for spin_4 in spin_projections:
                    spin_configurations.append((spin_1, spin_2, spin_3, spin_4))
                        
    return spin_configurations


def set_isospin_configurations():
    """Isospin projection configurations for \delta U^\dagger term."""
    
    isospin_projections = [1/2, -1/2]
    isospin_configurations = []
    for tau in isospin_projections:
        for taup in isospin_projections:
            for tau_1 in isospin_projections:
                for tau_2 in isospin_projections:
                    # Check that M_T is conserved
                    if tau + taup == tau_1 + tau_2:
                        isospin_configurations.append((tau, taup, tau_1, tau_2))
    
    return isospin_configurations


def compute_I_term(alpha, q_array, woods_saxon):
    """Compute the I term."""
        
    I_array = np.zeros_like(q_array, dtype=complex)
    theta_array = np.zeros_like(q_array)
    phi_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
        
        # Loop over isospin projections
        for tau in np.array([1/2, -1/2]):

            # Single-particle wave function with z-axis along q_vector
            psi_alpha_array = woods_saxon.psi(alpha, q_array, theta_array,
                                              phi_array, sigma, tau)
            
            I_array += psi_alpha_array

    print("Done with I term.")
                    
    return I_array


def compute_delta_U_dagger_term(
        alpha, q_array, woods_saxon, delta_U_matrix_element, neval
):
    """Compute the \delta U^\dagger term."""
    
    # Get sets of four spin projection configurations
    # (\sigma, \sigma', \sigma_1, \sigma_2)
    spin_configurations = set_spin_configurations()
    
    # Get sets of four isospin projection configurations
    # (tau, \tau', \tau_1, \tau_2)
    isospin_configurations = set_isospin_configurations()
    
    # Initialize \delta U^\dagger term
    delUdag = DeltaUDagger(
        alpha, woods_saxon, delta_U_matrix_element, spin_configurations,
        isospin_configurations
    )

    # Evaluate the \delta U^\dagger term for each q
    delta_U_array = np.zeros_like(q_array, dtype=complex)
    for i, q in enumerate(q_array):

        delta_U_array[i] = delUdag(q)
        
        percent = (i+1)/len(q_array)*100
        print(f"{percent:.2f} percent done with \delta U^\dagger.")

    return delta_U_array


def compute_overlap(
        nucleus_name, Z, N, alpha, kvnn, lamb, channels, kmax=15.0, kmid=3.0,
        ntot=120, generator='Wegner', neval=1e4, print_normalization=False,
        kvnn_hard=None, lambda_m=None, parametrization='match', ipm=False,
        save=False
):
    """Compute the spectroscopic overlap."""
    
    # Set table of Clebsch-Gordan coefficients
    jmax = 4  # This should cover nuclei as heavy as Ca48
    cg_table = compute_clebsch_gordan_table(jmax)
    
    # Set single-particle basis
    woods_saxon = WoodsSaxon(nucleus_name, Z, N, cg_table,
                             parametrization=parametrization)

    # Set points in q
    q_array, q_weights = momentum_mesh(10.0, 2.0, 120)
    
    # Compute the I term
    I_array = compute_I_term(alpha, q_array, woods_saxon)

    # IPM only
    if ipm:
        
        delta_U_array = np.zeros_like(q_array, dtype=complex)
    
    # Full overlap
    else:
        
        # Initialize \delta U matrix element class
        delta_U_matrix_element = DeltaUMatrixElement(
            cg_table, kvnn, kmax, kmid, ntot, generator, lamb, channels,
            kvnn_hard, lambda_m
        )
        
        # Compute the \delta U^\dagger term
        delta_U_array = compute_delta_U_dagger_term(
            alpha, q_array, woods_saxon, delta_U_matrix_element, neval
        )
    
    # Combine both terms for the total overlap [fm^3/2]
    overlap_array = I_array + delta_U_array

    # Option to print normalization of the total momentum distribution
    if print_normalization:
        normalization = compute_normalization(q_array, q_weights, overlap_array)
        print(f"Normalization = {normalization:.5f}.")
    
    # Option to save the momentum distribution as a .txt file
    if save:
        pass
    
    return q_array, q_weights, overlap_array


def compute_normalization(q_array, q_weights, overlap_array):
    """Compute the normalization of the overlap."""

    return 4 * np.pi * np.sum(
        q_weights * q_array ** 2 * np.abs(overlap_array) ** 2
    )


def save_overlap(
        nucleus_name, alpha, kvnn, lamb, q_array, q_weights, overlap_array,
        I_array, delta_U_array, delta_U_errors, kvnn_hard=None, lambda_m=None,
        parametrization='match'
):
    """Save the overlap along with the isolated contributions."""
    
    return None
    
    
def load_overlap(
        nucleus_name, alpha, kvnn, lamb, kvnn_hard=None, lambda_m=None,
        parametrization='match'
):
    """Load and return the overlap along with the isolated contributions."""

    return None


if __name__ == '__main__':
    
    # Nucleus
    nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'C12', 6, 6
    
    # Quantum state
    alpha = SingleParticleState(1, 0, 1/2, 1/2, 1/2)  # 1s_{1/2}
    # alpha = SingleParticleState(1, 1, 3/2, 1/2, 1/2)  # 1p_{3/2}

    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    
    # NN potential and momentum mesh
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120  # AV18
    
    # SRG \lambda value
    lamb = 1.5

    # Inverse-SRG evolution?
    kvnn_hard = None
    lambda_m = None
    
    # Woods-Saxon parametrization
    # prm = 'seminole'
    # prm = 'universal'
    prm = 'match'

    # Compute and save the overlap
    q_array, q_weights, overlap_array = compute_overlap(
        nucleus_name, Z, N, alpha, kvnn, lamb, channels, kvnn_hard=kvnn_hard,
        lambda_m=lambda_m, parametrization=prm, print_normalization=True,
        save=True
    )
    
    # # Testing IPM only
    # q_array, q_weights, overlap_array = compute_overlap(
    #     nucleus_name, Z, N, alpha, kvnn, lamb, channels,  kvnn_hard=kvnn_hard,
    #     lambda_m=lambda_m, parametrization=prm, print_normalization=True,
    #     ipm=True, save=False
    # )