#!/usr/bin/env python3

"""
File: test_momentum_distribution_script.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. In this version, we test several
`vegas` options to speed-up the code.

Last update: February 14, 2023

"""

# Python imports
import functools
import numpy as np
import numpy.linalg as la
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.special import spherical_jn, sph_harm
import shutil
from sympy.physics.quantum.cg import CG
import time
import vegas

# Imports from scripts
from scripts.integration import momentum_mesh, unattach_weights_from_matrix
from scripts.potentials import Potential
from scripts.srg import get_transformation
from scripts.tools import convert_l_to_string, coupled_channel
from scripts.woodsaxon import ws


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
    tau : float
        Isospin projection tau = 1/2 or -1/2.
    
    """
    
    
    def __init__(self, n, l, j, m_j, tau):
        
        # Check if m_j is valid
        if abs(m_j) > j:
            raise RuntimeError("m_j is not valid.")
            
        # Check that |\tau| = 1/2
        if abs(tau) != 1/2:
            raise RuntimeError("tau is not valid.")
            
        self.n = n
        self.l = l
        self.j = j
        self.m_j = m_j
        self.tau = tau
        
        if tau == 1/2:
            self.nucleon = 'proton'
        elif tau == -1/2:
            self.nucleon = 'neutron'
        
        
    def __eq__(self, sp_state):

        if (
            self.n == sp_state.n and self.l == sp_state.l
            and self.j == sp_state.j and self.m_j == sp_state.m_j
            and self.tau == sp_state.tau
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
    
    
class SingleParticleBasis:
    """
    Single-particle basis class. Handles the wave functions associated with the
    Wood-Saxon potential from the subroutine in woodsaxon.f90. Outputs wave
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
        Option to run the Wood-Saxon subroutine to generate orbital files.
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
        self.wood_saxon_directory = f"../data/wood_saxon/{nucleus_name}/"
        self.dr = rmax / ntab
        self.r_array = np.arange(self.dr, rmax + self.dr, self.dr)

        # Generate orbitals?
        if run_woodsaxon:
            
            self.run_wood_saxon_code(nucleus_name, Z, N, n_max, l_max, rmax,
                                     ntab)

            # Move output files to relevant directory
            shutil.move("ws_log", self.wood_saxon_directory + "ws_log")
            shutil.move("ws_pot", self.wood_saxon_directory + "ws_pot")
            shutil.move("ws_rho", self.wood_saxon_directory + "ws_rho")
                
        # Order single-particle states with lowest energy first
        self.order_sp_states(Z, N)
        
        # Organize wave functions in dictionary with the file name as the key
        self.sp_wfs = {}
        for sp_state in self.sp_states:
            # Wave functions are independent of m_j, so fix m_j=j
            if sp_state.m_j == sp_state.j:
                file_name = get_orbital_file_name(sp_state)
                if run_woodsaxon:
                    shutil.move(file_name,
                                self.wood_saxon_directory + file_name)
                data = np.loadtxt(self.wood_saxon_directory + file_name)
                # Use file name as the key
                self.sp_wfs[file_name] = data[:, 1]

            
    def run_wood_saxon_code(self, nucleus_name, Z, N, n_max, l_max, rmax, ntab):
        """Run Wood-Saxon code to generate data."""
        
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
    
        # Set parameters of the Wood-Saxon potential
        prm = np.zeros(shape=(2,9), order='F')
    
        # Starting with vws (p & n)
        if nucleus_name == 'He4':
            prm[:,0] = 76.8412
            # prm[:,0] = 90.0  # For bound 1p_{1/2} state
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
        ws_file = self.wood_saxon_directory + "ws_log"
    
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
                        )  # n, l, j, m_j, \tau
                    
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
                        )  # n, l, j, m_j, \tau
                    
                        self.sp_states.append(sp_state)
                    
                        if neutron_count < N:
                            self.occ_states.append(sp_state)
                            # Add up filled neutron states
                            neutron_count += 1
                        
                        
    def get_wf_rspace(self, sp_state, print_normalization=False):
        """Single-particle wave function in coordinate space."""
        
        # Orbital file name is the key
        u_array = self.sp_wfs[get_orbital_file_name(sp_state)]

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
            self, sp_state, print_normalization=False, interpolate=False,
            kmax=15.0, kmid=3.0, ntot=120
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


def compute_clebsch_gordan_table(j_max):
    """
    Calculate Clebsch-Gordan coefficients for all combinations of j and m_j up
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
        for m_1 in np.arange(-j_1, j_1+1, 1):
            for j_2 in j_array:
                for m_2 in np.arange(-j_2, j_2+1, 1):
                    for j_3 in j_array:
                        for m_3 in np.arange(-j_3, j_3+1, 1):
                            cg_table[(j_1,m_1,j_2,m_2,j_3,m_3)] = float(
                                CG(j_1,m_1,j_2,m_2,j_3,m_3).doit()
                            )
                                
    return cg_table


def get_orbital_file_name(sp_state):
    """Returns the file name of the orbital."""
        
    n, l, j = sp_state.n, sp_state.l, sp_state.j
    # Proton
    if sp_state.tau == 1/2:
        file_name = f"p.n{int(n-1)}.l{l}.j{int(2*j)}.orb"
    # Neutron
    elif sp_state.tau == -1/2:
        file_name = f"n.n{int(n-1)}.l{l}.j{int(2*j)}.orb"
        
    return file_name


def get_sp_wave_functions(sp_basis):
    """Set interpolating functions for s.p. wave functions \phi."""
    
    occ_states = sp_basis.occ_states

    phi_functions = {}
    for sp_state in occ_states: 
        file_name = get_orbital_file_name(sp_state)
        phi_functions[file_name] = sp_basis.get_wf_kspace(sp_state,
                                                          interpolate=True)
            
    return phi_functions


def psi(sp_state, q_vector, sigma, tau, cg_table, phi_functions):
    """Single-particle wave function including the Clebsch-Gordan coefficient 
    and spherical harmonic.
    """
        
    # Check that \tau = s.p. state \tau
    if sp_state.tau != tau:
        return 0
        
    # Unpack q_vector into magnitude and angles
    q = la.norm(q_vector)
    theta = np.arccos(q_vector[2]/q)
    phi = np.arctan2(q_vector[1], q_vector[0])
        
    # Calculate \phi_\alpha(q)
    phi_sp_wf = phi_functions[get_orbital_file_name(sp_state)](q)
    
    # m_l is determined by m_j and \sigma
    m_l = sp_state.m_j - sigma
        
    # Check that m_l is allowed
    if abs(m_l) > sp_state.l:
        return 0
        
    # Clebsch-Gordan coefficient
    cg = cg_table[(sp_state.l, m_l, 1/2, sigma, sp_state.j, sp_state.m_j)]
        
    # Spherical harmonic
    Y_lm = sph_harm(m_l, sp_state.l, phi, theta)

    return phi_sp_wf * cg * Y_lm


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
        
    # Set potential
    potential = Potential(kvnn, channel_arg, kmax, kmid, ntot)
    
    # Get momentum mesh
    k_array, k_weights = potential.load_mesh()
    
    # Initial and evolved Hamiltonians
    H_initial = potential.load_hamiltonian()
    if generator == 'Block-diag':
        H_evolved = potential.load_hamiltonian('srg', generator, 1.0,
                                               lambda_bd=lamb)
    else:
        H_evolved = potential.load_hamiltonian('srg', generator, lamb)
    
    # Get SRG transformation from Hamiltonians
    U_matrix_weights = get_transformation(H_initial, H_evolved)
    
    # Calculate \delta U = U - I
    I_matrix_weights = np.eye(len(U_matrix_weights), len(U_matrix_weights))
    if hermitian_conjugate:
        delU_matrix_weights = (U_matrix_weights - I_matrix_weights).T
    else:
        delU_matrix_weights = U_matrix_weights - I_matrix_weights

    # Get specific sub-block if coupled-channel
    n = ntot
    if channel in ['3S1-3D1', '3P2-3F2', '3D3-3G3']:
        delU_matrix = unattach_weights_from_matrix(
                k_array, k_weights, delU_matrix_weights[:n,n:2*n]
        )
    elif channel in ['3D1-3S1', '3F2-3P2', '3G3-3D3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[n:2*n,:n]
        )
    elif channel in ['3D1-3D1', '3F2-3F2', '3G3-3G3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[n:2*n,n:2*n]
        )
    else:
        delU_matrix = unattach_weights_from_matrix(k_array, k_weights,
                                                   delU_matrix_weights[:n,:n])
        
    # Interpolate \delta U(k, k')
    delU_func = RectBivariateSpline(k_array, k_array, delU_matrix)

    return delU_func


def get_delta_U_functions(
        channels, kvnn, kmax=15.0, kmid=3.0, ntot=120, generator='Wegner',
        lamb=1.35
):
    """Get \delta U and \delta U^\dagger functions."""

    delta_U_functions = {}
    delta_U_dagger_functions = {}
    for channel in channels:
        delta_U_functions[channel] = interpolate_delta_U(
            kvnn, channel, kmax, kmid, ntot, generator, lamb
        )
        delta_U_dagger_functions[channel] = (
            interpolate_delta_U(kvnn, channel, kmax, kmid, ntot, generator,
                                lamb, hermitian_conjugate=True)
        )
        
    return delta_U_functions, delta_U_dagger_functions


def get_channel_quantum_numbers(channel):
    """Gets the quantum numbers of a partial wave channel."""
    
    # Total orbital angular momentum L' = 0, 1, 2, ...
    if channel[1] == 'S':
        Lp = 0
    elif channel[1] == 'P':
        Lp = 1
    elif channel[1] == 'D':
        Lp = 2
    elif channel[1] == 'F':
        Lp = 3
    elif channel[1] == 'G':
        Lp = 4
    elif channel[1] == 'H':
        Lp = 5
    else:
        raise RuntimeError("Channel L' exceeds the range of the function.")
    
    # Total orbital angular momentum L = 0, 1, 2, ...
    if coupled_channel(channel[:3]):
        
        if channel[5] == 'S':
            L = 0
        elif channel[5] == 'P':
            L = 1
        elif channel[5] == 'D':
            L = 2
        elif channel[5] == 'F':
            L = 3
        elif channel[5] == 'G':
            L = 4
        elif channel[5] == 'H':
            L = 5
        else:
            raise RuntimeError("Channel L exceeds the range of the function.")
        
    # L = L' if the channel is not coupled
    else:
        
        L = Lp
        
    # Total spin S = 0 or 1
    S = int((int(channel[0])-1)/2)
    
    # Total angular momentum J = 0, 1, 2, ...
    J = int(channel[2])
            
    return Lp, L, S, J


def get_total_isospin(L, S):
    """Total isospin according to antisymmetry."""
    
    # Make sure [1-(-1)^(L+S+T)] is not zero
    if (1-(-1)**(L+S)) == 0:
        T = 1
    else:
        T = 0
        
    return T
    
    
def get_channel_str(Lp, L, S, J):
    """Gets the partial wave channel string given quantum numbers."""
    
    # Total orbital angular momentum L' = 0, 1, 2, ...
    if Lp == 0:
        Lp_str = 'S'
    elif Lp == 1:
        Lp_str = 'P'
    elif Lp == 2:
        Lp_str = 'D'
    elif Lp == 3:
        Lp_str = 'F'
    elif Lp == 4:
        Lp_str = 'G'
    elif Lp == 5:
        Lp_str = 'H'
    else:
        raise RuntimeError("Channel L' exceeds the range of the function.")
        
    channel = f"{2*S+1}{Lp_str}{J}"
    
    # Total orbital angular momentum L = 0, 1, 2, ...
    # L = L' if the channel is not coupled
    if channel in ['3S1', '3D1', '3P2', '3F2', '3D3', '3G3']:
        
        if L == 0:
            L_str = 'S'
        elif L == 1:
            L_str = 'P'
        elif L == 2:
            L_str = 'D'
        elif L == 3:
            L_str = 'F'
        elif L == 4:
            L_str = 'G'
        elif L == 5:
            L_str = 'H'
        else:
            raise RuntimeError("Channel L exceeds the range of the function.")
            
        channel += f"-{2*S+1}{L_str}{J}"

    return channel


# # Imposes L = L'
# def delta_U_matrix_element(
#         k_vector, kp_vector, sigma_1, sigma_2, sigma_3, sigma_4, tau_1, tau_2,
#         tau_3, tau_4, cg_table, channels, delta_U_functions,
#         delta_U_dagger_functions, hermitian_conjugate=False
# ):
#     """Returns the plane-wave matrix element of \delta U."""
        
#     # Check spin conservation
#     if sigma_1 + sigma_2 != sigma_3 + sigma_4:
#         return 0
        
#     # Check isospin conservation
#     if tau_1 + tau_2 != tau_3 + tau_4:
#         return 0
        
#     # Unpack k_vector and k'_vector into their magnitudes and angles
#     k = la.norm(k_vector)
#     theta_k = np.arccos(k_vector[2]/k)
#     phi_k = np.arctan2(k_vector[1], k_vector[0])
#     kp = la.norm(kp_vector)
#     theta_kp = np.arccos(kp_vector[2]/kp)
#     phi_kp = np.arctan2(kp_vector[1], kp_vector[0])

#     # Set total spin and total isospin projections
#     M_S = sigma_1 + sigma_2
#     M_T = tau_1 + tau_2
        
#     # Loop over partial wave channels
#     delta_U_matrix_element = 0
#     for channel in channels:
            
#         # Determine channel quantum numbers
#         L, Lp, S, J = get_channel_quantum_numbers(channel)
            
#         # Determine total isospin
#         T = get_total_isospin(L, S)
            
#         # Only proceed forward if |M_S| <= S, |M_T| <= T, and L=L'
#         if abs(M_S) <= S and abs(M_T) <= T and L == Lp:
            
#             # Calculate 1-(-1)^(L'+S+T) and 1-(-1)^(L+S+T) factors
#             factor = (1-(-1)**(Lp+S+T)) * (1-(-1)**(L+S+T))
            
#             # Spin CG's
#             spin_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
#             spin_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                    
#             # Isospin CG's
#             isospin_12_cg = cg_table[(1/2, tau_1, 1/2, tau_2, T, M_T)]
#             isospin_34_cg = cg_table[(1/2, tau_3, 1/2, tau_4, T, M_T)]

#             # \delta U(k, k') or \delta U^\dagger(k, k')
#             if hermitian_conjugate:
#                 delta_U_partial_wave = delta_U_dagger_functions[channel].ev(k,
#                                                                             kp)
#             else:
#                 delta_U_partial_wave = delta_U_functions[channel].ev(k, kp)

#             # Loop over possible M_J values
#             M_J_array = np.arange(-J, J+1, 1)
#             for M_J in M_J_array:
                
#                 # M_L and M_L' are fixed
#                 M_L = M_J - M_S
#                 M_Lp = M_J - M_S
                
#                 # Only proceed forward if |M_L| <= L and |M_L'| <= L'
#                 if abs(M_L) <= L and abs(M_Lp) <= Lp:
                    
#                     # L-S coupling CG's
#                     L_S_J_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
#                     Lp_S_J_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
                    
#                     # Spherical harmonics
#                     Y_L = sph_harm(M_L, L, phi_k, theta_k)
#                     Y_Lp = np.conj(sph_harm(M_Lp, Lp, phi_kp, theta_kp))

#                     delta_U_matrix_element += (
#                         1/2 * 2/np.pi * spin_12_cg * spin_34_cg * isospin_12_cg
#                         * isospin_34_cg * L_S_J_cg * Lp_S_J_cg * factor * Y_L
#                         * Y_Lp * delta_U_partial_wave
#                     )
                        
#     return delta_U_matrix_element


# Allows L != L'
def delta_U_matrix_element(
        k_vector, kp_vector, sigma_1, sigma_2, sigma_3, sigma_4, tau_1, tau_2,
        tau_3, tau_4, cg_table, channels, delta_U_functions,
        delta_U_dagger_functions, hermitian_conjugate=False
):
    """Returns the plane-wave matrix element of \delta U."""
        
    # Check spin conservation
    if sigma_1 + sigma_2 != sigma_3 + sigma_4:
        return 0
        
    # Check isospin conservation
    if tau_1 + tau_2 != tau_3 + tau_4:
        return 0
        
    # Unpack k_vector and k'_vector into their magnitudes and angles
    k = la.norm(k_vector)
    theta_k = np.arccos(k_vector[2]/k)
    phi_k = np.arctan2(k_vector[1], k_vector[0])
    kp = la.norm(kp_vector)
    theta_kp = np.arccos(kp_vector[2]/kp)
    phi_kp = np.arctan2(kp_vector[1], kp_vector[0])

    # Set total spin and total isospin projections
    M_S = sigma_1 + sigma_2
    M_T = tau_1 + tau_2
        
    # Loop over partial wave channels
    delta_U_matrix_element = 0
    for channel in channels:
            
        # Determine channel quantum numbers
        L, Lp, S, J = get_channel_quantum_numbers(channel)
            
        # Determine total isospin
        T = get_total_isospin(L, S)
            
        # Only proceed forward if |M_S| <= S and |M_T| <= T
        if abs(M_S) <= S and abs(M_T) <= T:
            
            # Calculate 1-(-1)^(L'+S+T) and 1-(-1)^(L+S+T) factors
            factor = (1-(-1)**(Lp+S+T)) * (1-(-1)**(L+S+T))
            
            # Spin CG's
            spin_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
            spin_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                    
            # Isospin CG's
            isospin_12_cg = cg_table[(1/2, tau_1, 1/2, tau_2, T, M_T)]
            isospin_34_cg = cg_table[(1/2, tau_3, 1/2, tau_4, T, M_T)]

            # \delta U(k, k') or \delta U^\dagger(k, k')
            if hermitian_conjugate:
                delta_U_partial_wave = delta_U_dagger_functions[channel].ev(k,
                                                                            kp)
            else:
                delta_U_partial_wave = delta_U_functions[channel].ev(k, kp)

            # Loop over possible M_J values
            M_J_array = np.arange(-J, J+1, 1)
            for M_J in M_J_array:
                
                # M_L and M_L' are fixed
                M_L = M_J - M_S
                M_Lp = M_J - M_S
                
                # Only proceed forward if |M_L| <= L and |M_L'| <= L'
                if abs(M_L) <= L and abs(M_Lp) <= Lp:
                    
                    # L-S coupling CG's
                    L_S_J_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
                    Lp_S_J_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
                    
                    # Spherical harmonics
                    Y_L = sph_harm(M_L, L, phi_k, theta_k)
                    Y_Lp = np.conj(sph_harm(M_Lp, Lp, phi_kp, theta_kp))

                    delta_U_matrix_element += (
                        1/2 * 2/np.pi * spin_12_cg * spin_34_cg * isospin_12_cg
                        * isospin_34_cg * L_S_J_cg * Lp_S_J_cg * factor * Y_L
                        * Y_Lp * delta_U_partial_wave
                    )
                        
    return delta_U_matrix_element


# # Imposes L = L''
# def delta_U2_matrix_elements(
#         k_vector, qK_vector, kp_vector, sigma_1, sigma_2, sigma, sigmap,
#         sigma_3, sigma_4, tau_1, tau_2, tau, taup, tau_3, tau_4, cg_table,
#         channels, delta_U_functions, delta_U_dagger_functions
# ):
#     """Returns the plane-wave matrix elements of \delta U^2."""
        
#     # Check spin conservation
#     if sigma_1 + sigma_2 != sigma + sigmap:
#         return 0
#     if sigma + sigmap != sigma_3 + sigma_4:
#         return 0
        
#     # Check isospin conservation
#     if tau_1 + tau_2 != tau + taup:
#         return 0
#     if tau + taup != tau_3 + tau_4:
#         return 0
        
#     # Unpack vectors into their magnitudes and angles
#     k = la.norm(k_vector)
#     theta_k = np.arccos(k_vector[2]/k)
#     phi_k = np.arctan2(k_vector[1], k_vector[0])
        
#     qK = la.norm(qK_vector)
#     theta_qK = np.arccos(qK_vector[2]/qK)
#     phi_qK = np.arctan2(qK_vector[1], qK_vector[0])
        
#     kp = la.norm(kp_vector)
#     theta_kp = np.arccos(kp_vector[2]/kp)
#     phi_kp = np.arctan2(kp_vector[1], kp_vector[0])

#     # Set total spin and total isospin projections
#     M_S = sigma + sigmap
#     M_T = tau + taup
        
#     # Loop over partial wave channels
#     delta_U2_matrix_element = 0
#     for channel in channels:
            
#         # Determine channel quantum numbers
#         L, Lp, S, J = get_channel_quantum_numbers(channel)
            
#         # Determine total isospin
#         T = get_total_isospin(L, S)
            
#         # Only proceed forward if |M_S| <= S and |M_T| <= T
#         if abs(M_S) <= S and abs(M_T) <= T:
            
#             # Calculate 1-(-1)^(L'+S+T) and 1-(-1)^(L+S+T) factors
#             factor = ((1-(-1)**(L+S+T)) * (1-(-1)**(Lp+S+T))
#                       * (1-(-1)**(Lp+S+T)) * (1-(-1)**(L+S+T)))
            
#             # Spin CG's
#             spin_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
#             spin_ssp_cg = cg_table[(1/2, sigma, 1/2, sigmap, S, M_S)]
#             spin_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                    
#             # Isospin CG's
#             isospin_12_cg = cg_table[(1/2, tau_1, 1/2, tau_2, T, M_T)]
#             isospin_ttp_cg = cg_table[(1/2, tau, 1/2, taup, T, M_T)]
#             isospin_34_cg = cg_table[(1/2, tau_3, 1/2, tau_4, T, M_T)]
                
#             # Get channel string for calling \delta U functions           
#             channel_1 = get_channel_str(L, Lp, S, J)
#             channel_2 = get_channel_str(Lp, L, S, J)
                
#             # \delta U(k, q-K/2) and \delta U^\dagger(q-K/2, k')
#             delta_U_partial_wave = delta_U_functions[channel_1].ev(k, qK)
#             delta_U_dag_partial_wave = (
#                 delta_U_dagger_functions[channel_2].ev(qK, kp)
#             )

#             # Loop over possible M_J values
#             M_J_array = np.arange(-J, J+1, 1)
#             for M_J in M_J_array:
                
#                 # M_L and M_L' are fixed
#                 M_L = M_J - M_S
#                 M_Lp = M_J - M_S
                
#                 # Only proceed forward if |M_L| <= L and |M_L'| <= L'
#                 if abs(M_L) <= L and abs(M_Lp) <= Lp:
                    
#                     # L-S coupling CG's
#                     L_S_J_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
#                     Lp_S_J_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
                    
#                     # Spherical harmonics
#                     Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
#                     Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
#                     Y_L_kp = np.conj(sph_harm(M_L, L, phi_kp, theta_kp))

#                     delta_U2_matrix_element += (
#                         (1/2 * 2/np.pi)**2 * spin_12_cg * spin_ssp_cg**2
#                         * spin_34_cg * isospin_12_cg * isospin_ttp_cg**2
#                         * isospin_34_cg * L_S_J_cg**2 * Lp_S_J_cg**2 * factor
#                         * Y_L_k * np.conj(Y_Lp_qK) * Y_Lp_qK * Y_L_kp
#                         * delta_U_partial_wave * delta_U_dag_partial_wave
#                     )
                        
#     return delta_U2_matrix_element


# # Allows L != L''
# def delta_U2_matrix_elements(
#         k_vector, qK_vector, kp_vector, sigma_1, sigma_2, sigma, sigmap,
#         sigma_3, sigma_4, tau_1, tau_2, tau, taup, tau_3, tau_4, cg_table,
#         channels, delta_U_functions, delta_U_dagger_functions
# ):
#     """Returns the plane-wave matrix elements of \delta U^2."""
        
#     # Check spin conservation
#     if sigma_1 + sigma_2 != sigma + sigmap:
#         return 0
#     if sigma + sigmap != sigma_3 + sigma_4:
#         return 0
        
#     # Check isospin conservation
#     if tau_1 + tau_2 != tau + taup:
#         return 0
#     if tau + taup != tau_3 + tau_4:
#         return 0
        
#     # Unpack vectors into their magnitudes and angles
#     k = la.norm(k_vector)
#     theta_k = np.arccos(k_vector[2]/k)
#     phi_k = np.arctan2(k_vector[1], k_vector[0])
        
#     qK = la.norm(qK_vector)
#     theta_qK = np.arccos(qK_vector[2]/qK)
#     phi_qK = np.arctan2(qK_vector[1], qK_vector[0])
        
#     kp = la.norm(kp_vector)
#     theta_kp = np.arccos(kp_vector[2]/kp)
#     phi_kp = np.arctan2(kp_vector[1], kp_vector[0])

#     # Set total spin and total isospin projections
#     M_S = sigma + sigmap
#     M_T = tau + taup
        
#     # Loop over partial wave channels
#     delta_U2_matrix_element = 0
#     # Partial wave channel for \delta U
#     for channel in channels:
            
#         # Determine channel quantum numbers
#         L, Lp, S, J = get_channel_quantum_numbers(channel)
            
#         # Determine total isospin
#         T = get_total_isospin(L, S)
        
#         # Only proceed forward if |M_S| <= S and |M_T| <= T
#         if abs(M_S) <= S and abs(M_T) <= T:
            
#             # Spin CG's
#             spin_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
#             spin_ssp_cg = cg_table[(1/2, sigma, 1/2, sigmap, S, M_S)]
#             spin_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                    
#             # Isospin CG's
#             isospin_12_cg = cg_table[(1/2, tau_1, 1/2, tau_2, T, M_T)]
#             isospin_ttp_cg = cg_table[(1/2, tau, 1/2, taup, T, M_T)]
#             isospin_34_cg = cg_table[(1/2, tau_3, 1/2, tau_4, T, M_T)]
            
#             # Get channel string for calling \delta U function          
#             channel_1 = get_channel_str(L, Lp, S, J)
            
#             # \delta U_{L,L'}(k, q-K/2)
#             delta_U_partial_wave = delta_U_functions[channel_1].ev(k, qK)
        
#             # Partial wave channel for \delta U^\dagger
#             for channel in channels:
            
#                 # Determine channel quantum numbers
#                 L_temp, Lpp, S_temp, J_temp = get_channel_quantum_numbers(
#                     channel
#                 )
            
#                 # Determine total isospin
#                 T_temp = get_total_isospin(L_temp, S_temp)
                
#                 # Make sure spin, isospin, and total angular momentum are
#                 # conserved, as well as the intermediate L
#                 if L_temp == Lp and S_temp == S and J_temp == J and T_temp == T:
                    
#                     # Calculate 1-(-1)^(L'+S+T) and 1-(-1)^(L+S+T) factors
#                     factor = ((1-(-1)**(L+S+T)) * (1-(-1)**(Lp+S+T))
#                               * (1-(-1)**(Lp+S+T)) * (1-(-1)**(Lpp+S+T)))
            
#                     # Get channel string for calling \delta U^\dagger function     
#                     channel_2 = get_channel_str(Lp, Lpp, S, J)

#                     # \delta U^\dagger(q-K/2, k')
#                     delta_U_dag_partial_wave = (
#                         delta_U_dagger_functions[channel_2].ev(qK, kp)
#                     )

#                     # Loop over possible M_J values
#                     M_J_array = np.arange(-J, J+1, 1)
#                     for M_J in M_J_array:
                
#                         # M_L, M_L', and M_L'' are fixed
#                         M_L = M_J - M_S
#                         M_Lp = M_J - M_S
#                         M_Lpp = M_J - M_S
                
#                         # Only proceed forward if |M_L| <= L, etc.
#                         if (abs(M_L) <= L and abs(M_Lp) <= Lp
#                             and abs(M_Lpp) <= Lpp):
                    
#                             # L-S coupling CG's
#                             L_S_J_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
#                             Lp_S_J_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
#                             Lpp_S_J_cg = cg_table[(Lpp, M_Lpp, S, M_S, J, M_J)]
                    
#                             # Spherical harmonics
#                             Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
#                             Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
#                             Y_Lpp_kp = np.conj(
#                                 sph_harm(M_Lpp, Lpp, phi_kp, theta_kp)
#                             )

#                             delta_U2_matrix_element += (
#                                 (1/2 * 2/np.pi)**2 * spin_12_cg * spin_ssp_cg**2
#                                 * spin_34_cg * isospin_12_cg * isospin_ttp_cg**2
#                                 * isospin_34_cg * L_S_J_cg * Lp_S_J_cg**2
#                                 * Lpp_S_J_cg * factor * Y_L_k * np.conj(Y_Lp_qK)
#                                 * Y_Lp_qK * Y_Lpp_kp * delta_U_partial_wave
#                                 * delta_U_dag_partial_wave
#                             )
                        
#     return delta_U2_matrix_element


# Allows L-L' L''-L''' structure
def delta_U2_matrix_elements(
        k_vector, qK_vector, kp_vector, sigma_1, sigma_2, sigma, sigmap,
        sigma_3, sigma_4, tau_1, tau_2, tau, taup, tau_3, tau_4, cg_table,
        channels, delta_U_functions, delta_U_dagger_functions
):
    """Returns the plane-wave matrix elements of \delta U^2."""
        
    # Check spin conservation
    if sigma_1 + sigma_2 != sigma + sigmap:
        return 0
    if sigma + sigmap != sigma_3 + sigma_4:
        return 0
        
    # Check isospin conservation
    if tau_1 + tau_2 != tau + taup:
        return 0
    if tau + taup != tau_3 + tau_4:
        return 0
        
    # Unpack vectors into their magnitudes and angles
    k = la.norm(k_vector)
    theta_k = np.arccos(k_vector[2]/k)
    phi_k = np.arctan2(k_vector[1], k_vector[0])
        
    qK = la.norm(qK_vector)
    theta_qK = np.arccos(qK_vector[2]/qK)
    phi_qK = np.arctan2(qK_vector[1], qK_vector[0])
        
    kp = la.norm(kp_vector)
    theta_kp = np.arccos(kp_vector[2]/kp)
    phi_kp = np.arctan2(kp_vector[1], kp_vector[0])

    # Set total spin and total isospin projections
    M_S = sigma + sigmap
    M_T = tau + taup
        
    # Loop over partial wave channels
    delta_U2_matrix_element = 0
    # Partial wave channel for \delta U
    for channel in channels:
            
        # Determine channel quantum numbers
        L, Lp, S, J = get_channel_quantum_numbers(channel)
            
        # Determine total isospin
        T = get_total_isospin(L, S)
        
        # Only proceed forward if |M_S| <= S and |M_T| <= T
        if abs(M_S) <= S and abs(M_T) <= T:
            
            # Spin CG's
            spin_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
            spin_ssp_cg = cg_table[(1/2, sigma, 1/2, sigmap, S, M_S)]
            spin_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                    
            # Isospin CG's
            isospin_12_cg = cg_table[(1/2, tau_1, 1/2, tau_2, T, M_T)]
            isospin_ttp_cg = cg_table[(1/2, tau, 1/2, taup, T, M_T)]
            isospin_34_cg = cg_table[(1/2, tau_3, 1/2, tau_4, T, M_T)]
            
            # Get channel string for calling \delta U function          
            channel_1 = get_channel_str(L, Lp, S, J)
            
            # \delta U_{L,L'}(k, q-K/2)
            delta_U_partial_wave = delta_U_functions[channel_1].ev(k, qK)
        
            # Partial wave channel for \delta U^\dagger
            for channel in channels:
            
                # Determine channel quantum numbers
                Lpp, Lppp, S_temp, J_temp = get_channel_quantum_numbers(
                    channel
                )
            
                # Determine total isospin
                T_temp = get_total_isospin(Lpp, S_temp)
                
                # Make sure spin, isospin, and total angular momentum are
                # conserved, as well as the intermediate L
                if S_temp == S and J_temp == J and T_temp == T:
                    
                    # Calculate 1-(-1)^(L'+S+T) and 1-(-1)^(L+S+T) factors
                    factor = ((1-(-1)**(L+S+T)) * (1-(-1)**(Lp+S+T))
                              * (1-(-1)**(Lpp+S+T)) * (1-(-1)**(Lppp+S+T)))
            
                    # Get channel string for calling \delta U^\dagger function     
                    channel_2 = get_channel_str(Lpp, Lppp, S, J)

                    # \delta U^\dagger(q-K/2, k')
                    delta_U_dag_partial_wave = (
                        delta_U_dagger_functions[channel_2].ev(qK, kp)
                    )

                    # Loop over possible M_J values
                    M_J_array = np.arange(-J, J+1, 1)
                    for M_J in M_J_array:
                
                        # M_L, M_L', and M_L'' are fixed
                        M_L = M_J - M_S
                        M_Lp = M_J - M_S
                        M_Lpp = M_J - M_S
                        M_Lppp = M_J - M_S
                
                        # Only proceed forward if |M_L| <= L, etc.
                        if (abs(M_L) <= L and abs(M_Lp) <= Lp
                            and abs(M_Lpp) <= Lpp and  abs(M_Lppp) <= Lppp):
                    
                            # L-S coupling CG's
                            L_S_J_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
                            Lp_S_J_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
                            Lpp_S_J_cg = cg_table[(Lpp, M_Lpp, S, M_S, J, M_J)]
                            Lppp_S_J_cg = cg_table[(Lppp, M_Lppp, S, M_S, J,
                                                    M_J)]
                    
                            # Spherical harmonics
                            Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
                            Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
                            Y_Lpp_qK = sph_harm(M_Lpp, Lpp, phi_qK, theta_qK)
                            Y_Lppp_kp = sph_harm(M_Lppp, Lppp, phi_kp, theta_kp)

                            delta_U2_matrix_element += (
                                (1/2 * 2/np.pi)**2 * spin_12_cg * spin_ssp_cg**2
                                * spin_34_cg * isospin_12_cg * isospin_ttp_cg**2
                                * isospin_34_cg * L_S_J_cg * Lp_S_J_cg
                                * Lpp_S_J_cg * Lppp_S_J_cg * factor * Y_L_k
                                * np.conj(Y_Lp_qK) * Y_Lpp_qK
                                * np.conj(Y_Lppp_kp) * delta_U_partial_wave
                                * delta_U_dag_partial_wave
                            )
                        
    return delta_U2_matrix_element


def delta_U_term_integrand(
        q, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions, momenta_array
):
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""
    
    # Spin projections
    spins = np.array([1/2, -1/2])
    
    # Choose z-axis to be along q_vector
    q_vector = np.array([0, 0, q])

    # Relative momenta k
    k, theta_k, phi_k = momenta_array[:3]
    k_vector = np.array([k * np.sin(theta_k) * np.cos(phi_k),
                          k * np.sin(theta_k) * np.sin(phi_k),
                          k * np.cos(theta_k)])
        
    # C.o.M. momenta K
    K, theta_K, phi_K = momenta_array[3:6]
    K_vector = np.array([K * np.sin(theta_K) * np.cos(phi_K),
                          K * np.sin(theta_K) * np.sin(phi_K),
                          K * np.cos(theta_K)])
        
    # Calculate the Jacobian determinant
    jacobian = k**2 * np.sin(theta_k) * K**2 * np.sin(theta_K)  
        
    # Sum over spin projections
    integrand = 0
    for sigma_1 in spins:
        for sigma_2 in spins:
            for sigma_3 in spins:
                sigma = sigma_1 + sigma_2 - sigma_3
                # Make sure \sigma = +1/2 or -1/2
                if sigma in spins:
                        
                    # Sum over occupied states \alpha and \beta
                    for alpha in occ_states:
                            
                        # \tau_\alpha, \tau, \tau, \tau_\alpha matrix element
                        matrix_element_ta_t_t_ta = delta_U_matrix_element(
                            k_vector, q_vector-K_vector/2, sigma_1, sigma_2,
                            sigma, sigma_3, alpha.tau, tau, tau, alpha.tau,
                            cg_table, channels, delta_U_functions,
                            delta_U_dagger_functions
                        )
                            
                        # \tau, \tau_\alpha, \tau_\alpha, \tau matrix element
                        matrix_element_dag_t_ta_ta_t = delta_U_matrix_element(
                            q_vector-K_vector/2, k_vector, sigma, sigma_3,
                            sigma_1, sigma_2, tau, alpha.tau, alpha.tau, tau,
                            cg_table, channels, delta_U_functions,
                            delta_U_dagger_functions, hermitian_conjugate=True
                        )
                            
                        for beta in occ_states:
                                
                            # \tau, \tau_\beta, \tau, \tau_\beta matrix element
                            matrix_element_t_tb_t_tb = delta_U_matrix_element(
                                k_vector, q_vector-K_vector/2, sigma_1, sigma_2,
                                sigma, sigma_3, tau, beta.tau, tau, beta.tau,
                                cg_table, channels, delta_U_functions,
                                delta_U_dagger_functions
                            )
                                
                            # \tau, \tau_\beta, \tau, \tau_\beta matrix element
                            matrix_element_dag_t_tb_t_tb = (
                                delta_U_matrix_element(
                                    q_vector-K_vector/2, k_vector, sigma,
                                    sigma_3, sigma_1, sigma_2, tau, beta.tau,
                                    tau, beta.tau, cg_table, channels,
                                    delta_U_functions, delta_U_dagger_functions,
                                    hermitian_conjugate=True
                                )
                            )
                                
                            # Full integrand of the \delta U and
                            # \delta U^\dagger linear terms
                            integrand += 1/2 * jacobian * (
                                matrix_element_t_tb_t_tb
                                * np.conj(
                                    psi(alpha, K_vector/2+k_vector,sigma_1, tau,
                                        cg_table, phi_functions)
                                )
                                * np.conj(
                                    psi(beta, K_vector/2-k_vector, sigma_2,
                                        beta.tau, cg_table, phi_functions)
                                )
                                * psi(beta, K_vector-q_vector, sigma_3,
                                      beta.tau, cg_table, phi_functions)
                                * psi(alpha, q_vector, sigma, tau, cg_table,
                                      phi_functions)
                                - matrix_element_ta_t_t_ta
                                * np.conj(
                                    psi(alpha, K_vector/2+k_vector, sigma_1,
                                        alpha.tau, cg_table, phi_functions)
                                )
                                * np.conj(
                                    psi(beta, K_vector/2-k_vector, sigma_2, tau,
                                        cg_table, phi_functions)
                                )
                                * psi(alpha, K_vector-q_vector, sigma_3,
                                      alpha.tau, cg_table, phi_functions)
                                * psi(beta, q_vector, sigma, tau, cg_table,
                                      phi_functions)
                                + matrix_element_dag_t_tb_t_tb
                                * np.conj(
                                    psi(alpha, q_vector, sigma, tau, cg_table,
                                        phi_functions)
                                )
                                * np.conj(
                                    psi(beta, K_vector-q_vector, sigma_3,
                                        beta.tau, cg_table, phi_functions)
                                )
                                * psi(beta, K_vector/2-k_vector, sigma_2,
                                      beta.tau, cg_table, phi_functions)
                                * psi(alpha, K_vector/2+k_vector, sigma_1, tau,
                                      cg_table, phi_functions)
                                - matrix_element_dag_t_ta_ta_t
                                * np.conj(
                                    psi(beta, q_vector, sigma, tau, cg_table,
                                        phi_functions)
                                )
                                * np.conj(
                                    psi(alpha, K_vector-q_vector, sigma_3,
                                        alpha.tau, cg_table, phi_functions)
                                )
                                * psi(beta, K_vector/2-k_vector, sigma_2, tau,
                                      cg_table, phi_functions)
                                * psi(alpha, K_vector/2+k_vector, sigma_1,
                                      alpha.tau, cg_table, phi_functions)
                            )

    return integrand.real


# Relies on calling delta U^2 matrix element function
def delta_U2_term_integrand(
        q, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions, momenta_array
):
    """Evaluate the integrand of the \delta U \delta U^\dagger term."""
    
    # Spin projections
    spins = np.array([1/2, -1/2])

    # Choose z-axis to be along q_vector
    q_vector = np.array([0, 0, q])

    # Relative momenta k
    k, theta_k, phi_k = momenta_array[:3]
    k_vector = np.array([k * np.sin(theta_k) * np.cos(phi_k),
                          k * np.sin(theta_k) * np.sin(phi_k),
                          k * np.cos(theta_k)])
        
    # Relative momenta k'
    kp, theta_kp, phi_kp = momenta_array[3:6]
    kp_vector = np.array([kp * np.sin(theta_kp) * np.cos(phi_kp),
                          kp * np.sin(theta_kp) * np.sin(phi_kp),
                          kp * np.cos(theta_kp)])
        
    # C.o.M. momenta K
    K, theta_K, phi_K = momenta_array[6:9]
    K_vector = np.array([K * np.sin(theta_K) * np.cos(phi_K),
                          K * np.sin(theta_K) * np.sin(phi_K),
                          K * np.cos(theta_K)])
        
    # Calculate the Jacobian determinant
    jacobian = (k**2 * np.sin(theta_k) * kp**2 * np.sin(theta_kp) * K**2
                * np.sin(theta_K))
        
    # Sum over spin projections
    integrand = 0 
    for sigma_1 in spins:
        for sigma_2 in spins:
            for sigma_3 in spins:
                sigma_4 = sigma_1 + sigma_2 - sigma_3
                for sigmap in spins:
                    sigma = sigma_1 + sigma_2 - sigmap
                        
                    # Make sure \sigma = +1/2 or -1/2
                    if sigma_4 in spins and sigma in spins:
                            
                        # Loop over \tau'
                        for taup in spins:
                        
                            # Sum over occupied states \alpha and \beta
                            for alpha in occ_states:
                                for beta in occ_states:
                                        
                                    # \tau_\alpha \tau_\beta \tau \tau'
                                    # \tau \tau' \tau_\alpha \tau_\beta
                                    # matrix elements
                                    matrix_elements_ta_tb = (
                                        delta_U2_matrix_elements(
                                            k_vector, q_vector-K_vector/2,
                                            kp_vector, sigma_1, sigma_2, sigma,
                                            sigmap, sigma_3, sigma_4, alpha.tau,
                                            beta.tau, tau, taup, alpha.tau,
                                            beta.tau, cg_table, channels,
                                            delta_U_functions,
                                            delta_U_dagger_functions
                                        )
                                    )
                                        
                                    # \tau_\alpha \tau_\beta \tau \tau'
                                    # \tau \tau' \tau_\beta \tau_\alpha
                                    # matrix elements
                                    matrix_elements_tb_ta = (
                                        delta_U2_matrix_elements(
                                            k_vector, q_vector-K_vector/2,
                                            kp_vector, sigma_1, sigma_2, sigma,
                                            sigmap, sigma_3, sigma_4, alpha.tau,
                                            beta.tau, tau, taup, beta.tau,
                                            alpha.tau, cg_table, channels,
                                            delta_U_functions,
                                            delta_U_dagger_functions
                                        )
                                    )

                                    # Full integrand of the \delta U
                                    # \delta U^\dagger term
                                    integrand += 1/4 * jacobian * (
                                        np.conj(
                                            psi(
                                                alpha, K_vector/2+k_vector,
                                                sigma_1, alpha.tau, cg_table,
                                                phi_functions
                                            )
                                        )
                                        * np.conj(
                                            psi(
                                                beta, K_vector/2-k_vector,
                                                sigma_2, beta.tau, cg_table,
                                                phi_functions
                                            )
                                        )
                                        * (
                                            matrix_elements_ta_tb
                                            * psi(
                                                beta, K_vector/2-kp_vector,
                                                sigma_4, beta.tau, cg_table,
                                                phi_functions
                                            )
                                            * psi(
                                                alpha, K_vector/2+kp_vector,
                                                sigma_3, alpha.tau, cg_table,
                                                phi_functions
                                            )
                                            - matrix_elements_tb_ta
                                            * psi(
                                                alpha, K_vector/2-kp_vector,
                                                sigma_4, alpha.tau, cg_table,
                                                phi_functions
                                            )
                                            * psi(
                                                beta, K_vector/2+kp_vector,
                                                sigma_3, beta.tau, cg_table,
                                                phi_functions
                                            )
                                        )
                                    )

    return integrand.real


def compute_I_term(q_array, tau, occ_states, cg_table, phi_functions):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in occ_states:
                
            # m_l is determined by m_j and \sigma
            m_l_alpha = alpha.m_j - sigma
                
            # Check that |m_l| <= l
            if abs(m_l_alpha) <= alpha.l:
                    
                # Loop over q
                psi_alpha_array = np.zeros_like(q_array, dtype='complex')
                for i, q in enumerate(q_array):
                        
                    # Choose z-axis to be along q_vector
                    q_vector = np.array([0, 0, q])
                    
                    # Single-particle wave function
                    psi_alpha_array[i] = psi(alpha, q_vector, sigma, tau,
                                             cg_table, phi_functions)
                    
                I_array += abs(psi_alpha_array)**2
                    
    return I_array


def compute_delta_U_terms(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions
):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
        
    # Set-up integration with vegas
    k_limits = [0, 10]  # Relative momenta up to 10 fm^-1
    K_limits = [0, 3]  # C.o.M. momenta up to 3 fm^-1
    theta_limits = [0, np.pi]
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = functools.partial(
            delta_U_term_integrand, q, tau, occ_states, cg_table, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
            
        # Train the integrator
        integ(integrand, nitn=5, neval=200)
        # Final result
        result = integ(integrand, nitn=10, neval=1e3)
        
        # # Train the integrator
        # integ(integrand, nitn=5, neval=1e3)
        # # Final result
        # result = integ(integrand, nitn=10, neval=3e3)

        delta_U_array[i] = result.mean
        delta_U_errors[i] = result.sdev
            
        t1 = time.time()

        percent = (i+1)/len(q_array)*100
        mins = (t1-t0)/60
        print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                    
    return delta_U_array, delta_U_errors


def compute_delta_U2_term(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions
):
    """Compute the \delta U * n(q) * \delta U^\dagger term."""
        
    delta_U2_array = np.zeros_like(q_array)
    delta_U2_errors = np.zeros_like(q_array)
        
    # Set-up integration with vegas
    k_limits = [0, 10]  # Relative momenta up to 10 fm^-1
    K_limits = [0, 3]  # C.o.M. momenta up to 3 fm^-1
    theta_limits = [0, np.pi]
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = functools.partial(
            delta_U2_term_integrand, q, tau, occ_states, cg_table, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
            
        # # Train the integrator
        # integ(integrand, nitn=5, neval=200)
        # # Final result
        # result = integ(integrand, nitn=10, neval=1e3)
        
        # Train the integrator
        integ(integrand, nitn=5, neval=1e3)
        # Final result
        result = integ(integrand, nitn=10, neval=3e3)

        delta_U2_array[i] = result.mean
        delta_U2_errors[i] = result.sdev
            
        t1 = time.time()

        percent = (i+1)/len(q_array)*100
        mins = (t1-t0)/60
        print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                    
    return delta_U2_array, delta_U2_errors


def compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, print_normalization=False,
        ipm_only=False, save=True
):
    """Compute the single-nucleon momentum distribution."""
    
    # Compute table of Clebsch-Gordan coefficients
    cg_table = compute_clebsch_gordan_table(4)
    print("Done calculating Clebsch-Gordan table.\n")
    
    # Set-up single-particle states
    sp_basis = SingleParticleBasis(nucleus_name, Z, N, run_woodsaxon=False)
    phi_functions = get_sp_wave_functions(sp_basis)
    
    # Set-up \delta U and \delta U^\dagger functions
    delta_U_functions, delta_U_dagger_functions = get_delta_U_functions(
        channels, kvnn
    )
    
    # Set momentum mesh
    q_array, q_weights = momentum_mesh(8.0, 2.0, 40)
    # q_array, q_weights = momentum_mesh(8.0, 2.0, 60)
    
    # Compute the I term
    I_array = compute_I_term(q_array, tau, sp_basis.occ_states, cg_table,
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
        delta_U_array, delta_U_errors = compute_delta_U_terms(
            q_array, tau, sp_basis.occ_states, cg_table, channels,
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
            q_array, tau, sp_basis.occ_states, cg_table, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
        t4 = time.time()
        print(f"Done after {(t4-t3)/3600:.3f} hours.\n")
        
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
            nucleus_name, tau, q_array, q_weights, n_array, n_errors, I_array,
            delta_U_array, delta_U_errors, delta_U2_array, delta_U2_errors
        )
    
    return q_array, q_weights, n_array, n_errors


def compute_normalization(q_array, q_weights, n_array):
    """Compute the normalization of the momentum distribution."""

    return 4*np.pi * np.sum(q_weights * q_array**2 * n_array)


def save_momentum_distribution(
        nucleus_name, tau, q_array, q_weights, n_array, n_errors, I_array,
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

    np.savetxt(f"{nucleus_name}_{nucleon}_momentum_distribution.txt",
               data, header=hdr)


def load_momentum_distribution(nucleus_name, nucleon):
    """Load and return the momentum distribution along with the isolated
    contributions.
    """
    
    data = np.loadtxt(f"{nucleus_name}_{nucleon}_momentum_distribution.txt")
    
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
    
    # Set nucleus, nucleon, partial wave channels, and potential
    nucleus_name, Z, N = 'He4', 2, 2
    tau = 1/2
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    kvnn = 6  # AV18

    # He4
    q_array, q_weights, n_array, n_errors = compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, print_normalization=True,
        save=True
    )
    
    # # O16
    # nucleus_name, Z, N = 'O16', 8, 8
    # q_array, q_weights, n_array, n_errors = compute_momentum_distribution(
    #     nucleus_name, Z, N, tau, channels, kvnn, print_normalization=True,
    #     save=True
    # )
    
    