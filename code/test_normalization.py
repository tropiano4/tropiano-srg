#!/usr/bin/env python3

"""
File: test_momentum_distribution_script.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. In this version, we test what affects
the normalization of the \delta U linear and quadratic terms, respectively.
These terms should cancel out after integrating the momentum distribution giving
a normalization of Z or N.

Last update: March 14, 2023

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
from scripts.tools import convert_l_to_string, coupled_channel, replace_periods
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


def load_H_no_coupling(kvnn, channel, generator, lamb, kmax, kmid, ntot):
    """Loads the evolved Hamiltonian assuming there are no coupled channels."""
    
    file_name = (
        f"H_evolved_kvnn_{kvnn}_{channel}_no_coupling_{generator}_lamb_{lamb}"
        f"_kmax_{kmax}_kmid_{kmid}_ntot_{ntot}_ode_BDF_wrt_lambda"
    )
    
    H_evolved = np.loadtxt("./test_srg/" + replace_periods(file_name) + ".txt")
    
    return H_evolved


def interpolate_delta_U(
        kvnn, channel, kmax, kmid, ntot, generator, lamb,
        hermitian_conjugate=False
):
    """Interpolate \delta U(k, k') for the given channel."""
    
    # Set channel argument to be compatible with potential functions
    if channel[:3] in ['3S1', '3D1']:  # 3S1-3D1
        channel_arg = '3S1'
    else:  # 1S0
        channel_arg = channel[:3]
        
    # Set potential
    potential = Potential(kvnn, channel_arg, kmax, kmid, ntot)
    
    # Get momentum mesh
    k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
    
    # Initial and evolved Hamiltonians
    if channel == '3D1-3D1':
        H_initial = potential.load_hamiltonian()[ntot:, ntot:]
        H_evolved = load_H_no_coupling(kvnn, channel_arg, generator, lamb, kmax,
                                       kmid, ntot)[ntot:, ntot:]
    else:
        H_initial = potential.load_hamiltonian()[:ntot, :ntot]
        H_evolved = load_H_no_coupling(kvnn, channel_arg, generator, lamb, kmax,
                                       kmid, ntot)[:ntot, :ntot]
        
    # Get SRG transformation from Hamiltonians
    U_matrix_weights = get_transformation(H_initial, H_evolved)
    
    # Calculate \delta U = U - I
    I_matrix_weights = np.eye(len(U_matrix_weights), len(U_matrix_weights))
    if hermitian_conjugate:
        delU_matrix_weights = (U_matrix_weights - I_matrix_weights).T
    else:
        delU_matrix_weights = U_matrix_weights - I_matrix_weights

    # Unattach integration factors
    delU_matrix = unattach_weights_from_matrix(k_array, k_weights,
                                               delU_matrix_weights)
        
    # Interpolate \delta U(k, k') using a cubic spline
    delU_func = RectBivariateSpline(k_array, k_array, delU_matrix)
                                        
    return delU_func


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


def compute_I_term(q_array, occ_states, cg_table, phi_functions):
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
                    psi_alpha_array[i] = psi(alpha, q_vector, sigma, alpha.tau,
                                             cg_table, phi_functions)
                    
                I_array += abs(psi_alpha_array)**2
                    
    return I_array


def delta_U_quantum_numbers(cg_table, occ_states, channel):
    """Returns a list of every combination of quantum numbers in the \delta U
    linear terms that gives a non-zero product of Clebsch-Gordan coefficients.
    """
    
    spins = np.array([1/2, -1/2])
    quantum_number_combinations = []
    
    L, Lp, S, J = get_channel_quantum_numbers(channel)
    T = get_total_isospin(L, S)
    
    # L S T factor
    lst_factor = 1-(-1)**(L+S+T)
    # L' S T factor
    lpst_factor = 1-(-1)**(Lp+S+T)

    # Many loops over quantum numbers
    for sigma_1 in spins:
        for sigma_2 in spins:
            
            M_S = sigma_1 + sigma_2
            if abs(M_S) <= S:
                
                # \sigma_1 \sigma_2 S M_S CG
                sig_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
            
                for sigma_3 in spins:
                    for sigma_4 in spins:
                        
                        if (sigma_3 + sigma_4) == M_S:
                        
                            # \sigma_3 \sigma_4 S M_S CG
                            sig_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                    
                            for alpha in occ_states:
                    
                                # \sigma_1 l_\alpha CG
                                if abs(alpha.m_j-sigma_1) > alpha.l:
                                    sig1_alpha_cg = 0
                                else:
                                    sig1_alpha_cg = cg_table[(
                                        alpha.l, alpha.m_j-sigma_1, 1/2, sigma_1, alpha.j, alpha.m_j
                                    )]
                                
                                # \sigma_3 l_\alpha CG
                                if abs(alpha.m_j-sigma_3) > alpha.l:
                                    sig3_alpha_cg = 0
                                else:
                                    sig3_alpha_cg = cg_table[(
                                        alpha.l, alpha.m_j-sigma_3, 1/2, sigma_3, alpha.j, alpha.m_j
                                    )]
                                    
                                for beta in occ_states:
                    
                                    # \sigma_2 l_\beta CG
                                    if abs(beta.m_j-sigma_2) > beta.l:
                                        sig2_beta_cg = 0
                                    else:
                                        sig2_beta_cg = cg_table[(
                                            beta.l, beta.m_j-sigma_2, 1/2, sigma_2, beta.j, beta.m_j
                                        )]
                                
                                    # \sigma_4 l_\beta CG
                                    if abs(beta.m_j-sigma_4) > beta.l:
                                        sig4_beta_cg = 0
                                    else:
                                        sig4_beta_cg = cg_table[(
                                            beta.l, beta.m_j-sigma_4, 1/2, sigma_4, beta.j, beta.m_j
                                        )]
                                        
                                    M_T = alpha.tau + beta.tau
                                    if abs(M_T) <= T:

                                        # \tau_\alpha \tau_\beta T M_T CG
                                        tau_ab_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, T, M_T)]
                                    
                                        for M_J in np.arange(-J, J+1, 1):
                                                
                                            M_L = M_J - M_S
                                            M_Lp = M_L
                                            
                                            # L S J CG
                                            if abs(M_L) > L:
                                                lsj_cg = 0
                                            else:
                                                lsj_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
                                            # L' S J CG
                                            if abs(M_Lp) > Lp:
                                                lpsj_cg = 0
                                            else:
                                                lpsj_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
 
                                            # product
                                            product = (
                                                lst_factor * lpst_factor
                                                * sig_12_cg * sig_34_cg
                                                * tau_ab_cg ** 2
                                                * lsj_cg * lpsj_cg
                                                * sig1_alpha_cg * sig2_beta_cg
                                                * sig3_alpha_cg * sig4_beta_cg
                                            )

                                            if product != 0:
                                                            
                                                d = {
                                                    'sigma_1': sigma_1,
                                                    'sigma_2': sigma_2,
                                                    'sigma_3': sigma_3,
                                                    'sigma_4': sigma_4,
                                                    'alpha': alpha,
                                                    'beta': beta,
                                                    'T': T, 'M_T': M_T,
                                                    'S': S, 'M_S': M_S,
                                                    'L': L, 'M_L': M_L,
                                                    'Lp': Lp, 'M_Lp': M_Lp,
                                                    'J': J, 'M_J': M_J
                                                }
                                                            
                                                quantum_number_combinations.append(d)
    
    return quantum_number_combinations


def delta_U_term_integrand(
        q, channel, quantum_numbers, cg_table, phi_functions, delta_U_function,
        delta_U_dagger_function, momenta_array
):
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""

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
    
    # Calculate vector q-K/2
    qK_vector = q_vector - K_vector/2
    qK = la.norm(qK_vector)
    theta_qK = np.arccos(qK_vector[2]/qK)
    phi_qK = np.arctan2(qK_vector[1], qK_vector[0])
        
    # Calculate the Jacobian determinant
    jacobian = k**2 * np.sin(theta_k) * K**2 * np.sin(theta_K)
    
    integrand = 0
    for d in quantum_numbers:
        
        # Unpack dictionary
        sigma_1 = d['sigma_1']
        sigma_2 = d['sigma_2']
        sigma_3 = d['sigma_3']
        sigma_4 = d['sigma_4']
        
        alpha = d['alpha']
        beta = d['beta']
        
        T, M_T = d['T'], d['M_T']
        S, M_S = d['S'], d['M_S']
        L, M_L = d['L'], d['M_L']
        Lp, M_Lp = d['Lp'], d['M_Lp']
        J, M_J = d['J'], d['M_J']
        
        # \sigma_1 \sigma_2 S M_S CG
        sig_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
        # \sigma_3 \sigma_4 S M_S CG
        sig_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
            
        # \tau_\alpha \tau_\beta T M_T CG
        tau_ab_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, T, M_T)]
        
        # L S J CG
        lsj_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
        # L' S J CG
        lpsj_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
        
        # L S T factor
        lst_factor = 1-(-1)**(L+S+T)
        # L' S T factor
        lpst_factor = 1-(-1)**(Lp+S+T)
        
        # Spherical harmonics
        Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
        Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
        # Y_L_k = 1/np.sqrt(4*np.pi)
        # Y_Lp_qK = 1/np.sqrt(4*np.pi)
        
        # \delta U_{L,L'}(k, q-K/2)
        delta_U_partial_wave = delta_U_function.ev(k, qK)
        # \delta U^\dagger_{L,L'}(q-K/2, k)
        delta_U_dag_partial_wave = delta_U_dagger_function.ev(qK, k)
        
        # Single-particle wave functions
        psi_alpha_k1 = psi(alpha, K_vector/2+k_vector, sigma_1, alpha.tau,
                           cg_table, phi_functions)
        psi_beta_k2 = psi(beta, K_vector/2-k_vector, sigma_2, beta.tau,
                          cg_table, phi_functions)
        psi_beta_Kq4 = psi(beta, K_vector-q_vector, sigma_4, beta.tau, cg_table,
                           phi_functions)
        psi_alpha_q3 = psi(alpha, q_vector, sigma_3, alpha.tau, cg_table,
                          phi_functions)
        
        psi_alpha_q1 = psi(alpha, q_vector, sigma_1, alpha.tau, cg_table,
                           phi_functions)
        psi_beta_Kq2 = psi(beta, K_vector-q_vector, sigma_2, beta.tau, cg_table,
                           phi_functions)
        psi_beta_k4 = psi(beta, K_vector/2-k_vector, sigma_4, beta.tau,
                          cg_table, phi_functions)
        psi_alpha_k3 = psi(alpha, K_vector/2+k_vector, sigma_3, alpha.tau,
                           cg_table, phi_functions)
    
        integrand += (
            (1/2*2/np.pi) * jacobian * sig_12_cg * sig_34_cg * tau_ab_cg ** 2
            * lsj_cg * lpsj_cg * lst_factor * lpst_factor * (
                Y_L_k * np.conj(Y_Lp_qK) * delta_U_partial_wave * psi_alpha_k1
                * psi_beta_k2 * psi_beta_Kq4 * psi_alpha_q3
                + Y_Lp_qK * np.conj(Y_L_k) * delta_U_dag_partial_wave
                * psi_alpha_q1 * psi_beta_Kq2 * psi_beta_k4 * psi_alpha_k3
            )
        )
    
    return integrand.real


def compute_delta_U_term(
        q_array, cg_table, occ_states, channel, phi_functions, delta_U_function,
        delta_U_dagger_function
):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    quantum_numbers = delta_U_quantum_numbers(cg_table, occ_states, channel)
        
    # Set-up integration with vegas
    k_array, _ = momentum_mesh(15.0, 3.0, 120)
    k_max, k_min = max(k_array), min(k_array)
    k_limits = [k_min, k_max]  # Relative momenta up to 15^-1
    K_limits = [0, 3]  # C.o.M. momenta up to 3 fm^-1
    theta_limits = [0, np.pi]
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):

        integrand = functools.partial(
            delta_U_term_integrand, q, channel, quantum_numbers, cg_table,
            phi_functions, delta_U_function, delta_U_dagger_function
        )
            
        # # Train the integrator
        # integ(integrand, nitn=5, neval=200)
        # # Final result
        # result = integ(integrand, nitn=10, neval=1e3)
        
        # Train the integrator
        integ(integrand, nitn=5, neval=1e3)
        # Final result
        result = integ(integrand, nitn=10, neval=3e3)

        delta_U_array[i] = result.mean
        delta_U_errors[i] = result.sdev
           
    return delta_U_array, delta_U_errors


def delta_U2_quantum_numbers(cg_table, occ_states, channel):
    """Returns a list of every combination of quantum numbers in the \delta U
    linear terms that gives a non-zero product of Clebsch-Gordan coefficients.
    """
    
    # UPDATE THIS
    spins = np.array([1/2, -1/2])
    quantum_number_combinations = []
    
    L, Lp, S, J = get_channel_quantum_numbers(channel)
    T = get_total_isospin(L, S)
    
    # L S T factor
    lst_factor = 1-(-1)**(L+S+T)
    # L' S T factor
    lpst_factor = 1-(-1)**(Lp+S+T)

    # Many loops over quantum numbers
    for sigma_1 in spins:
        for sigma_2 in spins:
            
            M_S = sigma_1 + sigma_2
            if abs(M_S) <= S:
                
                # \sigma_1 \sigma_2 S M_S CG
                sig_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
            
                for sigma_3 in spins:
                    for sigma_4 in spins:
                        
                        if (sigma_3 + sigma_4) == M_S:
                        
                            # \sigma_3 \sigma_4 S M_S CG
                            sig_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
                            
                            for sigma_5 in spins:
                                for sigma_6 in spins:
                                    
                                    if (sigma_5 + sigma_6) == M_S:
                                        
                                        # \sigma_5 \sigma_6 S M_S CG
                                        sig_56_cg = cg_table[(1/2, sigma_5, 1/2, sigma_6, S, M_S)]
                                        
                                        for alpha in occ_states:
                                            
                                            # \sigma_1 l_\alpha CG
                                            if abs(alpha.m_j-sigma_1) > alpha.l:
                                                sig1_alpha_cg = 0
                                            else:
                                                sig1_alpha_cg = cg_table[(
                                                    alpha.l, alpha.m_j-sigma_1, 1/2, sigma_1, alpha.j, alpha.m_j
                                                )]
                                            
                                            # \sigma_3 l_\alpha CG
                                            if abs(alpha.m_j-sigma_3) > alpha.l:
                                                sig3_alpha_cg = 0
                                            else:
                                                sig3_alpha_cg = cg_table[(
                                                    alpha.l, alpha.m_j-sigma_3, 1/2, sigma_3, alpha.j, alpha.m_j
                                                )]
                                                
                                            for beta in occ_states:
                                                
                                                # \sigma_2 l_\beta CG
                                                if abs(beta.m_j-sigma_2) > beta.l:
                                                    sig2_beta_cg = 0
                                                else:
                                                    sig2_beta_cg = cg_table[(
                                                        beta.l, beta.m_j-sigma_2, 1/2, sigma_2, beta.j, beta.m_j
                                                    )]
                                            
                                                # \sigma_4 l_\beta CG
                                                if abs(beta.m_j-sigma_4) > beta.l:
                                                    sig4_beta_cg = 0
                                                else:
                                                    sig4_beta_cg = cg_table[(
                                                        beta.l, beta.m_j-sigma_4, 1/2, sigma_4, beta.j, beta.m_j
                                                    )]

                                                M_T = alpha.tau + beta.tau
                                                if abs(M_T) <= T:

                                                    # \tau_\alpha \tau_\beta T M_T CG
                                                    tau_ab_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, T, M_T)]
                                                    
                                                    for tau in spins:
                                                        for taup in spins:
                                                            
                                                            if (tau+taup) == M_T:
                                                                
                                                                # \tau \tau' T M_T CG
                                                                tau_ttp_cg = cg_table[(1/2, tau, 1/2, taup, T, M_T)]
                                    
                                                                for M_J in np.arange(-J, J+1, 1):
                                                
                                                                    M_L = M_J - M_S
                                                                    M_Lp = M_L
                                            
                                                                    # L S J CG
                                                                    if abs(M_L) > L:
                                                                        lsj_cg = 0
                                                                    else:
                                                                        lsj_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
                                                                    # L' S J CG
                                                                    if abs(M_Lp) > Lp:
                                                                        lpsj_cg = 0
                                                                    else:
                                                                        lpsj_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
                                                                        
                                                                    for M_Jp in np.arange(-J, J+1, 1):
                                                                        
                                                                        M_Lpp = M_Jp - M_S
                                                                        M_Lppp = M_Lpp
                                                                        
                                                                        # L'' S J CG
                                                                        if abs(M_Lpp) > Lp:
                                                                            lppsj_cg = 0
                                                                        else:
                                                                            lppsj_cg = cg_table[(Lp, M_Lpp, S, M_S, J, M_J)]
                                                                        # L''' S J CG
                                                                        if abs(M_Lppp) > L:
                                                                            lpppsj_cg = 0
                                                                        else:
                                                                            lpppsj_cg = cg_table[(L, M_Lppp, S, M_S, J, M_J)]
 
                                                                        # product
                                                                        product = (
                                                                            lst_factor ** 2
                                                                            * lpst_factor ** 2
                                                                            * sig_12_cg 
                                                                            * sig_56_cg **2
                                                                            * sig_34_cg
                                                                            * tau_ab_cg ** 2
                                                                            * tau_ttp_cg ** 2
                                                                            * lsj_cg
                                                                            * lpsj_cg
                                                                            * lppsj_cg
                                                                            * lpppsj_cg
                                                                            * sig1_alpha_cg
                                                                            * sig2_beta_cg
                                                                            * sig3_alpha_cg
                                                                            * sig4_beta_cg
                                                                        )

                                                                        if product != 0:
                                                            
                                                                            d = {
                                                                                'sigma_1': sigma_1,
                                                                                'sigma_2': sigma_2,
                                                                                'sigma_3': sigma_3,
                                                                                'sigma_4': sigma_4,
                                                                                'sigma_5': sigma_5,
                                                                                'sigma_6': sigma_6,
                                                                                'alpha': alpha,
                                                                                'beta': beta,
                                                                                'tau': tau,
                                                                                'taup': taup,
                                                                                'T': T, 'M_T': M_T,
                                                                                'S': S, 'M_S': M_S,
                                                                                'L': L, 'M_L': M_L,
                                                                                'Lp': Lp, 'M_Lp': M_Lp,
                                                                                'M_Lpp': M_Lpp, 'M_Lppp': M_Lppp,
                                                                                'J': J, 'M_J': M_J, 'M_Jp': M_Jp
                                                                                }
                                                            
                                                                            quantum_number_combinations.append(d)
    
    return quantum_number_combinations


def delta_U2_term_integrand(
        q, channel, quantum_numbers, cg_table, phi_functions, delta_U_function,
        delta_U_dagger_function, momenta_array
):
    """Evaluate the integrand of the \delta U \delta U^\dagger term."""

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
    
    # Calculate vector q-K/2
    qK_vector = q_vector - K_vector/2
    qK = la.norm(qK_vector)
    theta_qK = np.arccos(qK_vector[2]/qK)
    phi_qK = np.arctan2(qK_vector[1], qK_vector[0])
        
    # Calculate the Jacobian determinant
    jacobian = (k**2 * np.sin(theta_k) * kp**2 * np.sin(theta_kp) * K**2
                * np.sin(theta_K))
    
    integrand = 0
    for d in quantum_numbers:
        
        # Unpack dictionary
        sigma_1 = d['sigma_1']
        sigma_2 = d['sigma_2']
        sigma_3 = d['sigma_3']
        sigma_4 = d['sigma_4']
        sigma_5 = d['sigma_5']
        sigma_6 = d['sigma_6']
        
        alpha = d['alpha']
        beta = d['beta']
        tau = d['tau']
        taup = d['taup']

        T, M_T = d['T'], d['M_T']
        S, M_S = d['S'], d['M_S']
        L, M_L, M_Lppp = d['L'], d['M_L'], d['M_Lppp']
        Lp, M_Lp, M_Lpp = d['Lp'], d['M_Lp'], d['M_Lpp']
        J, M_J, M_Jp = d['J'], d['M_J'], d['M_Jp']
        
        # \sigma_1 \sigma_2 S M_S CG
        sig_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
        # \sigma_3 \sigma_4 S M_S CG
        sig_34_cg = cg_table[(1/2, sigma_3, 1/2, sigma_4, S, M_S)]
        # \sigma_5 \sigma_6 S M_S CG
        sig_56_cg = cg_table[(1/2, sigma_5, 1/2, sigma_6, S, M_S)]
            
        # \tau_\alpha \tau_\beta T M_T CG
        tau_ab_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, T, M_T)]
        # \tau \tau' T M_T CG
        tau_ttp_cg = cg_table[(1/2, tau, 1/2, taup, T, M_T)]

        # L S J CG
        lsj_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
        lpppsj_cg = cg_table[(L, M_Lppp, S, M_S, J, M_Jp)]
        # L' S J CG
        lpsj_cg = cg_table[(Lp, M_Lp, S, M_S, J, M_J)]
        lppsj_cg = cg_table[(Lp, M_Lpp, S, M_S, J, M_Jp)]
        
        # L S T factor
        lst_factor = 1-(-1)**(L+S+T)
        # L' S T factor
        lpst_factor = 1-(-1)**(Lp+S+T)
        
        # Spherical harmonics
        Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
        Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
        Y_Lpp_qK = sph_harm(M_Lpp, Lp, phi_qK, theta_qK)
        Y_Lppp_kp = sph_harm(M_Lppp, L, phi_kp, theta_kp)
        # Y_L_k = 1/np.sqrt(4*np.pi)
        # Y_Lp_qK = 1/np.sqrt(4*np.pi)
        # Y_Lpp_qK = 1/np.sqrt(4*np.pi)
        # Y_Lppp_kp = 1/np.sqrt(4*np.pi)
        
        # \delta U_{L,L'}(k, q-K/2)
        delta_U_partial_wave = delta_U_function.ev(k, qK)
        # \delta U^\dagger_{L'',L'''}(q-K/2, k')
        delta_U_dag_partial_wave = delta_U_dagger_function.ev(qK, kp)
        
        # Single-particle wave functions
        psi_alpha_k1 = psi(alpha, K_vector/2+k_vector, sigma_1, alpha.tau,
                           cg_table, phi_functions)
        psi_beta_k2 = psi(beta, K_vector/2-k_vector, sigma_2, beta.tau,
                          cg_table, phi_functions)
        psi_beta_k4 = psi(beta, K_vector/2-kp_vector, sigma_4, beta.tau,
                          cg_table, phi_functions)
        psi_alpha_k3 = psi(alpha, K_vector/2+kp_vector, sigma_3, alpha.tau,
                           cg_table, phi_functions)

        integrand += (
            1/2 * (1/2*2/np.pi)**2 * jacobian * sig_12_cg * sig_34_cg * sig_56_cg ** 2
            * tau_ab_cg ** 2 * tau_ttp_cg ** 2
            * lsj_cg * lpsj_cg * lppsj_cg * lpppsj_cg
            * lst_factor ** 2 * lpst_factor ** 2
            * Y_L_k * np.conj(Y_Lp_qK) * Y_Lpp_qK * np.conj(Y_Lppp_kp)
            * delta_U_partial_wave * delta_U_dag_partial_wave
            * psi_alpha_k1 * psi_beta_k2 * psi_beta_k4 * psi_alpha_k3
        )
    
    return integrand.real


def compute_delta_U2_term(
        q_array, cg_table, occ_states, channel, phi_functions, delta_U_function,
        delta_U_dagger_function
):
    """Compute the \delta U * n(q) * \delta U^\dagger term."""
        
    delta_U2_array = np.zeros_like(q_array)
    delta_U2_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    quantum_numbers = delta_U2_quantum_numbers(cg_table, occ_states, channel)
        
    # Set-up integration with vegas
    k_array, _ = momentum_mesh(15.0, 3.0, 120)
    k_max, k_min = max(k_array), min(k_array)
    k_limits = [k_min, k_max]  # Relative momenta up to 15^-1
    K_limits = [0, 3]  # C.o.M. momenta up to 3 fm^-1
    theta_limits = [0, np.pi]
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):

        integrand = functools.partial(
            delta_U2_term_integrand, q, channel, quantum_numbers, cg_table,
            phi_functions, delta_U_function, delta_U_dagger_function
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
          
    return delta_U2_array, delta_U2_errors


def compute_distribution(
        nucleus_name, Z, N, kvnn, channel, kmax, kmid, ntot, generator, lamb
):
    """Compute \delta U + \delta U^\dagger + \sum \delta U \delta U^\dagger.
    This terms should cancel and give zero when integrating over q.
    """
    
    # Compute table of Clebsch-Gordan coefficients
    cg_table = compute_clebsch_gordan_table(4)
    print("Done calculating Clebsch-Gordan table.\n")
    
    # Set-up single-particle states
    sp_basis = SingleParticleBasis(nucleus_name, Z, N, run_woodsaxon=False)
    phi_functions = get_sp_wave_functions(sp_basis)
    
    # Set-up \delta U and \delta U^\dagger functions
    delta_U_function = interpolate_delta_U(kvnn, channel, kmax, kmid, ntot,
                                           generator, lamb, False)
    delta_U_dagger_function = interpolate_delta_U(kvnn, channel, kmax, kmid,
                                                  ntot, generator, lamb, True)

    # Set momentum mesh
    q_array, q_weights = momentum_mesh(8.0, 2.0, 40)

    t0 = time.time()

    # Compute the I term
    I_array = compute_I_term(q_array, sp_basis.occ_states, cg_table,
                             phi_functions)
    
    # Compute \delta U + \delta U^\dagger term using vegas
    print("Beginning \delta U linear terms.\n")
    t1 = time.time()
    delta_U_array, delta_U_errors = compute_delta_U_term(
        q_array, cg_table, sp_basis.occ_states, channel, phi_functions,
        delta_U_function, delta_U_dagger_function
    )
    t2 = time.time()
    print(f"Done after {(t2-t1)/3600:.3f} hours.\n")

    # Compute \delta U \delta U^\dagger term using vegas
    print("Beginning \delta U \delta U^\dagger term.\n")
    t3 = time.time()
    delta_U2_array, delta_U2_errors = compute_delta_U2_term(
        q_array, cg_table, sp_basis.occ_states, channel, phi_functions,
        delta_U_function, delta_U_dagger_function
    )
    t4 = time.time()
    print(f"Done after {(t4-t3)/3600:.3f} hours.\n")
        
    t5 = time.time()
    print(f"Total time elapsed: {(t5-t0)/3600:.3f} hours.\n")
    
    # Combine each term for the total momentum distribution [fm^3]
    n_array = I_array + delta_U_array + delta_U2_array
    n_errors = np.sqrt(delta_U_errors**2 + delta_U2_errors**2)

    # Overall normalization
    overall_norm = 4*np.pi * np.sum(q_weights * q_array**2 * n_array)
    # I normalization
    I_norm = 4*np.pi * np.sum(q_weights * q_array**2 * I_array)
    # \delta U normalization
    delU_norm = 4*np.pi * np.sum(q_weights * q_array**2 * delta_U_array)
    # \delta U^2 normalization
    delU2_norm = 4*np.pi * np.sum(q_weights * q_array**2 * delta_U2_array)
    # Relative error in cancelation
    rel_err = abs(delU2_norm + delU_norm) / delU2_norm * 100
    
    print(f"Overall normalization = {overall_norm}.")
    print(f"I normalization = {I_norm}.")
    print(f"\delta U + \delta U^\dagger normalization = {delU_norm}.")
    print(f"\delta U \delta U^\dagger normalization = {delU2_norm}.")
    print(f"Relative error = {rel_err}%.")
    
    return q_array, q_weights, n_array, n_errors


if __name__ == '__main__':
    
    # Set nucleus, nucleon, partial wave channels, and potential
    nucleus_name, Z, N = 'He4', 2, 2
    kvnn = 111
    kmax, kmid, ntot = 15.0, 3.0, 120
    generator = 'Wegner'
    lamb = 1.35
    
    # channel = '1S0'
    # channel = '3S1-3S1'
    channel = '3D1-3D1'

    # He4
    q_array, q_weights, n_array, n_errors = compute_distribution(
        nucleus_name, Z, N, kvnn, channel, kmax, kmid, ntot, generator, lamb
    )
