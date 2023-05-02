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
quantum numbers as an additional integral and uses batch mode in vegas.

Last update: May 2, 2023

"""

# Python imports
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
from scripts.srg import SRG
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
                file_name = get_orbital_file_name(sp_state)
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


# UPDATE V6 VERSION TO ACCOUNT FOR VECTORS INSTEAD OF SCALARS
@vegas.batchintegrand
class DeltaUIntegrand:
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""
    
    def __init__(
            self, q, tau, quantum_numbers, cg_table, phi_functions,
            delta_U_functions, delta_U_dagger_functions
    ):
        self.q = q
        self.tau = tau
        self.quantum_numbers = quantum_numbers
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
        for d in self.quantum_numbers:
            
            # Unpack dictionary
            sigma_1 = d['sigma_1']
            sigma_2 = d['sigma_2']
            sigma = d['sigma']
            sigmap = d['sigmap']
            tau_1 = d['tau_1']
            tau_2 = d['tau_2']
            taup = d['taup']
            alpha = d['alpha']
            beta = d['beta']
            
            channel = d['channel']
            S = d['S']
            M_S = d['M_S']
            M_Sp = d['M_Sp']                                           
            L = d['L']
            M_L = d['M_L']
            Lp = d['Lp']
            M_Lp = d['M_Lp']
            J = d['J']
            M_J = d['M_J']                                            
            T = d['T']
            M_T = d['M_T']
            
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
            
            # < k (L S) J T | \delta U | |q-K/2| (L' S) J T >
            delta_U_partial_wave = self.delta_U_functions[channel].ev(k, qK)
            
            # < |q-K/2| (L' S) J T | \delta U^\dagger | k (L S) J T >
            delta_U_dag_partial_wave = (
                self.delta_U_dagger_functions[channel].ev(qK, k)
            )
            
            # \psi_\alpha(K/2+k; \sigma_1, \tau_1)
            psi_alpha_1 = psi(alpha, k1, theta_k1, phi_k1, sigma_1, tau_1,
                              self.cg_table, self.phi_functions)
            
            # \psi_\beta(K/2-k; \sigma_2, \tau_2)
            psi_beta_2 = psi(beta, k2, theta_k2, phi_k2, sigma_2, tau_2,
                             self.cg_table, self.phi_functions)
            
            # \psi_\alpha(q; \sigma, \tau)
            psi_alpha_q = psi(alpha, self.q, 0, 0, sigma, self.tau,
                              self.cg_table, self.phi_functions)
            
            # \psi_\beta(K-q; \sigma', \tau')
            psi_beta_Kq = psi(beta, Kq, theta_Kq, phi_Kq, sigmap, taup,
                              self.cg_table, self.phi_functions)
            
            # \psi_\beta(q; \sigma, \tau)
            psi_beta_q = psi(beta, self.q, 0, 0, sigma, self.tau, self.cg_table,
                             self.phi_functions)
            
            # \psi_\alpha(K-q; \sigma', \tau')
            psi_alpha_Kq = psi(alpha, Kq, theta_Kq, phi_Kq, sigmap, taup,
                               self.cg_table, self.phi_functions)
            
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
            integrand += 1/2 * 1/2 * 2/np.pi * jacobian * (
                delta_U_term + delta_U_dag_term
            )

        return integrand.real
    
    
    def get_quantum_number_arrays(self):
        return None

    
@vegas.batchintegrand
class DeltaUIntegrand:
    """Evaluate the integrand of the \delta U + \delta U^\dagger terms."""
    
    
    def __init__(
        self, tau, q, nucleus_name, Z, N, sp_quantum_numbers, cg_table,
        channels, kvnn, kmax, kmid, ntot, lamb
    ):
        
        # Call Woods-Saxon class
        woods_saxon = WoodsSaxon(nucleus_name, Z, N, run_woodsaxon=False)
        
        # Set occupied states
        self.occupied_states = self.set_occupied_states(woods_saxon)

        # Set single-particle wave functions
        self.sp_wave_functions = self.get_sp_wave_functions(
            woods_saxon, nucleus_name, Z, N, 10.0, 2.0, 120
        )
        
        channels_list = []
        for channel in channels:
            channels_list.append(PartialWaveChannel(channel))
        self.channels = channels_list
        
        # Set Clebsch-Gordan table
        self.cg_table = cg_table
        
        # Vectorize Clebsch-Gordan coefficient function
        self.cg_func = np.vectorize(self.get_cg)
        
        # Set \delta U and \delta U^\dagger functions
        self.delta_U_functions, self.delta_U_dagger_functions = (
            self.get_delta_U_functions(kvnn, kmax, kmid, ntot, 'Wegner', lamb)
        )
        
        # Set instance attributes
        self.tau = tau
        self.q = q
        self.sp_quantum_numbers = sp_quantum_numbers
        self.N_sp = len(sp_quantum_numbers)
        
        # Use vectorized s.p. quantum numbers function
        self.vectorized_get_sp_quantum_numbers = np.vectorize(
            self.get_sp_quantum_numbers)
        
        # Use vectorized s.p. wave function
        # self.psi = np.vectorize(self.compute_psi, otypes=[complex])
        self.psi = np.vectorize(self.compute_psi)
        # GIVES ERROR
        # self.psi = self.compute_psi

        
        # Use vectorized \delta U matrix element function
        # self.delta_U_matrix_element = np.vectorize(
        #     self.compute_delta_U_matrix_element, otypes=[complex])
        self.delta_U_matrix_element = np.vectorize(
            self.compute_delta_U_matrix_element)
        # GIVES ERROR
        # self.delta_U_matrix_element = self.compute_delta_U_matrix_element
        
        
    def set_occupied_states(self, woods_saxon):
        """Create a N x 5 array of the occupied Woods-Saxon orbitals described
        by the quantum numbers n, l, j, m_j, and \tau. Here N corresponds to the
        number of occupied states.
        """
        
        N = np.size(woods_saxon.occ_states)
        occupied_states = np.zeros((N, 5))
        for i, sp_state in enumerate(woods_saxon.occ_states):
            occupied_states[i, 0] = sp_state.n
            occupied_states[i, 1] = sp_state.l
            occupied_states[i, 2] = sp_state.j
            occupied_states[i, 3] = sp_state.m_j
            occupied_states[i, 4] = sp_state.tau
            
        return occupied_states


    def get_sp_quantum_numbers(self, x):
        """Given an integration variable x, return a set of four spin
        projections and two s.p. states \alpha and \beta.
        """
        
        index = np.floor(x * self.N_sp).astype(int)
        return self.sp_quantum_numbers[index]

    
    def __call__(self, x_array):
        
        # Get s.p. quantum numbers
        (sigma_1, sigma_2, sigma, sigmap, n_alpha, l_alpha, j_alpha, m_j_alpha,
        tau_alpha, n_beta, l_beta, j_beta, m_j_beta, tau_beta) = (
            self.vectorized_get_sp_quantum_numbers(x_array[:, 0])
        )

        # Relative momenta k
        k = x_array[:, 1]
        theta_k = x_array[:, 2]
        phi_k = x_array[:, 3]
        k_vector = build_vector(k, theta_k, phi_k)
        
        # C.o.M. momenta K
        K = x_array[:, 4]
        theta_K = x_array[:, 5]
        phi_K = x_array[:, 6]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Choose z-axis to be along q_vector
        q_vector = np.zeros_like(k_vector)
        q_vector[-1, :] = self.q
    
        # Calculate vector q-K/2
        qK_vector = q_vector - K_vector/2
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
        
        # Calculate the Jacobian determinant
        jacobian = k ** 2 * np.sin(theta_k) * K ** 2 * np.sin(theta_K)
        
        # Get s.p. wave functions
        k1, theta_1, phi_1 = get_vector_components(K_vector/2 + k_vector) 
        k2, theta_2, phi_2 = get_vector_components(K_vector/2 - k_vector)
        Kq, theta_Kq, phi_Kq = get_vector_components(K_vector - q_vector)
        
        psi_alpha_k1 = self.psi(
            n_alpha, l_alpha, j_alpha, m_j_alpha, tau_alpha, k1, theta_1, phi_1,
            sigma_1
        )
        psi_beta_k2 = self.psi(n_beta, l_beta, j_beta, m_j_beta, tau_beta, k2,
                                theta_2, phi_2, sigma_2)
        psi_alpha_q = self.psi(n_alpha, l_alpha, j_alpha, m_j_alpha, tau_alpha,
                                self.q, 0.0, 0.0, sigma)
        psi_beta_Kq = self.psi(n_beta, l_beta, j_beta, m_j_beta, tau_beta, Kq,
                                theta_Kq, phi_Kq, sigmap)
        psi_alpha_Kq = self.psi(
            n_alpha, l_alpha, j_alpha, m_j_alpha, tau_alpha, Kq, theta_Kq,
            phi_Kq, sigmap
        )
        psi_beta_q = self.psi(
            n_beta, l_beta, j_beta, m_j_beta, tau_beta, self.q, 0.0, 0.0, sigma
        )
        
        # # TESTING
        # psi_alpha_k1_r, psi_alpha_k1_i = self.psi(
        #     n_alpha, l_alpha, j_alpha, m_j_alpha, tau_alpha, k1, theta_1, phi_1,
        #     sigma_1
        # )
        # psi_alpha_k1 = psi_alpha_k1_r + 1j * psi_alpha_k1_i
        
        # psi_beta_k2_r, psi_beta_k2_i = self.psi(
        #     n_beta, l_beta, j_beta, m_j_beta, tau_beta, k2, theta_2, phi_2,
        #     sigma_2
        # )
        # psi_beta_k2 = psi_beta_k2_r + 1j * psi_beta_k2_i
        
        # psi_alpha_q_r, psi_alpha_q_i = self.psi(
        #     n_alpha, l_alpha, j_alpha, m_j_alpha, tau_alpha, self.q, 0.0, 0.0,
        #     sigma
        # )
        # psi_alpha_q = psi_alpha_q_r + 1j * psi_alpha_q_i
        
        # psi_beta_Kq_r, psi_beta_Kq_i = self.psi(
        #     n_beta, l_beta, j_beta, m_j_beta, tau_beta, Kq, theta_Kq, phi_Kq,
        #     sigmap
        # )
        # psi_beta_Kq = psi_beta_Kq_r + 1j * psi_beta_Kq_i
        
        # psi_alpha_Kq_r, psi_alpha_Kq_i = self.psi(
        #     n_alpha, l_alpha, j_alpha, m_j_alpha, tau_alpha, Kq, theta_Kq,
        #     phi_Kq, sigmap
        # )
        # psi_alpha_Kq = psi_alpha_Kq_r + 1j * psi_alpha_Kq_i
        
        # psi_beta_q_r, psi_beta_q_i = self.psi(
        #     n_beta, l_beta, j_beta, m_j_beta, tau_beta, self.q, 0.0, 0.0, sigma
        # )
        # psi_beta_q = psi_beta_q_r + 1j * psi_beta_q_i

        # Calculate plane-wave matrix elements of \delta U and \delta U^\dagger
        delta_U_plane_wave_abab = self.delta_U_matrix_element(
            k, theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta, qK,
            theta_qK, phi_qK, sigma, tau_alpha, sigmap, tau_beta
        )
        delta_U_plane_wave_abba = self.delta_U_matrix_element(
            k, theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta, qK,
            theta_qK, phi_qK, sigma, tau_beta, sigmap, tau_alpha
        )
        delta_U_dag_plane_wave_abab = self.delta_U_matrix_element(
            qK, theta_qK, phi_qK, sigma, tau_alpha, sigmap, tau_beta, k,
            theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta,
            hermitian_conjugate=True
        )
        delta_U_dag_plane_wave_baab = self.delta_U_matrix_element(
            qK, theta_qK, phi_qK, sigma, tau_beta, sigmap, tau_alpha, k,
            theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta,
            hermitian_conjugate=True
        )
        
        # # TESTING
        # delta_U_abab_r, delta_U_abab_i = self.delta_U_matrix_element(
        #     k, theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta, qK,
        #     theta_qK, phi_qK, sigma, tau_alpha, sigmap, tau_beta
        # )
        # delta_U_plane_wave_abab = delta_U_abab_r + 1j * delta_U_abab_i
        
        # delta_U_abba_r, delta_U_abba_i = self.delta_U_matrix_element(
        #     k, theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta, qK,
        #     theta_qK, phi_qK, sigma, tau_beta, sigmap, tau_alpha
        # )
        # delta_U_plane_wave_abba = delta_U_abba_r + 1j * delta_U_abba_i
        
        # delta_U_dag_abab_r, delta_U_dag_abab_i = self.delta_U_matrix_element(
        #     qK, theta_qK, phi_qK, sigma, tau_alpha, sigmap, tau_beta, k,
        #     theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta,
        #     hermitian_conjugate=True
        # )
        # delta_U_dag_plane_wave_abab = delta_U_dag_abab_r + 1j * delta_U_dag_abab_i
        
        # delta_U_dag_baab_r, delta_U_dag_baab_i = self.delta_U_matrix_element(
        #     qK, theta_qK, phi_qK, sigma, tau_beta, sigmap, tau_alpha, k,
        #     theta_k, phi_k, sigma_1, tau_alpha, sigma_2, tau_beta,
        #     hermitian_conjugate=True
        # )
        # delta_U_dag_plane_wave_baab = delta_U_dag_baab_r + 1j * delta_U_dag_baab_i
        
        # Isospin Kronecker \delta's
        del_alpha = (tau_alpha == self.tau).astype(int)
        del_beta = (tau_beta == self.tau).astype(int)
            
        # Return integrand in plane-wave basis
        integrand = 1/2 * self.N_sp * jacobian * (
            np.conj(psi_alpha_k1) * np.conj(psi_beta_k2)
            * (
                del_alpha * psi_alpha_q * psi_beta_Kq * delta_U_plane_wave_abab
                - del_beta * psi_beta_q * psi_alpha_Kq * delta_U_plane_wave_abba
            )
            + psi_alpha_k1 * psi_beta_k2
            * (
                del_alpha * np.conj(psi_alpha_q) * np.conj(psi_beta_Kq)
                * delta_U_dag_plane_wave_abab
                - del_beta * np.conj(psi_beta_q) * np.conj(psi_alpha_Kq)
                * delta_U_dag_plane_wave_baab
            )
        )
        
        return integrand.real
    

# Use \delta U version to write this class
@vegas.batchintegrand
class DeltaU2Integrand:
    
    
    def __init__(self, tau, q):
        
        # Set instance attributes
        self.tau = tau
        self.q = q
        
    
    def __call__(self):
        pass
    
    
    def get_quantum_number_arrays(self):
        return None
    
        
def get_orbital_file_name(n, l, j, tau):
    """Returns the file name of the Woods-Saxon orbital."""
        
    # Proton
    if tau == 1/2:
        file_name = f"p.n{int(n-1)}.l{int(l)}.j{int(2*j)}.orb"
    # Neutron
    elif tau == -1/2:
        file_name = f"n.n{int(n-1)}.l{int(l)}.j{int(2*j)}.orb"
    
    return file_name


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
            j_3_array = np.arange(abs(j_1-j_2), j_1+j_2+1/2)
            for j_3 in j_3_array:
                for m_1 in np.arange(-j_1, j_1+1, 1):
                    for m_2 in np.arange(-j_2, j_2+1, 1):
                        m_3 = m_1 + m_2
                        if abs(m_3) <= j_3:
                            cg_table[(j_1,m_1,j_2,m_2,j_3,m_3)] = float(
                                CG(j_1,m_1,j_2,m_2,j_3,m_3).doit()
                            )
                                    
    return cg_table


def cg_func(j1, m1, j2, m2, j3, m3, cg_table):
    """Clebsch-Gordan coefficient < j1 m1 j2 m2 | j3 m3 >."""
        
    if (m1 + m2 == m3 and np.absolute(m1) <= j1 and np.absolute(m2) <= j2
        and np.absolute(m3) <= j3):

        try:
            return cg_table[(j1, m1, j2, m2, j3, m3)]
                
        except KeyError:
            return 0
                
    else:
        return 0
    
    
def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return int(x == y)

        
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

    k_vector = np.array([k * np.sin(theta) * np.cos(phi),
                         k * np.sin(theta) * np.sin(phi),
                         k * np.cos(theta)])

    return k_vector


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

    k = la.norm(k_vector, axis=0)
    theta = np.arccos(k_vector[2]/k)
    phi = np.arctan2(k_vector[1], k_vector[0])

    return k, theta, phi


def get_sp_wave_functions(sp_basis, kmax, kmid, ntot):
    """Set interpolating functions for s.p. wave functions \phi."""
    
    occ_states = sp_basis.occ_states

    phi_functions = {}
    for sp_state in occ_states: 
        file_name = get_orbital_file_name(sp_state)
        phi_functions[file_name] = sp_basis.get_wf_kspace(
            sp_state, kmax, kmid, ntot, interpolate=True)
            
    return phi_functions


def psi(sp_state, k, theta, phi, sigma, tau, cg_table, phi_functions):
    """Single-particle wave function."""

    # Calculate \phi_\alpha(q)
    phi_sp_wf = phi_functions[get_orbital_file_name(sp_state)](k)
    
    # Calculate spinor spherical harmonic
    Y_jml = spinor_spherical_harmonic(sp_state, theta, phi, sigma, cg_table)
    
    # Isospinor indexed by \tau \chi_{m_t}(\tau)
    chi_tau = kronecker_delta(tau, sp_state.m_t)

    return phi_sp_wf * Y_jml * chi_tau


def spinor_spherical_harmonic(sp_state, theta, phi, sigma, cg_table):
    """Spinor spherical harmonic for a s.p. state described by the quantum
    numbers j, m_j, l, and s=1/2.
    """
    
    Y_jml = 0+0j

    # Spinor indexed by \sigma \eta_{m_s}^(\sigma) = \delta_{m_s, \sigma}
    m_s = sigma
    
    # m_l must be fixed since m_j and m_s are determined
    m_l = sp_state.m_j - m_s
        
    # Check that |m_l| <= l
    if np.abs(m_l) <= sp_state.l:
        
        # Clebsch-Gordan coefficient for l-s coupling
        cg = cg_table[(sp_state.l, m_l, 1/2, m_s, sp_state.j, sp_state.m_j)]
        
        # Spherical harmonic
        Y_lm = sph_harm(m_l, sp_state.l, phi, theta)
            
        Y_jml += cg * Y_lm
            
    return Y_jml


def interpolate_delta_U(
        kvnn, channel, kmax, kmid, ntot, generator, lamb,
        hermitian_conjugate=False
):
    """Interpolate \delta U(k, k') for the given channel."""
    
    # Set channel argument to be compatible with potential functions
    channel_str = channel.channel
    if channel_str[:3] == '3D1':
        channel_arg = '3S1'
    elif channel_str[:3] == '3F2':
        channel_arg = '3P2'
    elif channel_str[:3] == '3G3':
        channel_arg == '3D3'
    else:
        channel_arg = channel_str[:3]

    # Get momentum mesh
    k_array, k_weights = momentum_mesh(kmax, kmid, ntot)
    
    srg = SRG(kvnn, channel_arg, kmax, kmid, ntot, generator)
    U_matrix_weights = srg.load_srg_transformation(lamb)

    # Calculate \delta U = U - I
    I_matrix_weights = np.eye(len(U_matrix_weights))
    if hermitian_conjugate:
        delU_matrix_weights = (U_matrix_weights - I_matrix_weights).T
    else:
        delU_matrix_weights = U_matrix_weights - I_matrix_weights

    # Get specific sub-block if coupled-channel
    if channel_str in ['3S1-3D1', '3P2-3F2', '3D3-3G3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[:ntot,ntot:]
        )
    elif channel_str in ['3D1-3S1', '3F2-3P2', '3G3-3D3']:
        delU_matrix = unattach_weights_from_matrix(
            k_array, k_weights, delU_matrix_weights[ntot:,:ntot]
        )
    elif channel_str in ['3D1-3D1', '3F2-3F2', '3G3-3G3']:
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
    for channel in self.channels:
        channel_str = channel.channel
        delta_U_functions[channel_str] = self.interpolate_delta_U(
            channel, kvnn, kmax, kmid, ntot, generator, lamb
        )
        delta_U_dagger_functions[channel_str] = self.interpolate_delta_U(
            channel, kvnn, kmax, kmid, ntot, generator, lamb,
            hermitian_conjugate=True
        )
        
    return delta_U_functions, delta_U_dagger_functions


def delta_U_quantum_numbers(tau, occ_states, cg_table, channels):
    """Returns a list of every combination of quantum numbers in the \delta U
    linear terms that gives a non-zero product of Clebsch-Gordan coefficients.
    """
    
    spins = np.array([1/2, -1/2])
    combinations = []
    
    # Partial wave channels
    for channel in channels:
        
        L, Lp, S, J = get_channel_quantum_numbers(channel)
        T = get_total_isospin(L, S)
        
        # 1 - (-1)^(L+S+T) factor
        lst_factor = 1-(-1)**(L+S+T)
        # 1 - (-1)^(L'+S+T) factor
        lpst_factor = 1-(-1)**(Lp+S+T)

        # Spin projections \sigma_1 and \sigma_2
        for sigma_1 in spins:
            for sigma_2 in spins:
                
                # < \sigma_1 \sigma_2 | S M_S >
                M_S = sigma_1 + sigma_2
                spin_12_cg = cg_func(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                     cg_table)
                
                # Spin projections \sigma and \sigma'
                for sigma in spins:
                    for sigmap in spins:
                        
                        # < S M_S' | \sigma \sigma' >
                        M_Sp = sigma + sigmap
                        spin_ssp_cg = cg_func(1/2, sigma, 1/2, sigmap, S, M_Sp,
                                              cg_table)
                        
                        # Isospin projections \tau_1 and \tau_2
                        for tau_1 in spins:
                            for tau_2 in spins:
                                
                                # < \tau_1 \tau_2 | T M_T >
                                M_T = tau_1 + tau_2
                                isospin_12_cg = cg_func(1/2, tau_1, 1/2, tau_2,
                                                        T, M_T, cg_table)
                                
                                # Isospin projection \tau'
                                for taup in spins:
                                    
                                    # < T M_T | \tau \tau' >
                                    M_Tp = tau + taup
                                    if M_Tp == M_T:
                                        isospin_ttp_cg = cg_func(
                                            1/2, tau, 1/2, taup, T, M_T,
                                            cg_table
                                        )
                                    else:
                                        isospin_ttp_cg = 0
                                    
                                    for M_J in np.arange(-J, J+1):
                                        
                                        # Channel orbital angular momentum
                                        # projection M_L  
                                        M_L = M_J - M_S
                                        M_Lp = M_J - M_Sp
                                        
                                        # < L M_L S M_S | J M_J >
                                        lsj_cg = cg_func(L, M_L, S, M_S, J, M_J,
                                                         cg_table)
                                        
                                        # < J M_J | L' M_L' S M_S' >
                                        lpsj_cg = cg_func(Lp, M_Lp, S, M_Sp, J,
                                                          M_J, cg_table)
                                        
                                        # Single-particle states
                                        for alpha in occ_states:
                                            
                                            # \delta_{\tau_1, \tau_\alpha}
                                            chi_1a = kronecker_delta(tau_1,
                                                                     alpha.m_t)
                                            
                                            # \delta_{\tau, \tau_\alpha}
                                            chi_ta = kronecker_delta(tau,
                                                                     alpha.m_t)
                                            
                                            # \delta_{\tau', \tau_\alpha}
                                            chi_tpa = kronecker_delta(taup,
                                                                      alpha.m_t)
                                            
                                            # l_\alpha \sigma_1 CG
                                            alpha_1_cg = cg_func(
                                                alpha.l, alpha.m_j-sigma_1, 1/2,
                                                sigma_1, alpha.j, alpha.m_j,
                                                cg_table
                                            )
                                            
                                            # l_\alpha \sigma CG
                                            alpha_s_cg = cg_func(
                                                alpha.l, alpha.m_j-sigma, 1/2,
                                                sigma, alpha.j, alpha.m_j,
                                                cg_table
                                            )
                                            
                                            # l_\alpha \sigma' CG
                                            alpha_sp_cg = cg_func(
                                                alpha.l, alpha.m_j-sigmap, 1/2,
                                                sigmap, alpha.j, alpha.m_j,
                                                cg_table
                                            )

                                            for beta in occ_states:
                                            
                                                # \delta_{\tau_2, \tau_\beta}
                                                chi_2b = kronecker_delta(
                                                    tau_2, beta.m_t)
                                                
                                                # \delta_{\tau', \tau_\beta}
                                                chi_tpb = kronecker_delta(
                                                    taup, beta.m_t)
                                                
                                                # \delta_{\tau, \tau_\beta}
                                                chi_tb = kronecker_delta(
                                                    tau, beta.m_t)
                                                
                                                # l_\beta \sigma_2 CG
                                                beta_2_cg = cg_func(
                                                    beta.l, beta.m_j-sigma_2,
                                                    1/2, sigma_2, beta.j,
                                                    beta.m_j, cg_table
                                                )
                                                
                                                # l_\beta \sigma' CG
                                                beta_sp_cg = cg_func(
                                                    beta.l, beta.m_j-sigmap,
                                                    1/2, sigmap, beta.j,
                                                    beta.m_j, cg_table
                                                )
                                                
                                                # l_\beta \sigma CG
                                                beta_s_cg = cg_func(
                                                    beta.l, beta.m_j-sigma, 1/2,
                                                    sigma, beta.j, beta.m_j,
                                                    cg_table
                                                )
                                                
                                                # Calculate the product of all
                                                # CG's and factors
                                                product = (
                                                    spin_12_cg
                                                    * spin_ssp_cg
                                                    * isospin_12_cg
                                                    * isospin_ttp_cg
                                                    * lsj_cg
                                                    * lpsj_cg
                                                    * lst_factor
                                                    * lpst_factor
                                                    * chi_1a * alpha_1_cg
                                                    * chi_2b * beta_2_cg
                                                    * (
                                                        chi_ta * alpha_s_cg
                                                        * chi_tpb * beta_sp_cg
                                                        - chi_tpa * alpha_sp_cg
                                                        * chi_tb * beta_s_cg
                                                    )
                                                )
                                                
                                                # Add this set of quantum
                                                # numbers to the list if the
                                                # product is nonzero
                                                if product != 0:
                                                                                                
                                                    d = {
                                                        'sigma_1': sigma_1,
                                                        'sigma_2': sigma_2,
                                                        'sigma': sigma,
                                                        'sigmap': sigmap,
                                                        'tau_1': tau_1,
                                                        'tau_2': tau_2,
                                                        'taup': taup,
                                                        'alpha': alpha,
                                                        'beta': beta,
                                                        'channel': channel,
                                                        'S': S,
                                                        'M_S': M_S,
                                                        'M_Sp': M_Sp,
                                                        'L': L,
                                                        'M_L': M_L,
                                                        'Lp': Lp,
                                                        'M_Lp': M_Lp,
                                                        'J': J,
                                                        'M_J': M_J,
                                                        'T': T,
                                                        'M_T': M_T
                                                    }
                                                                                                
                                                    combinations.append(d)
                                    
    return combinations


def delta_U2_quantum_numbers(tau, occ_states, cg_table, channels):
    """Returns a list of every combination of quantum numbers in the \delta U^2
    term that gives a non-zero product of Clebsch-Gordan coefficients.
    """
    
    spins = np.array([1/2, -1/2])
    combinations = []
    
    # Partial wave channel for \delta U
    for channel_1 in channels:
        
        L, Lp, S, J = get_channel_quantum_numbers(channel_1)
        T = get_total_isospin(L, S)
        
        # 1 - (-1)^(L+S+T) factor
        lst_factor = 1-(-1)**(L+S+T)
        # 1 - (-1)^(L'+S+T) factor
        lpst_factor = 1-(-1)**(Lp+S+T)
        
        # Partial wave channel for \delta U^\dagger
        for channel_2 in channels:
            
            Lpp, Lppp, Sp, Jp = get_channel_quantum_numbers(channel_2)
            
            # Evaluate Kronecker \delta function in S and S'
            if Sp == S:
                
                Tp = get_total_isospin(Lpp, S)
            
                # 1 - (-1)^(L''+S+T') factor
                lppstp_factor = 1-(-1)**(Lpp+S+Tp)
                # 1 - (-1)^(L'''+S+T') factor
                lpppstp_factor = 1-(-1)**(Lppp+S+Tp)

                # Spin projections \sigma_1 and \sigma_2
                for sigma_1 in spins:
                    for sigma_2 in spins:

                        # < \sigma_1 \sigma_2 | S M_S >
                        M_S = sigma_1 + sigma_2
                        spin_12_cg = cg_func(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                             cg_table)
                        
                        # Spin projections \sigma_3 and \sigma_4
                        for sigma_3 in spins:
                            for sigma_4 in spins:
                                
                                # < S M_S''' | \sigma_3 \sigma_4 >
                                M_Sppp = sigma_3 + sigma_4
                                spin_34_cg = cg_func(1/2, sigma_3, 1/2, sigma_4,
                                                     S, M_Sppp, cg_table)
                                
                                # Isospin projections \tau_1 and \tau_2
                                for tau_1 in spins:
                                    for tau_2 in spins:
                                        
                                        # < \tau_1 \tau_2 | T M_T >
                                        M_T = tau_1 + tau_2
                                        isospin_12_cg = cg_func(
                                            1/2, tau_1, 1/2, tau_2, T, M_T,
                                            cg_table
                                        )
                                        
                                        # Isospin projections \tau_3 and \tau_4
                                        for tau_3 in spins:
                                            for tau_4 in spins:
                                                
                                                # < T' M_T' | \tau_3 \tau_4 >
                                                M_Tp = tau_3 + tau_4
                                                isospin_34_cg = cg_func(
                                                    1/2, tau_3, 1/2, tau_4, Tp,
                                                    M_Tp, cg_table
                                                )
                                        
                                                # Isospin projection \tau'
                                                for taup in spins:
                                            
                                                    # < T M_T | \tau \tau' >
                                                    isospin_ttp_cg = cg_func(
                                                        1/2, tau, 1/2, taup, T,
                                                        M_T, cg_table
                                                    )
                                                    
                                                    # < \tau \tau' | T' M_T' >
                                                    isospin_ttpp_cg = cg_func(
                                                        1/2, tau, 1/2, taup, Tp,
                                                        M_Tp, cg_table
                                                    )
                                                    
                                                    for M_J in np.arange(-J, J+1):
                                                        
                                                        M_L = M_J - M_S
                                                        
                                                        # L S J CG
                                                        lsj_cg = cg_func(
                                                            L, M_L, S, M_S,
                                                            J, M_J, cg_table
                                                        )
                                                        
                                                        # Sum over M_S'
                                                        for M_Sp in np.arange(-S, S+1):
                                                        
                                                            M_Lp = M_J - M_Sp
                                                        
                                                            # L' S J CG
                                                            lpsj_cg = cg_func(
                                                                Lp, M_Lp, S,
                                                                M_Sp, J, M_J,
                                                                cg_table
                                                            )
                                                        
                                                            for M_Jp in np.arange(-Jp, Jp+1):
                                                            
                                                                M_Lpp = M_Jp - M_Sp
                                                                M_Lppp = M_Jp - M_Sppp
                                                                
                                                                # L'' S J' CG
                                                                lppsjp_cg = cg_func(
                                                                    Lpp, M_Lpp, S, M_Sp,
                                                                    Jp, M_Jp, cg_table
                                                                )
                                                                
                                                                # L''' S J' CG
                                                                lpppsjp_cg = cg_func(
                                                                    Lppp, M_Lppp, S, M_Sppp,
                                                                    Jp, M_Jp, cg_table
                                                                )
                                                                
                                                                # Single-particle states
                                                                for alpha in occ_states:
                                                                    
                                                                    # \delta_{\tau_1, \tau_\alpha}
                                                                    chi_1a = kronecker_delta(tau_1, alpha.m_t)
                                                                    
                                                                    # \delta_{\tau_3, \tau_\alpha}
                                                                    chi_3a = kronecker_delta(tau_3, alpha.m_t)
                                                                    
                                                                    # \delta_{\tau_4, \tau_\alpha}
                                                                    chi_4a = kronecker_delta(tau_4, alpha.m_t)
                                                        
                                                                    # l_\alpha \sigma_1 CG
                                                                    alpha_1_cg = cg_func(
                                                                        alpha.l, alpha.m_j-sigma_1, 1/2,
                                                                        sigma_1, alpha.j, alpha.m_j, cg_table
                                                                    )
                                                        
                                                                    # l_\alpha \sigma_3 CG
                                                                    alpha_3_cg = cg_func(
                                                                        alpha.l, alpha.m_j-sigma_3, 1/2,
                                                                        sigma_3, alpha.j, alpha.m_j, cg_table
                                                                    )
                                                                    
                                                                    # l_\alpha \sigma_4 CG
                                                                    alpha_4_cg = cg_func(
                                                                        alpha.l, alpha.m_j-sigma_4, 1/2,
                                                                        sigma_4, alpha.j, alpha.m_j, cg_table
                                                                    )
                                                                    
                                                                    for beta in occ_states:
                                                                        
                                                                        # \delta_{\tau_2, \tau_\beta}
                                                                        chi_2b = kronecker_delta(tau_2, beta.m_t)
                                                                        
                                                                        # \delta_{\tau_3, \tau_\beta}
                                                                        chi_3b = kronecker_delta(tau_3, beta.m_t)
                                                                        
                                                                        # \delta_{\tau_4, \tau_\beta}
                                                                        chi_4b = kronecker_delta(tau_4, beta.m_t)
                                                                        
                                                                        # l_\beta \sigma_2 CG
                                                                        beta_2_cg = cg_func(
                                                                            beta.l, beta.m_j-sigma_2, 1/2,
                                                                            sigma_2, beta.j, beta.m_j, cg_table
                                                                        )
                                                                        
                                                                        # l_\beta \sigma_3 CG
                                                                        beta_3_cg = cg_func(
                                                                            beta.l, beta.m_j-sigma_3, 1/2,
                                                                            sigma_3, beta.j, beta.m_j, cg_table
                                                                        )
                                                                        
                                                                        # l_\beta \sigma_4 CG
                                                                        beta_4_cg = cg_func(
                                                                            beta.l, beta.m_j-sigma_4, 1/2,
                                                                            sigma_4, beta.j, beta.m_j, cg_table
                                                                        )
                                        
                                                                        # Calculate the product of all
                                                                        # CG's and factors
                                                                        product = (
                                                                            spin_12_cg
                                                                            * spin_34_cg
                                                                            * isospin_12_cg
                                                                            * isospin_ttp_cg
                                                                            * isospin_ttpp_cg
                                                                            * isospin_34_cg
                                                                            * lsj_cg
                                                                            * lpsj_cg
                                                                            * lppsjp_cg
                                                                            * lpppsjp_cg
                                                                            * lst_factor
                                                                            * lpst_factor
                                                                            * lppstp_factor
                                                                            * lpppstp_factor
                                                                            * chi_1a * alpha_1_cg
                                                                            * chi_2b * beta_2_cg
                                                                            * (
                                                                                chi_3a * alpha_3_cg
                                                                                * chi_4b * beta_4_cg
                                                                                - chi_4a * alpha_4_cg
                                                                                * chi_3b * beta_3_cg
                                                                            )  
                                                                        )

                                                                        # Add this set of quantum
                                                                        # numbers to the list if the
                                                                        # product is nonzero
                                                                        if product != 0:
                                                                                                
                                                                            d = {
                                                                                'sigma_1': sigma_1,
                                                                                'sigma_2': sigma_2,
                                                                                'sigma_3': sigma_3,
                                                                                'sigma_4': sigma_4,
                                                                                'tau_1': tau_1,
                                                                                'tau_2': tau_2,
                                                                                'tau_3': tau_3,
                                                                                'tau_4': tau_4,
                                                                                'taup': taup,
                                                                                'alpha': alpha,
                                                                                'beta': beta,
                                                                                'channel_1': channel_1,
                                                                                'S': S,
                                                                                'M_S': M_S,
                                                                                'M_Sp': M_Sp,
                                                                                'L': L,
                                                                                'M_L': M_L,
                                                                                'Lp': Lp,
                                                                                'M_Lp': M_Lp,
                                                                                'J': J,
                                                                                'M_J': M_J,
                                                                                'T': T,
                                                                                'M_T': M_T,
                                                                                'channel_2': channel_2,
                                                                                'M_Sppp': M_Sppp,
                                                                                'Lpp': Lpp,
                                                                                'M_Lpp': M_Lpp,
                                                                                'Lppp': Lppp,
                                                                                'M_Lppp': M_Lppp,
                                                                                'Jp': Jp,
                                                                                'M_Jp': M_Jp,
                                                                                'Tp': Tp,
                                                                                'M_Tp': M_Tp
                                                                            }
                                                                                                
                                                                            combinations.append(d)
                                    
    return combinations


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
                psi_alpha_array[i] = psi(alpha, q, 0, 0, sigma, tau, cg_table,
                                         phi_functions)
                    
            I_array += abs(psi_alpha_array)**2
                    
    return I_array


def compute_delta_U_term(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions,
):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    quantum_numbers = delta_U_quantum_numbers(tau, occ_states, cg_table,
                                              channels)
    print("Done organizing \delta U quantum numbers.\n")
    delU_Ntot = len(quantum_numbers)
        
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
        
        integrand = delta_U_term_integrand(
            q, tau, quantum_numbers, cg_table, phi_functions, delta_U_functions,
            delta_U_dagger_functions
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


def compute_delta_U2_term(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions
):
    """Compute the \delta U * n(q) * \delta U^\dagger term."""
        
    delta_U2_array = np.zeros_like(q_array)
    delta_U2_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    quantum_numbers = delta_U2_quantum_numbers(tau, occ_states, cg_table,
                                               channels)
    print("Done organizing \delta U \delta U^\dagger quantum numbers.\n")
    delU2_Ntot = len(quantum_numbers)
        
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
                              k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits,
                              sum_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = delta_U2_term_integrand(
            q, tau, quantum_numbers, cg_table, phi_functions, delta_U_functions,
            delta_U_dagger_functions
        )
        
        # Train the integrator
        integ(integrand, nitn=5, neval=5e4)
        # Final result
        result = integ(integrand, nitn=10, neval=5e4)

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
    # q_array, q_weights = momentum_mesh(8.0, 2.0, 40)
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
        print(f"Done after {(t2-t1)/60:.3f} minutes.\n")
        
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
        print(f"Done after {(t4-t3)/60:.3f} minutes.\n")
        
        # # TESTING
        # delta_U2_array = np.zeros_like(I_array)
        # delta_U2_errors = np.zeros_like(I_array)
        
        t5 = time.time()
        print(f"Total time elapsed: {(t5-t0)/60:.3f} minutes.\n")
    
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
    # kvnn, kmax, kmid, ntot = 6, 30.0, 4.0, 120  # AV18
    # kvnn, kmax, kmid, ntot = 111, 15.0, 3.0, 120  # SMS N4LO 450 MeV
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120
    
    # SRG \lambda value
    lamb = 1.35
    # lamb = 2.0
    # lamb = 3.0
    # lamb = 6.0

    # Compute and save the momentum distribution
    q_array, q_weights, n_array, n_errors = compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, kmax, kmid, ntot, lamb,
        print_normalization=True, save=True
    )