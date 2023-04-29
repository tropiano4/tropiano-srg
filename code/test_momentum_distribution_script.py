#!/usr/bin/env python3

"""
File: test_momentum_distribution_script.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. Testing how to sum and use batch mode
with vegas.

Last update: April 20, 2023

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
    

class PartialWaveChannel:

    
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
        if (1-(-1)**(L+S)) == 0:
            T = 1
        else:
            T = 0
        
        return T
    
    
class WoodsSaxon:
    """
    Woods-Saxon class. Handles the wave functions associated with the
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
                                                  sp_state.j, sp_state.tau)
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
        u_array = self.sp_wfs[get_orbital_file_name(sp_state.n, sp_state.l,
                                                    sp_state.j, sp_state.tau)]

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
    
    
    def get_cg(self, j1, m1, j2, m2, j3, m3):
        """Evaluates the Clebsch-Gordan coefficient < j1 m1 j2 m2 | j3 m3 >."""
        
        if (m1 + m2 == m3 and np.absolute(m1) <= j1 and np.absolute(m2) <= j2
            and np.absolute(m3) <= j3):
            
            # TESTING
            try:
                return self.cg_table[(j1, m1, j2, m2, j3, m3)]
                
            except KeyError:
                return 0
                
        else:
            return 0
        
        
    def get_sp_wave_functions(
            self, woods_saxon, nucleus_name, Z, N, kmax, kmid, ntot
    ):
        """Set interpolating functions for s.p. wave functions \phi."""
        
        sp_wave_functions = {}
        for sp_state in woods_saxon.occ_states:
            file_name = get_orbital_file_name(sp_state.n, sp_state.l,
                                              sp_state.j, sp_state.tau)
            sp_wave_functions[file_name] = woods_saxon.get_wf_kspace(
                sp_state, kmax, kmid, ntot, interpolate=True)
                    
        return sp_wave_functions
    
    
    def compute_psi(self, n, l, j, m_j, tau, k, theta, phi, sigma):
        """Single-particle wave function including the Clebsch-Gordan
        coefficient and spherical harmonic.
        """
        
        # m_l is determined by m_j and \sigma
        m_l = m_j - sigma
        
        # Check that |m_l| <= l
        if np.abs(m_l) > l:
            # return 0+0j
            return 0

        # Calculate \phi_{n l j \tau}(k) [fm^3/2]
        key = get_orbital_file_name(n, l, j, tau)
        phi_sp_wf = self.sp_wave_functions[key](k)
        
        # Clebsch-Gordan coefficient < l m_l 1/2 \sigma | j m_j >
        cg = self.cg_func(l, m_l, 1/2, sigma, j, m_j)
        
        # Spherical harmonic Y_{l m_l}(\theta, \phi)
        Y_lm = sph_harm(m_l, l, phi, theta)

        # return phi_sp_wf * cg * Y_lm
        psi = phi_sp_wf * cg * Y_lm
        # return psi.real, psi.imag
        return psi.real
    
    
    def vectorized_psi(self):
        pass


    def interpolate_delta_U(
            self, channel, kvnn, kmax, kmid, ntot, generator, lamb,
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
    
    
    def get_delta_U_functions(
            self, kvnn, kmax, kmid, ntot, generator, lamb
    ):
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
    
    
    def get_sp_quantum_numbers(self, x):
        """Given an integration variable x, return a set of four spin
        projections and two s.p. states \alpha and \beta.
        """
        
        index = np.floor(x * self.N_sp).astype(int)
        return self.sp_quantum_numbers[index]
    
    
    def compute_delta_U_matrix_element(
            self, k, theta_k, phi_k, sigma_1, tau_1, sigma_2, tau_2, kp,
            theta_kp, phi_kp, sigma_3, tau_3, sigma_4, tau_4,
            hermitian_conjugate=False
    ):
        """Plane-wave matrix element of \delta U or \delta U^\dagger."""
            
        matrix_element = np.zeros_like(k, dtype='complex')

        # Total spin projections
        M_S = sigma_1 + sigma_2
        M_Sp = sigma_3 + sigma_4
            
        # Total isospin projection
        M_T = tau_1 + tau_2
         
        # Check that total isospin projection is conserved
        if M_T != tau_3 + tau_4:
            # return matrix_element
            return matrix_element.real
            
        # Loop over partial wave channels
        for channel in self.channels:
                
            # Spin CG's
            spin_12_cg = self.cg_func(1/2, sigma_1, 1/2, sigma_2, channel.S,
                                      M_S)
            spin_34_cg = self.cg_func(1/2, sigma_3, 1/2, sigma_4, channel.S,
                                      M_Sp)
                
            # Isospin CG's
            isospin_12_cg = self.cg_func(1/2, tau_1, 1/2, tau_2, channel.T,
                                            M_T)
            isospin_34_cg = self.cg_func(1/2, tau_3, 1/2, tau_4, channel.T,
                                            M_T)
            
            # 1-(-1)^(L+S+T) factors
            lst_factor = 1 - (-1) ** (channel.L + channel.S + channel.T)
            lpst_factor = 1 - (-1) ** (channel.Lp + channel.S + channel.T)
            
            # Partial wave channel matrix element of \delta U
            # ( k J (L S) T | \delta U | k' J (L' S) T )
            channel_str = channel.channel
            if hermitian_conjugate:
                delta_U_partial_wave = (
                        self.delta_U_dagger_functions[channel_str].ev(k, kp)
                )
            else:
                delta_U_partial_wave = (
                    self.delta_U_functions[channel_str].ev(k, kp)
                )
                
            # Loop over M_J from -J to J
            for M_J in np.arange(-channel.J, channel.J + 1):
                
                # Total orbital angular momentum projections are fixed
                M_L = M_J - M_S
                M_Lp = M_J - M_Sp
                
                # Check |M_L| <= L conditions
                if np.abs(M_L) <= channel.L and np.abs(M_Lp) <= channel.Lp:

                    # L-S coupling CG's
                    lsj_cg = self.cg_func(channel.L, M_L, channel.S, M_S,
                                             channel.J, M_J)
                    lpsj_cg = self.cg_func(channel.Lp, M_Lp, channel.S, M_Sp,
                                              channel.J, M_J)
                
                    # Calculate spherical harmonics
                    Y_k = sph_harm(M_L, channel.L, phi_k, theta_k)
                    Y_kp = sph_harm(M_Lp, channel.Lp, phi_kp, theta_kp)
                
                    # Add to matrix element
                    matrix_element += (
                        1/2 * 2/np.pi * spin_12_cg * spin_34_cg * isospin_12_cg
                        * isospin_34_cg * lst_factor * lpst_factor
                        * delta_U_partial_wave * lsj_cg * lpsj_cg * Y_k
                        * np.conj(Y_kp)
                    )
                
        # return matrix_element
        # return matrix_element.real, matrix_element.imag
        return matrix_element.real
    
    
    def vectorized_delta_U_matrix_element():
        pass

    
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
    
    
@vegas.batchintegrand
class DeltaU2Integrand:
    
    
    def __init__(self, tau, q):
        
        # Set instance attributes
        self.tau = tau
        self.q = q
        
    
    def __call__(self):
        pass
    
    
    def get_sp_quantum_numbers(self):
        return None
    
    
class SingleNucleonMomentumDistribution:
    
    
    def __init__(self, nucleus_name, Z, N, tau, channels, kvnn, kmax, kmid,
                 ntot, lamb, j_max=4):
        
        # Call Woods-Saxon class
        self.woods_saxon = WoodsSaxon(nucleus_name, Z, N, run_woodsaxon=False)
    
        # Set Clebsch-Gordan table
        self.cg_table = compute_clebsch_gordan_table(j_max)
        print("Done calculating Clebsch-Gordan table.\n")
        
        # Set instance attributes
        self.channels = channels
        self.kvnn = kvnn
        self.kmax = kmax
        self.kmid = kmid
        self.ntot = ntot
        self.lamb = lamb
        
        
    def compute_I_term(self, tau, q_array):
        """Compute the I * n(q) * I term."""

        I_array = np.zeros_like(q_array)
            
        # Loop over spin projections
        for sigma in np.array([1/2, -1/2]):
                
            # Loop over occupied s.p. states
            for alpha in self.woods_saxon.occ_states:

                # m_l is determined by m_j and \sigma
                m_l = alpha.m_j - sigma

                # Check that \tau_\alpha = \tau and |m_l| <= l
                if alpha.tau == tau and np.abs(m_l) <= alpha.l:
                    
                    # Calculate \phi_{n l j \tau}(k) [fm^3/2]
                    phi_sp_wf = self.woods_saxon.get_wf_kspace(
                        alpha, 10.0, 2.0, 120, interpolate=True
                    )
                    
                    # Clebsch-Gordan coefficient < l m_l 1/2 \sigma | j m_j >
                    cg = self.cg_table[(alpha.l, m_l, 1/2, sigma, alpha.j,
                                        alpha.m_j)]
                    
                    # Spherical harmonic Y_{l m_l}(0, 0)
                    Y_lm = sph_harm(m_l, alpha.l, 0, 0)
                        
                    # \psi(q_vector) choosing z-axis to be along q_vector
                    psi_alpha = phi_sp_wf(q_array) * cg * Y_lm

                    I_array += np.abs(psi_alpha) ** 2
                        
        return I_array
    
    
    def compute_delta_U_term(self, tau, q_array):
        """Compute the sum of the \delta U * n(q) * I term and the
        I * n(q) * \delta U^\dagger term.
        """
        
        delta_U_array = np.zeros_like(q_array)
        delta_U_errors = np.zeros_like(q_array)
        
        # Create a list of all possible combinations of s.p. quantum numbers
        sp_quantum_numbers = []
        for sigma_1 in [1/2, -1/2]:
            for sigma_2 in [1/2, -1/2]:
                for sigma_3 in [1/2, -1/2]:
                    for sigma_4 in [1/2, -1/2]:
                        for alpha in self.woods_saxon.occ_states:
                            
                            n_alpha = alpha.n
                            l_alpha = alpha.l
                            j_alpha = alpha.j
                            m_j_alpha = alpha.m_j
                            tau_alpha = alpha.tau
                            
                            for beta in self.woods_saxon.occ_states:
                                
                                n_beta = beta.n
                                l_beta = beta.l
                                j_beta = beta.j
                                m_j_beta = beta.m_j
                                tau_beta = beta.tau
                                
                                sp_quantum_numbers.append(
                                    (sigma_1, sigma_2, sigma_3, sigma_4,
                                     n_alpha, l_alpha, j_alpha, m_j_alpha,
                                     tau_alpha, n_beta, l_beta, j_beta,
                                     m_j_beta, tau_beta)
                                )
           
        # Integration limits for spin projections and s.p. states
        discrete_limits = [0, 1]
        # Integration limits for relative momenta
        k_limits = [0, 10]
        # Integration limits for C.o.M. momenta
        K_limits = [0, 3]
        # Integration limits for polar angle
        theta_limits = [0, np.pi]
        # Integration limits for azimuthal angle
        phi_limits = [0, 2*np.pi]

        # Set-up integrator with multiple processors
        integ = vegas.Integrator(
            [discrete_limits, k_limits, theta_limits, phi_limits, K_limits,
             theta_limits, phi_limits], nproc=8
        )
            
        # Loop over q
        for i, q in enumerate(q_array):
            
            t0 = time.time()
        
            # integrand = DeltaUIntegrand(tau, q, sp_quantum_numbers, self.md)
            # TESTING
            integrand = DeltaUIntegrand(
                tau, q, nucleus_name, Z, N, sp_quantum_numbers, self.cg_table,
                self.channels, self.kvnn, self.kmax, self.kmid, self.ntot,
                self.lamb
            )

            # Train the integrator
            integ(integrand, nitn=5, neval=1e4)
            # Final result
            result = integ(integrand, nitn=10, neval=1e4)
            
            delta_U_array[i] = result.mean
            delta_U_errors[i] = result.sdev
        
            t1 = time.time()

            percent = (i+1)/len(q_array)*100
            mins = (t1-t0)/60
            print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                        
        return delta_U_array, delta_U_errors
        
        
    def compute_delta_U2_term(self):
        return None
    
    
    def compute_momentum_distribution(
            self, kmax=10.0, kmid=4.0, ntot=100, nmod=70, print_normalization=False,
            ipm_only=False, save=False
    ):
        """Compute the single-nucleon momentum distribution."""

        # Set momentum mesh
        q_array, q_weights = momentum_mesh(kmax, kmid, ntot, nmod)

        # Compute the I term
        I_array = self.compute_I_term(tau, q_array)
        print("Done calculating IPM term.\n")

        # Independent particle model has no \delta U contributions
        if ipm_only:
            
            delta_U_array = np.zeros_like(I_array)
            delta_U_error = np.zeros_like(I_array)
            delta_U2_array = np.zeros_like(I_array)
            delta_U2_error = np.zeros_like(I_array)
        
        # Include \delta U, \delta U^\dagger, and \delta U \delta U^\dagger
        else:
            
            t0 = time.time()

            # Compute \delta U + \delta U^\dagger term
            print("Beginning \delta U linear terms.\n")
            t1 = time.time()
            delta_U_array, delta_U_error = self.compute_delta_U_term(tau,
                                                                     q_array)
            t2 = time.time()
            print(f"Done after {(t2-t1)/3600:.3f} hours.\n")

            # # TESTING
            # delta_U_array = np.zeros_like(I_array)
            # delta_U_error = np.zeros_like(I_array)

            # Compute \delta U \delta U^\dagger term using vegas
            print("Beginning \delta U \delta U^\dagger term.\n")
            t3 = time.time()
            # ...
            t4 = time.time()
            print(f"Done after {(t4-t3)/3600:.3f} hours.\n")

            # TESTING
            delta_U2_array = np.zeros_like(I_array)
            delta_U2_error = np.zeros_like(I_array)

            t5 = time.time()
            print(f"Total time elapsed: {(t5-t0)/3600:.3f} hours.\n")
                
        # Combine each term for the total momentum distribution [fm^3]
        n_array = I_array + delta_U_array + delta_U2_array
        n_error = np.sqrt(delta_U_error ** 2 + delta_U2_error ** 2)

        if print_normalization:
            normalization = self.compute_normalization(q_array, q_weights,
                                                       n_array)
            print(f"Normalization = {normalization:.5f}.")
            
        if save and not(ipm_only):  # Do not save IPM-only data
            self.save_momentum_distribution(
                q_array, q_weights, n_array, n_error, I_array, delta_U_array,
                delta_U_error, delta_U2_array, delta_U2_error
            )
        
        return q_array, q_weights, n_array, n_error
    
    
    def compute_normalization(self, q_array, q_weights, n_array):
        """Compute the normalization of the momentum distribution."""

        return 4*np.pi * np.sum(q_weights * q_array**2 * n_array)
    
    
    def get_file_name(self, nucleus_name, nucleon, lamb):
        """File name of the single-nucleon momentum distribution."""
        
        file_name = replace_periods(
            f"{nucleus_name}_{nucleon}_momentum_distribution_{lamb}"
        )
        
        return file_name
    
        
    def save_momentum_distribution(
            self, q_array, q_weights, n_array, n_error, I_array, delta_U_array,
            delta_U_error, delta_U2_array, delta_U2_error
    ):
        """Save the momentum distribution along with the isolated contributions."""
        
        data = np.vstack(
            (q_array, q_weights, n_array, n_error, I_array, delta_U_array,
             delta_U_error, delta_U2_array, delta_U2_error)
        ).T
       
        hdr = (
            "q, q weight, n(q), n(q) error, I, \delta U + \delta U^\dagger,"
            " \delta U + \delta U^\dagger error, \delta U^2, \delta U^2 error\n"
        )
        
        if tau == 1/2:
            nucleon = 'proton'
        elif tau == -1/2:
            nucleon = 'neutron'

        file_name = self.get_file_name(nucleus_name, nucleon, lamb)
        
        np.savetxt(file_name + '.txt', data, header=hdr)

        
    def load_momentum_distribution(self, nucleus_name, nucleon, lamb):
        """Load and return the momentum distribution along with the isolated
        contributions.
        """
        
        file_name = self.get_file_name(nucleus_name, nucleon, lamb)
        
        data = np.loadtxt(file_name + '.txt')
        
        q_array = data[:, 0]
        q_weights = data[:, 1]
        n_array = data[:, 2]
        n_error = data[:, 3]
        I_array = data[:, 4]
        delta_U_array = data[:, 5]
        delta_U_error = data[:, 6]
        delta_U2_array = data[:, 7]
        delta_U2_error = data[:, 8]
        
        return (q_array, q_weights, n_array, n_error, I_array, delta_U_array,
                delta_U_error, delta_U2_array, delta_U2_error)
    
        
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
    
    # Call single-nucleon momentum distribution class
    snmd = SingleNucleonMomentumDistribution(nucleus_name, Z, N, tau, channels,
                                             kvnn, kmax, kmid, ntot, lamb)

    # Compute and save the momentum distribution
    q_array, q_weights, n_array, n_error = snmd.compute_momentum_distribution(
        print_normalization=True, save=True
    )