#!/usr/bin/env python3

"""
File: test_momentum_distribution_script.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: February 7, 2023

This script serves as a testbed for calculating momentum distributions using
mean field approximations for initial and final states and applying SRG
transformations to the operator. This differs from the previous momentum
distribution calculations by directly utilizing single-particle wave functions
instead of a local density approximation. Testing how to do batch mode with
vegas.

Last update: April 17, 2023

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
    
    
class SingleParticleBasis:
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
                file_name = get_orbital_file_name(sp_state)
                if run_woodsaxon:
                    shutil.move(file_name,
                                self.woods_saxon_directory + file_name)
                data = np.loadtxt(self.woods_saxon_directory + file_name)
                # Use file name as the key
                self.sp_wfs[file_name] = data[:, 1]

            
    def run_woods_saxon_code(self, nucleus_name, Z, N, n_max, l_max, rmax, ntab):
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

    k = la.norm(k_vector, axis=0)
    theta = np.arccos(k_vector[2]/k)
    phi = np.arctan2(k_vector[1], k_vector[0])

    return k, theta, phi


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


def get_sp_wave_functions(sp_basis, kmax, kmid, ntot):
    """Set interpolating functions for s.p. wave functions \phi."""
    
    occ_states = sp_basis.occ_states

    phi_functions = {}
    for sp_state in occ_states: 
        file_name = get_orbital_file_name(sp_state)
        phi_functions[file_name] = sp_basis.get_wf_kspace(
            sp_state, kmax, kmid, ntot, interpolate=True)
            
    return phi_functions


def psi(sp_state, q_vector, sigma, cg_table, phi_functions):
    """Single-particle wave function including the Clebsch-Gordan coefficient 
    and spherical harmonic.
    """
        
    # Unpack q_vector into magnitude and angles
    q, theta, phi = get_vector_components(q_vector)
        
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
    

def compute_I_term(q_array, tau, occ_states, cg_table, phi_functions):
    """Compute the I * n(q) * I term."""
        
    I_array = np.zeros_like(q_array)
        
    # Loop over spin projections
    for sigma in np.array([1/2, -1/2]):
            
        # Loop over occupied s.p. states
        for alpha in occ_states:
                
            # m_l is determined by m_j and \sigma
            m_l_alpha = alpha.m_j - sigma
                
            # Check that \tau_\alpha = \tau and |m_l| <= l:
            if alpha.tau == tau and abs(m_l_alpha) <= alpha.l:
            # # TESTING
            # # if alpha.tau == tau and abs(m_l_alpha) <= alpha.l and alpha.l == 0:
            # if alpha.tau == tau and abs(m_l_alpha) <= alpha.l and alpha.l == 1 and alpha.j == 3/2:
                    
                # Loop over q
                psi_alpha_array = np.zeros_like(q_array, dtype='complex')
                for i, q in enumerate(q_array):
                        
                    # Choose z-axis to be along q_vector
                    q_vector = np.array([0, 0, q])
                    
                    # Single-particle wave function
                    psi_alpha_array[i] = psi(alpha, q_vector, sigma, cg_table,
                                             phi_functions)
                    
                I_array += abs(psi_alpha_array)**2
                    
    return I_array


def delta_U_quantum_numbers(tau, occ_states, cg_table, channels):
  """Returns a list of every combination of quantum numbers in the \delta U
  linear terms that gives a non-zero product of Clebsch-Gordan coefficients.
  """
    
  spins = np.array([1/2, -1/2])
  quantum_number_combinations = []
    
  # Many loops over quantum numbers
  for channel in channels:
        
    L, Lp, S, J = get_channel_quantum_numbers(channel)
    T = get_total_isospin(L, S)
    
    # L S T factor
    lst_factor = 1-(-1)**(L+S+T)
    # L' S T factor
    lpst_factor = 1-(-1)**(Lp+S+T)
        
    # Channel isospin projection M_T
    for M_T in np.arange(-T, T+1):
      
      # Occupied single-particle state \alpha
      for alpha in occ_states:
                
        # Kronecker \delta function \delta_{\tau_\alpha \tau}
        del_alpha = int(alpha.tau == tau)
        
        # TESTING
        # del_swave_alpha = int(alpha.l == 0)
        # del_swave_alpha = int(alpha.l == 1 and alpha.j == 3/2)

        # Occupied single-particle state \beta
        for beta in occ_states:
                    
          # Kronecker \delta function \delta_{\tau_\beta \tau}
          del_beta = int(beta.tau == tau)
          
          # TESTING
          # del_swave_beta = int(beta.l == 0)
          # del_swave_beta = int(beta.l == 1 and beta.j == 3/2)
                    
          # \tau_\alpha \tau_\beta T M_T CG
          if alpha.tau + beta.tau == M_T:
            isospin_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, T, M_T)]
          else:
            isospin_cg = 0
          
          # Channel angular momentum projection M_J
          for M_J in np.arange(-J, J+1):

            # Channel spin projection M_S
            for M_S in np.arange(-S, S+1):
                 
              # Channel orbital angular momentum projection M_L  
              M_L = M_J - M_S
                    
              # L S J CG
              if abs(M_L) <= L:
                lsj_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
              else:
                lsj_cg = 0
              
              # Channel spin projection M_S'   
              for M_Sp in np.arange(-S, S+1):
                
                # Channel orbital angular momentum projection M_L'  
                M_Lp = M_J - M_Sp
                        
                # L' S J CG
                if abs(M_Lp) <= Lp:
                  lpsj_cg = cg_table[(Lp, M_Lp, S, M_Sp, J, M_J)]
                else:
                  lpsj_cg = 0
                
                # Single-particle spin projection \sigma_1   
                for sigma_1 in spins:
                                    
                  # \sigma_1 l_\alpha CG
                  if abs(alpha.m_j-sigma_1) <= alpha.l:
                    sig1_alpha_cg = cg_table[(alpha.l, alpha.m_j-sigma_1, 1/2,
                                              sigma_1, alpha.j, alpha.m_j)]
                  else:
                    sig1_alpha_cg = 0
                   
                  # Single-particle spin projection \sigma_2   
                  for sigma_2 in spins:
                                        
                    # \sigma_2 l_\beta CG
                    if abs(beta.m_j-sigma_2) <= beta.l:
                      sig2_beta_cg = cg_table[(beta.l, beta.m_j-sigma_2, 1/2,
                                               sigma_2, beta.j, beta.m_j)]
                    else:
                      sig2_beta_cg = 0
                                        
                    # \sigma_1 \sigma_2 S M_S CG
                    if sigma_1 + sigma_2 == M_S:
                      spin_12_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S,
                                             M_S)]
                    else:
                      spin_12_cg = 0
                                            
                    # \sigma_1 \sigma_2 S M_S' CG
                    if sigma_1 + sigma_2 == M_Sp:
                      spin_12p_cg = cg_table[(1/2, sigma_1, 1/2, sigma_2, S,
                                              M_Sp)]
                    else:
                      spin_12p_cg = 0
                    
                    # Single-particle spin projection \sigma   
                    for sigma in spins:
                                            
                      # \sigma l_\alpha CG
                      if abs(alpha.m_j-sigma) <= alpha.l:
                        sig_alpha_cg = cg_table[(
                          alpha.l, alpha.m_j-sigma, 1/2, sigma, alpha.j,
                          alpha.m_j
                        )]
                      else:
                        sig_alpha_cg = 0
                                        
                      # \sigma l_\beta CG
                      if abs(beta.m_j-sigma) <= beta.l:
                        sig_beta_cg = cg_table[(beta.l, beta.m_j-sigma, 1/2,
                                                sigma, beta.j, beta.m_j)]
                      else:
                        sig_beta_cg = 0
                      
                      # Single-particle spin projection \sigma'    
                      for sigmap in spins:
                                                
                        # \sigma' l_\alpha CG
                        if abs(alpha.m_j-sigmap) <= alpha.l:
                          sigp_alpha_cg = cg_table[(
                            alpha.l, alpha.m_j-sigmap, 1/2, sigmap, alpha.j,
                            alpha.m_j
                          )]
                        else:
                          sigp_alpha_cg = 0
                                            
                        # \sigma' l_\beta CG
                        if abs(beta.m_j-sigmap) <= beta.l:
                          sigp_beta_cg = cg_table[(
                            beta.l, beta.m_j-sigmap, 1/2, sigmap, beta.j,
                            beta.m_j
                          )]
                        else:
                          sigp_beta_cg = 0
                                                    
                        # \sigma \sigma' S M_S CG
                        if sigma + sigmap == M_S:
                          spin_ssp_cg = cg_table[(1/2, sigma, 1/2, sigmap, S,
                                                  M_S)]
                        else:
                          spin_ssp_cg = 0
                                                    
                        # \sigma \sigma' S M_S' CG
                        if sigma + sigmap == M_Sp:
                          spin_sspp_cg = cg_table[(1/2, sigma, 1/2, sigmap, S,
                                                   M_Sp)]
                        else:
                          spin_sspp_cg = 0
                        
                        # Product of CG's and factors    
                        product = (
                          isospin_cg ** 2 * lsj_cg * lpsj_cg * lst_factor
                          * lpst_factor * sig1_alpha_cg * sig2_beta_cg * (
                            del_alpha * sig_alpha_cg * sigp_beta_cg
                            - (-1) ** (T-1) * del_beta * sigp_alpha_cg
                            * sig_beta_cg
                          ) * (
                            spin_12_cg * spin_sspp_cg + spin_12p_cg
                            * spin_ssp_cg
                          )
                        ) # * del_swave_alpha * del_swave_beta
                        # TESTING
                        
                        # Add this set of quantum numbers to the list if the
                        # product is nonzero
                        if product != 0:
                                                                        
                          d = {
                            'channel': channel,
                            'T': T,
                            'M_T': M_T,
                            'alpha': alpha,
                            'beta': beta,
                            'J': J,
                            'M_J': M_J,
                            'S': S, 
                            'M_S': M_S,
                            'M_Sp': M_Sp,
                            'L': L,
                            'M_L': M_L,
                            'Lp': Lp,
                            'M_Lp': M_Lp,
                            'sigma_1': sigma_1,
                            'sigma_2': sigma_2,
                            'sigma': sigma,
                            'sigmap': sigmap
                          }
                                                                        
                          quantum_number_combinations.append(d)
            
  return quantum_number_combinations


@vegas.batchintegrand
class delta_U_term_integrand:
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
            
        # Calculate the Jacobian determinant
        jacobian = k**2 * np.sin(theta_k) * K**2 * np.sin(theta_K)
        
        integrand = np.zeros_like(jacobian, dtype='complex')
        for d in self.quantum_numbers:
            
            # Unpack dictionary
            sigma_1 = d['sigma_1']
            sigma_2 = d['sigma_2']
            sigma = d['sigma']
            sigmap = d['sigmap']
            alpha = d['alpha']
            beta = d['beta']
            channel = d['channel']
            T = d['T']
            M_T = d['M_T']
            S = d['S']
            M_S = d['M_S']
            M_Sp = d['M_Sp']
            L = d['L']
            M_L = d['M_L']
            Lp = d['Lp']
            M_Lp = d['M_Lp']
            J = d['J']
            M_J = d['M_J']
            
            # \sigma_1 \sigma_2 S M_S CG
            if sigma_1 + sigma_2 == M_S:
                spin_12_cg = self.cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
            else:
                spin_12_cg = 0
            # \sigma \sigma' S M_S CG
            if sigma + sigmap == M_S:
                spin_ssp_cg = self.cg_table[(1/2, sigma, 1/2, sigmap, S, M_S)]
            else:
                spin_ssp_cg = 0
            # \sigma_1 \sigma_2 S M_S' CG
            if sigma_1 + sigma_2 == M_Sp:
                spin_12p_cg = self.cg_table[(1/2, sigma_1, 1/2, sigma_2, S,
                                             M_Sp)]
            else:
                spin_12p_cg = 0
            # \sigma \sigma' S M_S' CG
            if sigma + sigmap == M_Sp:
                spin_sspp_cg = self.cg_table[(1/2, sigma, 1/2, sigmap, S, M_Sp)]
            else:
                spin_sspp_cg = 0
                
            # \tau_\alpha \tau_\beta T M_T CG
            isospin_cg = self.cg_table[(1/2, alpha.tau, 1/2, beta.tau, T, M_T)]
            
            # L S J CG
            lsj_cg = self.cg_table[(L, M_L, S, M_S, J, M_J)]
            # L' S J CG
            lpsj_cg = self.cg_table[(Lp, M_Lp, S, M_Sp, J, M_J)]
            
            # L S T factor
            lst_factor = 1-(-1)**(L+S+T)
            # L' S T factor
            lpst_factor = 1-(-1)**(Lp+S+T)
            
            # Spherical harmonics
            Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
            Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
            Y_L_qK = sph_harm(M_L, L, phi_qK, theta_qK)
            Y_Lp_k = sph_harm(M_Lp, Lp, phi_k, theta_k)
            
            # \delta U_{L,L'}(k, q-K/2)
            delta_U_partial_wave = self.delta_U_functions[channel].ev(k, qK)
            # \delta U^\dagger_{L,L'}(q-K/2, k)
            delta_U_dag_partial_wave = (
                self.delta_U_dagger_functions[channel].ev(qK, k)
            )
            
            # Single-particle wave functions
            psi_alpha_k1 = psi(alpha, K_vector/2+k_vector, sigma_1,
                               self.cg_table, self.phi_functions)
            psi_beta_k2 = psi(beta, K_vector/2-k_vector, sigma_2, self.cg_table,
                              self.phi_functions)
            psi_alpha_q = psi(alpha, q_vector, sigma, self.cg_table,
                              self.phi_functions)
            psi_beta_Kq = psi(beta, K_vector-q_vector, sigmap, self.cg_table,
                              self.phi_functions)
            psi_alpha_Kq = psi(alpha, K_vector-q_vector, sigmap, self.cg_table,
                               self.phi_functions)
            psi_beta_q = psi(beta, q_vector, sigma, self.cg_table,
                             self.phi_functions)
            
            # Kronecker \delta functions
            del_alpha = int(alpha.tau == self.tau)
            del_beta = int(beta.tau == self.tau)
            
            integrand += 1/2 * (1/2 * 2/np.pi) * jacobian * (
                isospin_cg ** 2 * lsj_cg * lpsj_cg * lst_factor * lpst_factor
                * (
                    np.conj(psi_alpha_k1) * np.conj(psi_beta_k2)
                    * spin_12_cg * spin_sspp_cg * Y_L_k * np.conj(Y_Lp_qK)
                    * delta_U_partial_wave * (
                        del_alpha * psi_beta_Kq * psi_alpha_q 
                        - (-1) ** (T-1) * del_beta * psi_alpha_Kq * psi_beta_q
                    )
                    + psi_alpha_k1 * psi_beta_k2 * spin_ssp_cg * spin_12p_cg
                    * Y_L_qK * np.conj(Y_Lp_k) * delta_U_dag_partial_wave * (
                        del_alpha * np.conj(psi_beta_Kq) * np.conj(psi_alpha_q)
                        - (-1) ** (T-1) * del_beta * np.conj(psi_alpha_Kq)
                        * np.conj(psi_beta_q)
                    )
                )  
            )

        return integrand.real
    

def compute_delta_U_term(
        q_array, tau, occ_states, cg_table, channels, phi_functions,
        delta_U_functions, delta_U_dagger_functions
):
    """Compute the sum of the \delta U * n(q) * I term and the 
    I * n(q) * \delta U^\dagger term.
    """
        
    delta_U_array = np.zeros_like(q_array)
    delta_U_errors = np.zeros_like(q_array)
    
    # Get an organized list of quantum numbers to sum over
    quantum_numbers = delta_U_quantum_numbers(tau, occ_states, cg_table,
                                              channels)
        
    # Relative momenta from 0 to maximum of momentum mesh
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    # Polar angle
    theta_limits = [0, np.pi]
    # Azimuthal angle
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = delta_U_term_integrand(
            q, tau, quantum_numbers, cg_table, phi_functions, delta_U_functions,
            delta_U_dagger_functions
        )

        # # Train the integrator
        # integ(integrand, nitn=5, neval=200)
        # # Final result
        # result = integ(integrand, nitn=10, neval=1e3)
        
        # Train the integrator
        integ(integrand, nitn=5, neval=1e3)
        # Final result
        result = integ(integrand, nitn=10, neval=1e3)

        delta_U_array[i] = result.mean
        delta_U_errors[i] = result.sdev

        t1 = time.time()

        percent = (i+1)/len(q_array)*100
        mins = (t1-t0)/60
        print(f"{percent:.2f}% done after {mins:.3f} minutes.")
                    
    return delta_U_array, delta_U_errors


def delta_U2_quantum_numbers(tau, occ_states, cg_table, channels):
  """Returns a list of every combination of quantum numbers in the \delta U^2
  term that gives a non-zero product of Clebsch-Gordan coefficients.
  """
    
  spins = np.array([1/2, -1/2])
  quantum_number_combinations = []
    
  # Many loops over quantum numbers
  for channel_1 in channels:
  
    L, Lp, S, J = get_channel_quantum_numbers(channel_1)
    T = get_total_isospin(L, S)
        
    # L S T factor
    lst_factor = 1-(-1)**(L+S+T)
    # L' S T factor
    lpst_factor = 1-(-1)**(Lp+S+T)
        
    for channel_2 in channels:
            
      Lpp, Lppp, Sp, Jp = get_channel_quantum_numbers(channel_2)
      Tp = get_total_isospin(Lpp, Sp)
            
      # Evaluate Kronecker \delta in S and S'
      if Sp == S:
            
        # L'' S T' factor
        lppstp_factor = 1-(-1)**(Lpp+S+Tp)
        # L''' S T' factor
        lpppstp_factor = 1-(-1)**(Lppp+S+Tp)
        
        # Channel isospin projection M_T
        for M_T in np.arange(-T, T+1):
          # Channel isospin projection M_T'
          for M_Tp in np.arange(-Tp, Tp+1):
            # Occupied single-particle state \alpha
            for alpha in occ_states:
                
              # TESTING
              # del_swave_alpha = int(alpha.l == 0)
              # del_swave_alpha = int(alpha.l == 1 and alpha.j == 3/2)
                
              # Occupied single-particle state \beta
              for beta in occ_states:
                  
                # TESTING
                # del_swave_beta = int(beta.l == 0)
                # del_swave_beta = int(beta.l == 1 and beta.j == 3/2)
                                
                # \tau_\alpha \tau_\beta T M_T CG
                if alpha.tau + beta.tau == M_T:
                  isospin_ab_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, T,
                                            M_T)]
                else:
                  isospin_ab_cg = 0
                                    
                # \tau_\alpha \tau_\beta T' M_T' CG
                if alpha.tau + beta.tau == M_Tp:
                  isospin_abp_cg = cg_table[(1/2, alpha.tau, 1/2, beta.tau, Tp,
                                             M_Tp)]
                else:
                  isospin_abp_cg = 0
                
                # Single-particle isospin projection \tau'    
                for taup in spins:
                                    
                  # \tau \tau' T M_T CG
                  if tau + taup == M_T:
                    isospin_ttp_cg = cg_table[(1/2, tau, 1/2, taup, T, M_T)]
                  else:
                    isospin_ttp_cg = 0
                                        
                  # \tau \tau' T' M_T' CG
                  if tau + taup == M_Tp:
                    isospin_ttpp_cg = cg_table[(1/2, tau, 1/2, taup, Tp, M_Tp)]
                  else:
                    isospin_ttpp_cg = 0
                  
                  # Channel angular momentum projection M_J    
                  for M_J in np.arange(-J, J+1):
                    # Channel angular momentum projection M_J'
                    for M_Jp in np.arange(-Jp, Jp+1):
                      # Channel spin projection M_S
                      for M_S in np.arange(-S, S+1):
                        
                        # Channel orbital angular momentum projection M_L  
                        M_L = M_J - M_S
                                                
                        # L S J CG
                        if abs(M_L) <= L:
                          lsj_cg = cg_table[(L, M_L, S, M_S, J, M_J)]
                        else:
                          lsj_cg = 0
                        
                        # Channel spin projection M_S'    
                        for M_Sp in np.arange(-S, S+1):
                          
                          # Channel orbital angular momentum projection M_L'  
                          M_Lp = M_J - M_Sp
                                                
                          # L' S J CG
                          if abs(M_Lp) <= Lp:
                            lpsj_cg = cg_table[(Lp, M_Lp, S, M_Sp, J, M_J)]
                          else:
                            lpsj_cg = 0
                                                    
                          # Evaluate Kronecker \delta in M_S' and M_S''
                          M_Spp = M_Sp
                          
                          # Channel orbital angular momentum projection M_L''
                          M_Lpp = M_Jp - M_Spp
                                                    
                          # L'' S' J' CG
                          if abs(M_Lpp) <= Lpp and abs(M_Spp) <= Sp:
                            lppsj_cg = cg_table[(Lpp, M_Lpp, Sp, M_Spp, Jp,
                                                 M_Jp)]
                          else:
                            lppsj_cg = 0
                          
                          # Channel spin projection M_S'''    
                          for M_Sppp in np.arange(-Sp, Sp+1):
                            
                            # Channel orbital angular momentum projection M_L'''  
                            M_Lppp = M_Jp - M_Sppp
                                                        
                            # L''' S' J' CG
                            if abs(M_Lppp) <= Lppp and abs(M_Sppp) <= Sp:
                              lpppsj_cg = cg_table[(Lppp, M_Lppp, Sp, M_Sppp,
                                                    Jp, M_Jp)]
                            else:
                              lpppsj_cg = 0
                            
                            # Single-particle spin projection \sigma_1    
                            for sigma_1 in spins:
                                                            
                              # \sigma_1 l_\alpha CG
                              if abs(alpha.m_j-sigma_1) <= alpha.l:
                                sig1_alpha_cg = cg_table[(
                                  alpha.l, alpha.m_j-sigma_1, 1/2, sigma_1,
                                  alpha.j, alpha.m_j
                                )]
                              else:
                                sig1_alpha_cg = 0
                                
                              # Single-particle spin projection \sigma_2                           
                              for sigma_2 in spins:
                                                                
                                # \sigma_2 l_\beta CG
                                if abs(beta.m_j-sigma_2) <= beta.l:
                                  sig2_beta_cg = cg_table[(
                                    beta.l, beta.m_j-sigma_2, 1/2, sigma_2,
                                    beta.j, beta.m_j
                                  )]
                                else:
                                  sig2_beta_cg = 0
            
                                # \sigma_1 \sigma_2 S M_S CG
                                if sigma_1 + sigma_2 == M_S:
                                  spin_12_cg = cg_table[(1/2, sigma_1, 1/2,
                                                         sigma_2, S, M_S)]
                                else:
                                  spin_12_cg = 0
                                
                                # Single-particle spin projection \sigma_3  
                                for sigma_3 in spins:
                                                                    
                                  # \sigma_3 l_\alpha CG
                                  if abs(alpha.m_j-sigma_3) <= alpha.l:
                                    sig3_alpha_cg = cg_table[(
                                      alpha.l, alpha.m_j-sigma_3, 1/2, sigma_3,
                                      alpha.j, alpha.m_j
                                    )]
                                  else:
                                    sig3_alpha_cg = 0
                                                                        
                                  # \sigma_3 l_\beta CG
                                  if abs(beta.m_j-sigma_3) <= beta.l:
                                    sig3_beta_cg = cg_table[(
                                      beta.l, beta.m_j-sigma_3, 1/2, sigma_3,
                                      beta.j, beta.m_j
                                    )]
                                  else:
                                    sig3_beta_cg = 0
                                  
                                  # Single-particle spin projection \sigma_4    
                                  for sigma_4 in spins:
                                                                        
                                    # \sigma_4 l_\alpha CG
                                    if abs(alpha.m_j-sigma_4) <= alpha.l:
                                      sig4_alpha_cg = cg_table[(
                                        alpha.l, alpha.m_j-sigma_4, 1/2,
                                        sigma_4, alpha.j, alpha.m_j
                                      )]
                                    else:
                                      sig4_alpha_cg = 0
                                                                        
                                    # \sigma_4 l_\beta CG
                                    if abs(beta.m_j-sigma_4) <= beta.l:
                                      sig4_beta_cg = cg_table[(
                                        beta.l, beta.m_j-sigma_4, 1/2, sigma_4,
                                        beta.j, beta.m_j
                                      )]
                                    else:
                                      sig4_beta_cg = 0
                                                                        
                                    # \sigma_3 \sigma_4 S' M_S''' CG
                                    if sigma_3 + sigma_4 == M_Sppp:
                                      spin_34_cg = cg_table[(
                                        1/2, sigma_3, 1/2, sigma_4, Sp, M_Sppp
                                      )]
                                    else:
                                      spin_34_cg = 0
                                    
                                    # Product of CG's and factors  
                                    product = (
                                      sig1_alpha_cg * sig2_beta_cg * spin_12_cg
                                      * spin_34_cg * isospin_ab_cg
                                      * isospin_ttp_cg * isospin_ttpp_cg
                                      * isospin_abp_cg * lsj_cg * lpsj_cg
                                      * lppsj_cg * lpppsj_cg * lst_factor
                                      * lpst_factor * lppstp_factor
                                      * lpppstp_factor * (
                                        sig3_alpha_cg * sig4_beta_cg
                                        - (-1) ** (Tp-1) * sig3_beta_cg
                                        * sig4_alpha_cg
                                      )
                                    ) # * del_swave_alpha * del_swave_beta
                                    # TESTING
                                    
                                    # Add this set of quantum numbers to the
                                    # list if the product is nonzero
                                    if product != 0:
                                                                            
                                      d = {
                                        'channel_1': channel_1,
                                        'channel_2': channel_2,
                                        'T': T,
                                        'M_T': M_T,
                                        'Tp': Tp,
                                        'M_Tp': M_Tp,
                                        'alpha': alpha,
                                        'beta': beta,
                                        'J': J,
                                        'M_J': M_J,
                                        'Jp': Jp,
                                        'M_Jp': M_Jp,
                                        'S': S,
                                        'M_S': M_S,
                                        'M_Sp': M_Sp,
                                        'Sp': Sp,
                                        'M_Spp': M_Spp,
                                        'M_Sppp': M_Sppp,
                                        'L': L,
                                        'M_L': M_L,
                                        'Lp': Lp,
                                        'M_Lp': M_Lp,
                                        'Lpp': Lpp,
                                        'M_Lpp': M_Lpp,
                                        'Lppp': Lppp,
                                        'M_Lppp': M_Lppp,
                                        'sigma_1': sigma_1,
                                        'sigma_2': sigma_2,
                                        'sigma_3': sigma_3,
                                        'sigma_4': sigma_4,
                                        'taup': taup
                                      }
                                                                                                
                                      quantum_number_combinations.append(d)
                                                                        
  return quantum_number_combinations


@vegas.batchintegrand
class delta_U2_term_integrand:
    """Evaluate the integrand of the \delta U \delta U^\dagger term."""
    
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
        
        # Relative momenta k
        k = momenta_array[:,0]
        theta_k = momenta_array[:,1]
        phi_k = momenta_array[:,2]
        k_vector = build_vector(k, theta_k, phi_k)
            
        # Relative momenta k'
        kp = momenta_array[:,3]
        theta_kp = momenta_array[:,4]
        phi_kp = momenta_array[:,5]
        kp_vector = build_vector(kp, theta_kp, phi_kp)
            
        # C.o.M. momenta K
        K = momenta_array[:,6]
        theta_K = momenta_array[:,7]
        phi_K = momenta_array[:,8]
        K_vector = build_vector(K, theta_K, phi_K)
        
        # Choose z-axis to be along q_vector
        q_vector = np.zeros_like(k_vector)
        q_vector[-1, :] = self.q
        
        # Calculate vector q-K/2
        qK_vector = q_vector - K_vector/2
        qK, theta_qK, phi_qK = get_vector_components(qK_vector)
            
        # Calculate the Jacobian determinant
        jacobian = (k**2 * np.sin(theta_k) * kp**2 * np.sin(theta_kp) * K**2
                    * np.sin(theta_K))
        
        integrand = np.zeros_like(jacobian, dtype='complex')
        for d in self.quantum_numbers:
            
            # Unpack dictionary
            channel_1 = d['channel_1']
            channel_2 = d['channel_2']
            T = d['T']
            M_T = d['M_T']
            Tp = d['Tp']
            M_Tp = d['M_Tp']
            alpha = d['alpha']
            beta = d['beta']
            J = d['J']
            M_J = d['M_J']
            Jp = d['Jp']
            M_Jp = d['M_Jp']
            S = d['S']
            M_S = d['M_S']
            M_Sp = d['M_Sp']
            Sp = d['Sp']
            M_Spp = d['M_Spp']
            M_Sppp = d['M_Sppp']
            L = d['L']
            M_L = d['M_L']
            Lp = d['Lp']
            M_Lp = d['M_Lp']
            Lpp = d['Lpp']
            M_Lpp = d['M_Lpp']
            Lppp = d['Lppp']
            M_Lppp = d['M_Lppp']
            sigma_1 = d['sigma_1']
            sigma_2 = d['sigma_2']
            sigma_3 = d['sigma_3']
            sigma_4 = d['sigma_4']
            taup = d['taup']
            
            # \sigma_1 \sigma_2 S M_S CG
            spin_12_cg = self.cg_table[(1/2, sigma_1, 1/2, sigma_2, S, M_S)]
            # \sigma_3 \sigma_4 S' M_S''' CG
            spin_34_cg = self.cg_table[(1/2, sigma_3, 1/2, sigma_4, Sp, M_Sppp)]
                
            # \tau_\alpha \tau_\beta T M_T CG
            isospin_ab_cg = self.cg_table[(1/2, alpha.tau, 1/2, beta.tau, T,
                                           M_T)]
            # \tau \tau' T M_T CG
            isospin_ttp_cg = self.cg_table[(1/2, self.tau, 1/2, taup, T, M_T)]
            # \tau \tau' T' M_T' CG
            isospin_ttpp_cg = self.cg_table[(1/2, self.tau, 1/2, taup, Tp,
                                             M_Tp)]
            # \tau_\alpha \tau_\beta T' M_T' CG
            isospin_abp_cg = self.cg_table[(1/2, alpha.tau, 1/2, beta.tau, Tp,
                                            M_Tp)]

            # L S J CG
            lsj_cg = self.cg_table[(L, M_L, S, M_S, J, M_J)]
            # L' S J CG
            lpsj_cg = self.cg_table[(Lp, M_Lp, S, M_Sp, J, M_J)]
            # L'' S J' CG
            lppsj_cg = self.cg_table[(Lpp, M_Lpp, Sp, M_Spp, Jp, M_Jp)]
            # L''' S J' CG
            lpppsj_cg = self.cg_table[(Lppp, M_Lppp, Sp, M_Sppp, Jp, M_Jp)]
            
            # L S T factor
            lst_factor = 1-(-1)**(L+S+T)
            # L' S T factor
            lpst_factor = 1-(-1)**(Lp+S+T)
            # L'' S T' factor
            lppst_factor = 1-(-1)**(Lpp+Sp+Tp)
            # L''' S T' factor
            lpppst_factor = 1-(-1)**(Lppp+Sp+Tp)
            
            # Spherical harmonics
            Y_L_k = sph_harm(M_L, L, phi_k, theta_k)
            Y_Lp_qK = sph_harm(M_Lp, Lp, phi_qK, theta_qK)
            Y_Lpp_qK = sph_harm(M_Lpp, Lpp, phi_qK, theta_qK)
            Y_Lppp_kp = sph_harm(M_Lppp, Lppp, phi_kp, theta_kp)
            
            # \delta U_{L,L'}(k, q-K/2)
            delta_U_partial_wave = self.delta_U_functions[channel_1].ev(k, qK)
            # \delta U^\dagger_{L'',L'''}(q-K/2, k')
            delta_U_dag_partial_wave = (
                self.delta_U_dagger_functions[channel_2].ev(qK, kp)
            )
            
            # Single-particle wave functions
            psi_alpha_k1 = psi(alpha, K_vector/2+k_vector, sigma_1,
                               self.cg_table, self.phi_functions)
            psi_beta_k2 = psi(beta, K_vector/2-k_vector, sigma_2,
                              self.cg_table, self.phi_functions)
            psi_alpha_k3 = psi(alpha, K_vector/2+kp_vector, sigma_3,
                               self.cg_table, self.phi_functions)
            psi_beta_k4 = psi(beta, K_vector/2-kp_vector, sigma_4,
                              self.cg_table, self.phi_functions)
            psi_alpha_k4 = psi(alpha, K_vector/2-kp_vector, sigma_4,
                               self.cg_table, self.phi_functions)
            psi_beta_k3 = psi(beta, K_vector/2+kp_vector, sigma_3,
                              self.cg_table, self.phi_functions)
        
            integrand += 1/4 * (1/2 * 2/np.pi) ** 2 * jacobian * (
                spin_12_cg * spin_34_cg * isospin_ab_cg * isospin_ttp_cg 
                * isospin_ttpp_cg * isospin_abp_cg * lsj_cg * lpsj_cg * lppsj_cg
                * lpppsj_cg * lst_factor * lpst_factor * lppst_factor
                * lpppst_factor * Y_L_k * np.conj(Y_Lp_qK) * Y_Lpp_qK
                * np.conj(Y_Lppp_kp) * delta_U_partial_wave
                * delta_U_dag_partial_wave * np.conj(psi_alpha_k1)
                * np.conj(psi_beta_k2) * (
                    psi_alpha_k3 * psi_beta_k4
                    - (-1) ** (Tp-1) * psi_alpha_k4 * psi_beta_k3
                )
            )
        
        return integrand.real


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
        
    # Relative momenta from 0 to maximum of momentum mesh
    k_limits = [0, 10]
    # C.o.M. momenta up to 3 fm^-1
    K_limits = [0, 3]
    # Polar angle
    theta_limits = [0, np.pi]
    # Azimuthal angle
    phi_limits = [0, 2*np.pi]

    # Set-up integrator with multiple processors
    integ = vegas.Integrator([k_limits, theta_limits, phi_limits,
                              k_limits, theta_limits, phi_limits,
                              K_limits, theta_limits, phi_limits], nproc=8)

    # Loop over q_vector
    for i, q in enumerate(q_array):
            
        t0 = time.time()
        
        integrand = delta_U2_term_integrand(
            q, tau, quantum_numbers, cg_table, phi_functions, delta_U_functions,
            delta_U_dagger_functions
        )

        # # Train the integrator
        # integ(integrand, nitn=5, neval=200)
        # # Final result
        # result = integ(integrand, nitn=10, neval=1e3)
        
        # Train the integrator
        integ(integrand, nitn=5, neval=1e3)
        # Final result
        result = integ(integrand, nitn=10, neval=1e3)

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
    sp_basis = SingleParticleBasis(nucleus_name, Z, N, run_woodsaxon=False)
    phi_functions = get_sp_wave_functions(sp_basis, 10.0, 2.0, 120)
    
    # Set-up \delta U and \delta U^\dagger functions
    delta_U_functions, delta_U_dagger_functions = get_delta_U_functions(
        channels, kvnn, kmax, kmid, ntot, generator, lamb
    )
    
    # Set momentum mesh
    # q_array, q_weights = momentum_mesh(8.0, 2.0, 40)
    q_array, q_weights = momentum_mesh(10.0, 4.0, 100, nmod=60)
    
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
        delta_U_array, delta_U_errors = compute_delta_U_term(
            q_array, tau, sp_basis.occ_states, cg_table, channels,
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
            q_array, tau, sp_basis.occ_states, cg_table, channels,
            phi_functions, delta_U_functions, delta_U_dagger_functions
        )
        t4 = time.time()
        print(f"Done after {(t4-t3)/60:.3f} minutes.\n")
        
        # # # TESTING
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
    nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'O16', 8, 8
    
    # Nucleon
    tau = 1/2
    
    # Partial wave channels for expansion of plane-wave \delta U matrix elements
    channels = ('1S0', '3S1-3S1', '3S1-3D1', '3D1-3S1', '3D1-3D1')
    
    # NN potential and momentum mesh
    # kvnn, kmax, kmid, ntot = 6, 30.0, 4.0, 120  # AV18
    # kvnn, kmax, kmid, ntot = 111, 15.0, 3.0, 120  # SMS N4LO 450 MeV
    kvnn, kmax, kmid, ntot = 6, 15.0, 3.0, 120
    
    # SRG \lambda value
    # lamb = 1.35
    # lamb = 1.5
    lamb = 1.7

    # Compute and save the momentum distribution
    q_array, q_weights, n_array, n_errors = compute_momentum_distribution(
        nucleus_name, Z, N, tau, channels, kvnn, kmax, kmid, ntot, lamb,
        print_normalization=True, save=True
    )