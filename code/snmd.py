"""
File: snmd.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 17, 2022

The SingleNucleon class calculates single-nucleon momentum distributions (snmd)
with SRG-evolved operators assuming the evolved wave function is given by HF 
treated in LDA. This class is a sub-class of the MomentumDistribution class
from momentum_distributions.py

Last update: April 8, 2022

"""

# To-do: Make sure R_array parameter is described correctly.

# Python imports
import numpy as np

# Imports from A.T. codes
from modules.integration import gaussian_quadrature_mesh
from momentum_distributions import MomentumDistribution


class SingleNucleon(MomentumDistribution):
        
    def n_I(self, q_array, R_array, dR, kF1_array):
        """
        Evaluates the I term in U n(q) U^\dagger ~ I n(q) I.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from I term as a function of q.

        """
        
        # Initialize 2-D meshgrids (q, R)
        q_grid, R_grid = np.meshgrid(q_array, R_array, indexing='ij')
        
        # Get 2-D kF1(R) meshgrid
        _, kF1_grid = np.meshgrid(q_array, kF1_array, indexing='ij')
        
        # Evaluate the Heaviside step function in I term
        # theta_grid = self.theta_I(q_grid, kF1_grid)
        theta_grid = None
        
        # Calculate R integrand (ntot_q, ntot_R)
        integrand_R = theta_grid * R_grid**2 * dR

        # Integrate over R leaving a (ntot_q, 1) shape array
        # Factor of 2 is overall factor for summation over spin projections
        return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def n_deltaU(self, q_array, R_array, dR, kF1_array, kF2_array):
        """
        Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
        we are combining \delta U and \delta U^\dagger, hence the factor of 2.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U term as a function of q.

        """
        
        # Set number of k points for integration over k
        ntot_k = 40  # Typical integration goes up to kF ~ 1.3 fm^-1
        
        # Initialize 3-D meshgrids (q, R, k) with k values equal to 0
        q_grid, R_grid, k_grid = np.meshgrid(
            q_array, R_array, np.zeros(ntot_k), indexing='ij')
        
        # Get 3-D kF1(R) and kF2(R) meshgrids and initialize k weights mesh
        _, kF1_grid, dk_grid = np.meshgrid(
            q_array, kF1_array, np.zeros(ntot_k), indexing='ij')
        _, kF2_grid, _ = np.meshgrid(
            q_array, kF2_array, np.zeros(ntot_k), indexing='ij')

        # Loop over q and kF2 to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iq, q in enumerate(q_array):
            for ikF2, kF2 in enumerate(kF2_array):

                # Create integration meshgrid k_array up to (kF2 + q)/2 which
                # corresponds to the upper limit of \theta(kF2(R) - |q-2k|)
                k_max = (kF2 + q)/2
 
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max, ntot_k)
                
                # Fill in k_grid and dk_grid given the specific k_array
                k_grid[iq, ikF2, :] = k_array
                dk_grid[iq, ikF2, :] = k_weights
        
        # Evaluate angle-average with Heaviside step functions in \delta U 
        # term for \tau and \tau'
        # theta_pp_grid = self.theta_deltaU(q_grid, kF1_grid, kF1_grid, k_grid)
        # theta_pn_grid = self.theta_deltaU(q_grid, kF1_grid, kF2_grid, k_grid)
        theta_pp_grid = None
        theta_pn_grid = None
        
        # Evaluate \delta U(k,k) for \tau and \tau'
        deltaU_pp_grid = self.deltaU_pp_func(k_grid)
        deltaU_pn_grid = self.deltaU_pn_func(k_grid)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # Factor of 2 from \delta U + \delta U^\dagger
        # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        deltaU_factor = 2 * 8 * 2/np.pi * 2
        
        # Calculate the k integrand where we split terms according to pp and
        # pn (ntot_q, ntot_R, ntot_k)
        integrand_k = (deltaU_factor * k_grid**2 * dk_grid * R_grid**2 * dR
                       * (deltaU_pp_grid*theta_pp_grid
                          + deltaU_pn_grid*theta_pn_grid))
        
        # Integrate over k leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R leaving a (ntot_q, 1) shape array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
    
    def n_deltaU2(self, q_array, R_array, dR, kF1_array, kF2_array):
        """
        Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U \delta U^\dagger term as
            a function of q.

        """
        
        ntot_q = len(q_array)
        # Set number of K and k points for integration over K and k
        ntot_K = 40
        ntot_k = 40
        
        # y = cos(\theta) angles for averaging integration over angle between
        # K/2 and k
        ntot_y = 7
        y_array, y_weights = leggauss(ntot_y)
        
        # Initialize 4-D meshgrids (q, R, K, k) with k and K equal to 0
        q_grid, R_grid, K_grid, k_grid = np.meshgrid(
            q_array, R_array, np.zeros(ntot_K), np.zeros(ntot_k),
            indexing='ij')
        
        # Get 4-D kF1(R) and kF2(R) meshgrids and initialize K and k weights
        # meshgrids
        _, kF1_grid, dK_grid, dk_grid = np.meshgrid(
            q_array, kF1_array, np.zeros(ntot_K), np.zeros(ntot_k),
            indexing='ij')
        _, kF2_grid, _, _ = np.meshgrid(q_array, kF2_array, np.zeros(ntot_K),
                                        np.zeros(ntot_k), indexing='ij')

        # Loop over q, R, and k to find limits of K integration and then
        # create K integration mesh using Gaussian quadrature
        for iR, R in enumerate(R_array):
                
            kF1, kF2 = kF1_array[iR], kF2_array[iR]
    
            # K integration goes from 0 to kF1+kF2
            K_max = kF1 + kF2
                
            # Get Gaussian quadrature mesh for K integration
            K_array, K_weights = gaussian_quadrature_mesh(K_max, ntot_K)
              
            # Loop over remaining variables and fill in 4-D array
            for iq in range(ntot_q):
                for ik in range(ntot_k):
                    
                    # Fill in K_grid and dK_grid given the specific K_array
                    K_grid[iq, iR, :, ik] = K_array
                    dK_grid[iq, iR, :, ik] = K_weights

        # Loop over q, R, and K to find limits of k integration and then
        # create k_array using Gaussian quadrature
        for iR, R in enumerate(R_array):
        
            kF1, kF2 = kF1_array[iR], kF2_array[iR]
        
            # Get minimum kF value
            kF_min = min(kF1, kF2)
            
            # K_array only depends on R, so loop over K_grid[0, iR, :, 0]
            for iK, K in enumerate(K_grid[0, iR, :, 0]):
                
                # Lower limit of k integration
                k_min = max(K/2 - kF_min, 0)
                
                # Upper limit of k integration
                if K**2/4 < (kF1**2 + kF2**2)/2:
                    k_max = min(np.sqrt((kF1**2+kF2**2)/2-K**2/4), kF_min+K/2)
                else:
                    k_max = kF_min + K/2

                # Get Gaussian quadrature mesh for k integration
                k_array, k_weights = gaussian_quadrature_mesh(
                    k_max, ntot_k, xmin=k_min)
                
                # Loop over remaining variables and fill in 4-D meshgrid
                for iq in range(ntot_q):
                    
                    # Fill in k_grid and dk_grid given the specific k_array
                    k_grid[iq, iR, iK, :] = k_array
                    dk_grid[iq, iR, iK, :] = k_weights

        # Evaluate angle-average with Heaviside step functions in \delta U^2
        # term for \tau and \tau'
        # theta_pp_grid = self.theta_deltaU2(kF1_grid, kF1_grid, K_grid, k_grid)
        # theta_pn_grid = self.theta_deltaU2(kF1_grid, kF2_grid, K_grid, k_grid)
        theta_pp_grid = None
        theta_pn_grid = None

        # Set-up 4-D \delta U \delta U^\dagger(k, |q-K/2|) for \tau and \tau'
        deltaU2_pp_grid = np.zeros_like(theta_pp_grid)  # shape (q, R, K, k)
        deltaU2_pn_grid = np.zeros_like(theta_pn_grid)
        
        # Integrate over y
        for y, dy in zip(y_array, y_weights):
            
            # Evaluate |q-K/2| meshgrid
            q_K_grid = np.sqrt(q_grid**2 + K_grid**2/4 - q_grid*K_grid*y)
            
            deltaU2_pp_grid += self.deltaU2_pp_func.ev(k_grid, q_K_grid) * dy/2
            deltaU2_pn_grid += self.deltaU2_pn_func.ev(k_grid, q_K_grid) * dy/2
        
        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2

        # Calculate the k integrand where we split terms according to pp and
        # pn leaving a (ntot_q, ntot_R, ntot_K, ntot_k) shape array
        integrand_k = (deltaU2_factor * k_grid**2 * dk_grid * K_grid**2 
                       * dK_grid * R_grid**2 * dR
                       * (deltaU2_pp_grid*theta_pp_grid
                          +deltaU2_pn_grid*theta_pn_grid))
                      
        # Integrate over k leaving K integrand (ntot_q, ntot_R, ntot_K)
        integrand_K = np.sum(integrand_k, axis=-1)
        
        # Integrate over K leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_K, axis=-1)
        
        # Integrate over R leaving a (ntot_q, 1) shape array
        return 4*np.pi * np.sum(integrand_R, axis=-1)

#------------------------------------------------------------------------------
# File: snmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 25, 2021
# 
# Calculates SRG-evolved single-nucleon momentum distributions for nuclei
# assuming the evolved wave function is given by HF treated in LDA. (Note,
# 'snmd' stands for single-nucleon momentum distribution.)
#
# Warning: High momentum matrix elements of 3P2-3F2 and 3D3-3G3 channels are
# screwed up even with kmax=30 fm^-1 mesh. Ignore these channels for now!
#
# Revision history:
#   03/29/21 --- Finalized and produces same results (up to factors of 2*\pi)
#                as single_particle_momentum_dist.py. Also, fixed a bug
#                involving the integration over total momentum K.
#   03/31/21 --- Moved channel_L_value function to vnn.py. See vnn.py for
#                details of the function.
#   04/02/21 --- Added option to return pp and pn (or nn and np) contributions
#                to \delta U \delta U^\dagger term along with total.
#   04/22/21 --- Correcting overall factors in front of \delta U and
#                \delta U^2 terms.
#   04/27/21 --- Fixed issue with evaluating at maximum q-point in q_array.
#   05/25/21 --- Saving older version as snmd_v1.py with new updates: angle-
#                averaging, etc.
#   06/04/21 --- Including interpolating option to speed up code.
#   06/10/21 --- Replacing UnivariateSpline with interp1d for better accuracy,
#                though RectBivariateSpline works well for 2-D interpolations.
#   06/24/21 --- Speeding up code by switching from loops to np.sum() to do
#                integrations. Saved old version as snmd_v1.py in Old_codes.
#   07/15/21 --- Switched interpolation functions interp1d and 
#                RectBivariateSpline from cubic to linear since \delta U^2(k,q)
#                was returning negative values.
#   08/18/21 --- Added SRG block-diagonal option.
#   02/15/22 --- Added optional lambda_init argument to class. This allows one
#                to use an SRG-evolved potential as the initial interaction.
#   02/24/22 --- Added optional kvnn_hard argument to class. This allows one
#                to SRG-evolve the initial potential back to a harder one
#                using transformations of the kvnn_hard potential.
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d, RectBivariateSpline
# Scripts made by A.T.
# from densities import load_density
# from figures import figures_functions as ff
# from misc.integration import gaussian_quadrature_mesh
# from potentials.vsrg_macos import vnn
# from srg.srg_unitary_transformation import SRG_unitary_transformation
from .figure_labels import replace_periods
from integration import gaussian_quadrature_mesh
from .vnn import Potential
from srg_transformation import get_transformation
from tools import channel_L_value, coupled_channel


def load_density(nucleus_name, nucleon, Z, N, edf='SLY4'):
    """
    Loads a nucleonic density for the given nucleus. Densities are normalized
    according to
        4*\pi \int_0^\infty dR R^2 \rho_A(R) = Z or N.
    
    Parameters
    ----------
    nucleus_name : str
        Specify the nucleus (e.g., 'O16', 'Ca40', etc.)
    nucleon : str
        Specify 'proton' or 'neutron'.
    Z : int
        Proton number of the nucleus.
    N : int
        Neutron number of the nucleus.
    edf : str, optional
        Name of EDF (e.g., 'SLY4').
        
    Returns
    -------
    R_array : 1-D ndarray
        C.o.M. coordinates [fm].
    rho_array : 1-D ndarray
        Nucleonic density as a function of R [# of nucleons / vol].
                                              
    Notes
    -----
    Momentum distributions code compute intermediate integration arrays in 
    relative k and C.o.M. K which rely on kF(R) values. These values can be
    zero if the density \rho(R) = 0. We must replace zeros in \rho(R) with an
    extremely small number so the codes run correctly to avoid zero division
    errors. (This only happens for edf = 'AV18' densities.)
    
    """


    # Go to directory corresponding to specified nucleus and EDF
    if edf == 'SLY4':
        
        densities_directory = f'../../densities/HFBRAD_{edf}/{nucleus_name}/'
        file_name = f'{nucleon}_{N:d}_{Z:d}.dens'
        column_number = 1
        
    elif edf == 'Gogny':
        
        densities_directory = f'../../densities/{edf}/{nucleus_name}/'
        file_name = 'DensityQP.dat'
        if nucleon == 'proton':
            column_number = 1
        elif nucleon == 'neutron':
            column_number = 2
    
    # Technically it doesn't make sense to have a case edf == AV18 since
    # AV18 does not use an EDF. It would also make more sense to call it VMC.
    elif edf == 'AV18':
        
        densities_directory = '../../densities/{edf}/'
        file_name = '{nucleus_name}_densities_{N:d}_{Z:d}.txt'
        
        # AV18 files either have single \rho column for N=Z nuclei or
        # two columns for proton (1) and neutron (3)
        if N == Z:
            column_number = 1
        else:
            if nucleon == 'proton':
                column_number = 1
            elif nucleon == 'neutron':
                column_number = 3 
        
    # Load file
    table = np.loadtxt(densities_directory + file_name)
    
    R_array = table[:, 0]
    rho_array = table[:, column_number]
    
    # Avoiding zero division errors
    zero_case = rho_array == 0
    rho_array[zero_case] = 1e-30
    
    return R_array, rho_array


class single_nucleon_momentum_distributions(object):
    
    
    def __init__(
            self, kvnn, channels, lamb, kmax, kmid, ntot, generator='Wegner',
            interp=False, lambda_init=np.inf, kvnn_hard=0):
        """
        Evaluates and saves the pp and pn matrix elements of \delta U and
        \delta U^{\dagger} given the input potential and SRG \lambda.
        
        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        channels : tuple
            Partial wave channels to include in the calculation.
        lamb : float
            SRG evolution parameter \lambda [fm^-1].
        kmax : float
            Maximum value in the momentum mesh [fm^-1]. (Default of zero
            automatically selects default mesh based on kvnn.)
        kmid : float
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int
            Number of momentum points in mesh.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        interp : bool, optional
            Option to use interpolated n_\lambda(q) functions.
        lambda_init : float, optional
            SRG evolution parameter \lambda for initial Hamiltonian [fm^-1].
            This allows one to use an SRG-evolved potential as the starting
            point.
        kvnn_hard : int, optional
            Inputing a nonzero argument here will evolve the initial potential
            (corresponding to kvnn) back to a harder scale using 
            transformations from a harder potential (kvnn_hard).
            
        """
        
        # Get relevant info for file and directory names
        # Part of data directory name
        self.kvnn = kvnn 
        # Part of file name
        self.channels = channels
        self.lamb = lamb
        self.kmax = kmax
        self.generator = generator
        
        # Set-up for calculation
        if interp == False:

            # Save highest allowed L based on input channels
            highest_L = 0
            for channel in channels:
                next_L = channel_L_value(channel)
                if next_L > highest_L:
                    highest_L = next_L

            # Load and save momentum arrays
            k_array, k_weights = Potential(kvnn, '1S0', kmax, kmid, ntot).load_mesh()
            # Save k_array for writing files
            self.k_array = k_array
            
            # For dividing out momenta/weights
            factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
            # For coupled-channel matrices
            factor_array_cc = np.concatenate( (factor_array, factor_array) )

            # Initialize pp and pn matrix elements
            deltaU_pp = np.zeros( (ntot, ntot) ) # \delta U linear term
            deltaU_pn = np.zeros( (ntot, ntot) )
            deltaU2_pp = np.zeros( (ntot, ntot) ) # \delta U \delta U^\dagger
            deltaU2_pn = np.zeros( (ntot, ntot) )
        
            # Allowed channels for pp (and nn) up through the D-waves
            pp_channels = ('1S0', '3P0', '3P1', '3P2', '1D2')
        
            # Loop over channels and evaluate matrix elements
            for channel in channels:
                
                potential = Potential(kvnn, channel, kmax, kmid, ntot)

                # Load SRG transformation
                
                # If lambda_init = np.inf, then take the initial Hamiltonian
                # of kvnn as the starting point
                if lambda_init == np.inf:
                    
                    H_initial = potential.load_hamiltonian()
                    
                # Otherwise, this splits into two cases
                else:
                    
                    # If kvnn_hard is zero, take an SRG-evolved Hamiltonian
                    # corresponding to kvnn and lambda_init as the initial
                    # Hamiltonian
                    if kvnn_hard == 0:
                        
                        # Block-diagonal
                        if generator == 'Block-diag':
                            H_initial = potential.load_hamiltonian(
                                'srg', generator, 1.0, lambda_init)
                            
                        # Band-diagonal
                        else:
                            H_initial = potential.load_hamiltonian(
                                'srg', generator, lambda_init)
                    
                    # If kvnn_hard is nonzero, SRG-evolve the initial 
                    # Hamiltonian of kvnn back using transformations from
                    # kvnn_hard at lambda_initial
                    else:
                        
                        potential_hard = Potential(kvnn_hard, channel, kmax, kmid, ntot)
                        
                        H_hard_initial = potential_hard.load_hamiltonian()
                        
                        # Block-diagonal
                        if generator == 'Block-diag':
                            H_hard_evolved = potential_hard.load_hamiltonian(
                                'srg', generator, 1.0, lambda_init)
                            
                        # Band-diagonal
                        else:
                            H_hard_evolved = potential_hard.load_hamiltonian(
                                'srg', generator, lambda_init)
                            
                        # Get SRG transformation from hard potential
                        U_hard = get_transformation(H_hard_initial,
                                                    H_hard_evolved)
                        
                        # Get initial Hamiltonian for kvnn
                        H_matrix = potential.load_hamiltonian()
                        
                        # Do inverse transformation on softer Hamiltonian
                        H_initial = U_hard.T @ H_matrix @ U_hard
                    
                if generator == 'Block-diag':
                    # Take \lambda = 1 fm^-1 and set \Lambda_BD = input \lambda
                    H_evolved = potential.load_hamiltonian('srg', generator,
                                                           1.0, lamb)
                else:
                    H_evolved = potential.load_hamiltonian('srg', generator,
                                                           lamb)
                # Load U(k, k') [unitless]
                U_matrix_unitless = get_transformation(H_initial, H_evolved)

                # Isolate 2-body term and convert to fm^3
                if coupled_channel(channel):
                    I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
                    row, col = np.meshgrid(factor_array_cc, factor_array_cc)
                else:
                    I_matrix_unitless = np.eye(ntot, ntot)
                    row, col = np.meshgrid(factor_array, factor_array)
                delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
                delta_U_matrix = delta_U_matrix_unitless / row / col # fm^3
            
                # 2J+1 factor
                J = int( channel[-1] )
            
                # Add to the pp and pn terms
                # Coupled-channel
                if coupled_channel(channel):
                    
                    # First L of coupled-channel
                    # Isospin CG's=1/\sqrt(2) for pn
                    deltaU_pn += (2*J+1)/2 * delta_U_matrix[:ntot, :ntot]
                    deltaU2_pn += (2*J+1)/2 * ( \
                                      delta_U_matrix[:ntot, :ntot]**2 + \
                                      delta_U_matrix[:ntot, ntot:]**2 )

                    # Isospin CG's=1 for pp
                    if channel in pp_channels:
                        deltaU_pp += (2*J+1) * delta_U_matrix[:ntot, :ntot]
                        deltaU2_pp += (2*J+1) * ( \
                                          delta_U_matrix[:ntot, :ntot]**2 \
                                        + delta_U_matrix[:ntot, ntot:]**2 )
                    
                    # Decide whether to add second L based on highest allowed
                    # L value (e.g., 0 + 2 <= 2 meaning we include the 3D1-3D1
                    # part of the coupled 3S1-3D1 channel if we input D-waves
                    # in channels)
                    if channel_L_value(channel) + 2 <= highest_L:
                        deltaU_pn += (2*J+1)/2 * delta_U_matrix[ntot:, ntot:]
                        deltaU2_pn += (2*J+1)/2 * ( \
                                          delta_U_matrix[ntot:, :ntot]**2 + \
                                          delta_U_matrix[ntot:, ntot:]**2 )
                        
                        if channel in pp_channels:
                            deltaU_pp += (2*J+1) * delta_U_matrix[ntot:, ntot:]
                            deltaU2_pp += (2*J+1) * ( \
                                            delta_U_matrix[ntot:, :ntot]**2 + \
                                            delta_U_matrix[ntot:, ntot:]**2 )
            
                else:
                
                    # Isospin CG's=1/\sqrt(2) for pn
                    deltaU_pn += (2*J+1)/2 * delta_U_matrix
                    deltaU2_pn += (2*J+1)/2 * delta_U_matrix**2
                
                    # Isospin CG's=1 for pp
                    if channel in pp_channels:
                        deltaU_pp += (2*J+1) * delta_U_matrix
                        deltaU2_pp += (2*J+1) * delta_U_matrix**2

            # Interpolate pp and pn < k | \delta U | k >
            self.deltaU_pp_func = interp1d( k_array, np.diag(deltaU_pp),
                                            kind='linear', bounds_error=False,
                                            fill_value='extrapolate' )
            self.deltaU_pn_func = interp1d( k_array, np.diag(deltaU_pn),
                                            kind='linear', bounds_error=False,
                                            fill_value='extrapolate' )
        
            # Interpolate pp and pn < k | \delta U \delta U^{\dagger} | k' > 
            self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                       deltaU2_pp, kx=1, ky=1)
            self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                       deltaU2_pn, kx=1, ky=1)


    def theta_I(self, q_mesh, kF1_mesh):
        """
        Evaluates \theta( kF1(R) - q ). This function appears in the I term.

        Parameters
        ----------
        q_mesh : 2-D ndarray
            Momentum values [fm^-1].
        kF1_mesh : 2-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.

        Returns
        -------
        theta_mesh : 2-D ndarray
            \theta function [unitless] evaluated for each q and kF1(R).

        """
        
        # Initialize 2-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_R) )
        
        # Gives 1 if q < kF1(R)
        theta_mesh[ q_mesh < kF1_mesh ] = 1
        
        # This is a (ntot_q, ntot_R) size array
        return theta_mesh
        

    def theta_deltaU(self, q_mesh, kF1_mesh, kF2_mesh, k_mesh):
        """
        Evaluates angle-average of \theta( kF2(R) - \abs(q - 2k) ). This
        function appears in the \delta U term.

        Parameters
        ----------
        q_mesh : 3-D ndarray
            Momentum values [fm^-1].
        kF1_mesh : 3-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_mesh : 3-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.
        k_mesh : 3-D ndarray
            Relative momentum values [fm^-1].

        Returns
        -------
        theta_mesh : 3-D ndarray
            \theta function [unitless] evaluated for each q, kF1(R), kF2(R),
            and k.
            
        Notes
        -----
        Not sure why the cases had to be computed in reverse order. Does that
        mean there is an overlap of truth values in case 3 and case 1 when
        there shouldn't be?

        """
        
        # Initialize 3-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_R, self.ntot_k) )
        
        # Evaluate each boolean case and use these to fill in the theta_mesh
        
        # This applies to each case: q < kF1
        case_all = q_mesh < kF1_mesh
        
        # Case 3: q-kF < 2k < kF+q
        case_3 = case_all * ( q_mesh - kF2_mesh < 2*k_mesh ) * \
                 ( 2*k_mesh < kF2_mesh + q_mesh )
        theta_mesh[case_3] = ( ( kF2_mesh**2 - ( q_mesh - 2*k_mesh )**2 ) / \
                               ( 8*k_mesh*q_mesh ) )[case_3]
            
        # Case 2: q < kF and kF-q <= 2k < kF+q
        case_2 = case_all * ( q_mesh < kF2_mesh ) * \
                 ( kF2_mesh - q_mesh <= 2*k_mesh ) * \
                 ( 2*k_mesh < kF2_mesh + q_mesh )
        theta_mesh[case_2] = ( ( kF2_mesh**2 - ( q_mesh - 2*k_mesh )**2 ) / \
                               ( 8*k_mesh*q_mesh ) )[case_2]
        
        # Case 1: q < kF and 2k < kF-q
        case_1 = case_all * ( q_mesh < kF2_mesh ) * \
                 ( 2*k_mesh < kF2_mesh - q_mesh )
        theta_mesh[case_1] = 1

        # This is a (ntot_q, ntot_R, ntot_k) size array
        return theta_mesh


    def theta_deltaU2(self, kF1_mesh, kF2_mesh, K_mesh, k_mesh):
        """
        Evaluates angle-average of \theta( kF1(R) - \abs(K/2 + k) ) x
        \theta( kF2(R) - \abs(K/2 - k) ). This function appears in the
        \delta U \delta U^\dagger term.

        Parameters
        ----------
        kF1_mesh : 4-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_mesh : 4-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.
        K_mesh : 4-D ndarray
            C.o.M. momentum values [fm^-1].
        k_mesh : 4-D ndarray
            Relative momentum values [fm^-1].

        Returns
        -------
        theta_mesh : 4-D ndarray
            \theta function [unitless] evaluated for each q, kF1(R), kF2(R), K,
            and k. Note, this function does not depend on q but we need the
            array to match the size of \delta U \delta U^\dagger matrix
            elements.

        """
        
        # Initialize 4-D array
        theta_mesh = np.zeros( (self.ntot_q, self.ntot_R, self.ntot_K,
                                self.ntot_k) )
        
        # Evaluate each boolean case and use these to fill in the theta_mesh
        
        # Case 4: 2k+K < 2kF1 and -4 < (4k^2 - 4kF2^2 + K^2)/(kK) < 4
        case_4 = ( 2*k_mesh + K_mesh <= 2*kF1_mesh ) * \
                 ( -4 < (4*k_mesh**2-4*kF2_mesh**2+K_mesh**2) / \
                 (k_mesh*K_mesh) ) * \
                 ( (4*k_mesh**2-4*kF2_mesh**2+K_mesh**2) / \
                 (k_mesh*K_mesh) <= 4 )
        theta_mesh[case_4] = ( ( 4*kF2_mesh**2 - (K_mesh-2*k_mesh)**2 ) / \
                               (8*k_mesh*K_mesh) )[case_4]
            
        # Case 3: 2k+K < 2kF2 and -4 < (4k^2 - 4kF1^2 + K^2)/(kQ) < 4
        case_3 = ( 2*k_mesh + K_mesh <= 2*kF2_mesh ) * \
                 ( -4 < (4*k_mesh**2-4*kF1_mesh**2+K_mesh**2) / \
                 (k_mesh*K_mesh) ) * \
                 ( (4*k_mesh**2-4*kF1_mesh**2+K_mesh**2) / \
                 (k_mesh*K_mesh) <= 4 )
        theta_mesh[case_3] = ( ( 4*kF1_mesh**2 - (K_mesh-2*k_mesh)**2 ) / \
                               (8*k_mesh*K_mesh) )[case_3]
            
        # Case 2: 2k+K > 2kF1 and 2k+K > 2kF2 and 4k^2+K^2 < 2(kF1^2+kF2^2)
        case_2 = ( 2*k_mesh + K_mesh > 2*kF1_mesh ) * \
                 ( 2*k_mesh + K_mesh > 2*kF2_mesh ) * \
                 ( 4*k_mesh**2 + K_mesh**2 <= 2*(kF1_mesh**2 + kF2_mesh**2) )    
        theta_mesh[case_2] = ( ( 2*(kF1_mesh**2+kF2_mesh**2) - 4*k_mesh**2 - \
                                 K_mesh**2 ) / (4*k_mesh*K_mesh) )[case_2]

        # Case 1: 2k+K < 2kF1 and 2k+K < 2kF2
        case_1 = ( 2*k_mesh + K_mesh <= 2*kF1_mesh ) * \
                 ( 2*k_mesh + K_mesh <= 2*kF2_mesh )
        theta_mesh[case_1] = 1
        
        # This is a (ntot_q, ntot_R, ntot_K, ntot_k) size array
        return theta_mesh
    
    
    def n_I(self, q_array, R_array, dR, kF1_array):
        """
        Evaluates the I term in U n(q) U^\dagger ~ I n(q) I.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from I term as a function of q.

        """
        
        # Initialize 2-D meshgrids (q, R)
        q_mesh, R_mesh = np.meshgrid(q_array, R_array, indexing='ij')
        
        # Get 2-D kF1(R) mesh
        _, kF1_mesh = np.meshgrid(q_array, kF1_array, indexing='ij')
        
        # Evaluate the \theta-function in I term
        theta_mesh = self.theta_I(q_mesh, kF1_mesh)
        
        # Calculate R integrand (ntot_q, ntot_R)
        integrand_R = theta_mesh * R_mesh**2 * dR

        # Integrate over R
        # This is a (ntot_q, 1) size array
        # Factor of 2 is overall factor
        return 2 * 4*np.pi * np.sum(integrand_R, axis=-1)
    
    
    def n_deltaU(self, q_array, R_array, dR, kF1_array, kF2_array):
        """
        Evaluates second and third terms in U n(q) U^\dagger ~ \delta U. Here
        we are combining \delta U and \delta U^\dagger, hence the factor of 2.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U term as a function of q.

        """
        
        # Set number of k points for integration over k
        self.ntot_k = 40 # Typical integration goes up to kF ~ 1.3 fm^-1
        
        # Initialize 3-D meshgrids (q, R, k) with k values equal to 0
        q_mesh, R_mesh, k_mesh = np.meshgrid(q_array, R_array,
                                             np.zeros(self.ntot_k),
                                             indexing='ij')
        
        # Get 3-D kF1(R) and kF2(R) meshes and initialize k weights mesh
        _, kF1_mesh, dk_mesh = np.meshgrid(q_array, kF1_array,
                                           np.zeros(self.ntot_k),
                                           indexing='ij')
        _, kF2_mesh, _ = np.meshgrid(q_array, kF2_array, np.zeros(self.ntot_k),
                                     indexing='ij')

        # Loop over q and R to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iq, q in enumerate(q_array):
            for iR, R in enumerate(R_array):
                
                kF2 = kF2_array[iR]

                # Create integration mesh k_array up to (kF2 + q)/2 which
                # corresponds to the upper limit of \theta( kF2(R) - |q-2k| )
                k_max = (kF2 + q)/2
 
                # Get Gaussian quadrature mesh
                k_array, k_weights = gaussian_quadrature_mesh(k_max,
                                                              self.ntot_k)
                
                # Fill in k_mesh and dk_mesh given the specific k_array
                k_mesh[iq, iR, :] = k_array
                dk_mesh[iq, iR, :] = k_weights
        
        # Evaluate angle-average of \theta-functions in \delta U term for \tau
        # and \tau'
        theta_pp_mesh = self.theta_deltaU(q_mesh, kF1_mesh, kF1_mesh, k_mesh)
        theta_pn_mesh = self.theta_deltaU(q_mesh, kF1_mesh, kF2_mesh, k_mesh)
        
        # Evaluate < k | \delta U | k > for \tau and \tau'
        deltaU_pp_mesh = self.deltaU_pp_func(k_mesh)
        deltaU_pn_mesh = self.deltaU_pn_func(k_mesh)

        # Contractions of a's, 1/4 factors, and [1-(-1)^(L+S+T)] factors
        # combine to give 2
        # Factor of 2 from \delta U + \delta U^\dagger
        # Factor of 8 is from evaluating \int d^3K \delta(K/2 - ...)
        # 2/\pi for two | k_vec > -> | k J L S ... > changes
        deltaU_factor = 2 * 8 * 2/np.pi * 2
        
        # Calculate the k integrand where we split terms according to pp and
        # pn (or nn and np if kF_1 corresponds to a neutron)
        # (ntot_q, ntot_R, ntot_k)
        integrand_k = deltaU_factor * k_mesh**2 * dk_mesh * R_mesh**2 * dR * (\
                      deltaU_pp_mesh * theta_pp_mesh + \
                      deltaU_pn_mesh * theta_pn_mesh )
        
        # Integrate over k leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_k, axis=-1)

        # Integrate over R
        # This is a (ntot_q, 1) size array
        return 4*np.pi * np.sum(integrand_R, axis=-1)


    def n_deltaU2(self, q_array, R_array, dR, kF1_array, kF2_array):
        """
        Evaluates fourth term in U n(q) U^\dagger ~ \delta U \delta U^\dagger.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        kF1_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau with
            respect to R.
        kF2_array : 1-D ndarray
            Fermi momentum [fm^-1] for the nucleon corresponding to \tau' with
            respect to R.

        Returns
        -------
        output : 1-D ndarray
            Momentum distribution [fm^3] from \delta U \delta U^\dagger term as
            a function of q.

        """
        
        # Set number of K and k points for integration over K and k
        self.ntot_K = 40
        self.ntot_k = 40
        
        # y = cos(\theta) angles for averaging integration
        ntot_y = 7
        y_array, y_weights = leggauss(ntot_y)
        
        # Initialize 4-D meshgrids (q, R, K, k) with k and K equal to 0
        q_mesh, R_mesh, K_mesh, k_mesh = np.meshgrid(q_array, R_array,
                                                     np.zeros(self.ntot_K),
                                                     np.zeros(self.ntot_k),
                                                     indexing='ij')
        
        # Get 4-D kF1(R) and kF2(R) meshes and initialize K and k weights mesh
        _, kF1_mesh, dK_mesh, dk_mesh = np.meshgrid(q_array, kF1_array,
                                                    np.zeros(self.ntot_K),
                                                    np.zeros(self.ntot_k),
                                                    indexing='ij')
        _, kF2_mesh, _, _ = np.meshgrid(q_array, kF2_array,
                                        np.zeros(self.ntot_K),
                                        np.zeros(self.ntot_k), indexing='ij')

        # Loop over q, R, and k to find limits of K integration and then create
        # K_array using Gaussian quadrature
        for iR, R in enumerate(R_array):
                
            kF1, kF2 = kF1_array[iR], kF2_array[iR]
    
            # K integration goes from 0 to kF1+kF2
            K_max = kF1 + kF2
                
            # Get Gaussian quadrature mesh for K integration
            K_array, K_weights = gaussian_quadrature_mesh(K_max, self.ntot_K)
              
            # Loop over remaining variables and fill in 4-D array
            for iq in range(self.ntot_q):
                for ik in range(self.ntot_k):
                    
                    # Fill in K_mesh and dK_mesh given the specific K_array
                    K_mesh[iq, iR, :, ik] = K_array
                    dK_mesh[iq, iR, :, ik] = K_weights

        # Loop over q, R, and K to find limits of k integration and then create
        # k_array using Gaussian quadrature
        for iR, R in enumerate(R_array):
        
            kF1, kF2 = kF1_array[iR], kF2_array[iR]
        
            # Get minimum kF value
            kF_min = min(kF1, kF2)
            
            # K_array only depends on R, so loop over K_mesh[0, iR, :, 0]
            for iK, K in enumerate( K_mesh[0, iR, :, 0] ):
                
                # Lower limit of k integration
                k_min = max(K/2 - kF_min, 0)
                
                # Upper limit of k integration
                if K**2/4 < (kF1**2 + kF2**2)/2:
                    k_max = min( np.sqrt( ( kF1**2 + kF2**2 )/2 - K**2/4 ),
                                 kF_min + K/2 )
                else:
                    k_max = kF_min + K/2

                # Get Gaussian quadrature mesh for k integration
                k_array, k_weights = gaussian_quadrature_mesh(k_max,
                                                              self.ntot_k,
                                                              xmin=k_min)
                
                # Loop over remaining variables and fill in 4-D array    
                for iq in range(self.ntot_q):
                    
                    # Fill in k_mesh and dk_mesh given the specific k_array
                    k_mesh[iq, iR, iK, :] = k_array
                    dk_mesh[iq, iR, iK, :] = k_weights
                        
        # Evaluate angle-average of \theta-functions in \delta U^2 term for
        # \tau and \tau'
        theta_pp_mesh = self.theta_deltaU2(kF1_mesh, kF1_mesh, K_mesh, k_mesh)
        theta_pn_mesh = self.theta_deltaU2(kF1_mesh, kF2_mesh, K_mesh, k_mesh)

        # Evaluate 4-D < k | \delta U | |q_vec-K_vec/2| >^2 for \tau and \tau'
        # while also averaging over angle y
        deltaU2_pp_mesh = np.zeros_like(theta_pp_mesh) # (q, R, K, k)
        deltaU2_pn_mesh = np.zeros_like(theta_pn_mesh)
        
        # Integrate over y
        for y, dy in zip(y_array, y_weights):
            
            # Evaluate |q-K/2| mesh
            q_K_mesh = np.sqrt( q_mesh**2 + K_mesh**2/4 - q_mesh*K_mesh*y )
            
            deltaU2_pp_mesh += self.deltaU2_pp_func.ev(k_mesh, q_K_mesh) * dy/2
            deltaU2_pn_mesh += self.deltaU2_pn_func.ev(k_mesh, q_K_mesh) * dy/2
        
        # Contractions of a's, 1/4 factors, and [ 1 - (-1)^(L+S+T) ] factors
        # combine to give 2
        # (2/\pi)^2 for four | k_vec > -> | k J L S ... > changes
        deltaU2_factor = 2 * (2/np.pi)**2

        # Calculate the k integrand where we split terms according to pp and
        # pn (or nn and np if kF_1 corresponds to a neutron)
        # (ntot_q, ntot_R, ntot_K, ntot_k)
        integrand_k = deltaU2_factor * k_mesh**2 * dk_mesh * K_mesh**2 * \
                      dK_mesh * R_mesh**2 * dR * ( \
                      deltaU2_pp_mesh * theta_pp_mesh + \
                      deltaU2_pn_mesh * theta_pn_mesh )

        # Integrate over k leaving K integrand (ntot_q, ntot_R, ntot_K)
        integrand_K = np.sum(integrand_k, axis=-1)
        
        # Integrate over K leaving R integrand (ntot_q, ntot_R)
        integrand_R = np.sum(integrand_K, axis=-1)
        
        # Integrate over R
        # This is a (ntot_q, 1) size array
        return 4*np.pi * np.sum(integrand_R, axis=-1)
            
    
    def n_total(self, q_array, R_array, dR, rho_1_array, rho_2_array):
        """
        Single-nucleon momentum distribution where the nucleonic densities
        specify the nucleus and distribution type (e.g., O16 and p).

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the nucleon corresponding
            to \tau.
        rho_2_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the nucleon corresponding
            to \tau'.

        Returns
        -------
        n_total : 1-D ndarray
            Single-nucleon momentum distribution [fm^3] for each q.

        """
        
        # Save lengths of q_array and R_array
        self.ntot_q = len(q_array)
        self.ntot_R = len(R_array)

        # Evaluate kF values at each point in R_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
            
        # Get each contribution with respect to q
        n_I = self.n_I(q_array, R_array, dR, kF1_array)
        n_deltaU = self.n_deltaU(q_array, R_array, dR, kF1_array, kF2_array)
        n_deltaU2 = self.n_deltaU2(q_array, R_array, dR, kF1_array, kF2_array)

        # Return total (ntot_q, 1)
        return n_I + n_deltaU + n_deltaU2
    
    
    def n_contributions(self, q_array, R_array, dR, rho_1_array, rho_2_array):
        """
        Contributions to the single-nucleon momentum distribution where the
        nucleonic densities specify the nucleus and distribution type (e.g.,
        O16 and p). This function isolates the I, \delta U, and
        \delta U \delta U^\dagger terms.

        Parameters
        ----------
        q_array : 1-D ndarray
            Momentum values [fm^-1].
        R_array : 1-D ndarray
            C.o.M. coordinates [fm].
        dR : float
            C.o.M. coordinates step-size (assuming linearly-spaced array) [fm].
        rho_1_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the nucleon corresponding
            to \tau.
        rho_2_array : 1-D ndarray
            Densities as a function of R [fm^-3] for the nucleon corresponding
            to \tau'.

        Returns
        -------
        n_contributions : tuple
            Tuple of 1-D ndarrays corresponding to the I, \delta U, and
            \delta U^\dagger terms of the pair momentum distribution [fm^3] for
            each q.

        """
        
        # Save lengths of q_array and R_array
        self.ntot_q = len(q_array)
        self.ntot_R = len(R_array)

        # Evaluate kF values at each point in R_array
        kF1_array = (3*np.pi**2 * rho_1_array)**(1/3)
        kF2_array = (3*np.pi**2 * rho_2_array)**(1/3)
            
        # Get each contribution with respect to q
        n_I = self.n_I(q_array, R_array, dR, kF1_array)
        n_deltaU = self.n_deltaU(q_array, R_array, dR, kF1_array, kF2_array)
        n_deltaU2 = self.n_deltaU2(q_array, R_array, dR, kF1_array, kF2_array)
        
        # Return tuple of contributions ( (ntot_q, ntot_Q), ... )
        return n_I, n_deltaU, n_deltaU2
    
    
    def write_file(self, nucleus, nucleon, Z, N, edf='SLY4'):
        """
        Write single-nucleon momentum distribution file for interpolation
        purposes. Split things into total, I, \delta U, and \delta U^2
        contributions.

        Parameters
        ----------
        nucleus : str
            Specify nucleus (e.g., 'O16', 'Ca40', etc.)
        nucleon : str
            Specify 'proton' or 'neutron'.
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
        edf : str, optional
            Name of EDF (e.g., 'SLY4').

        """
        
        # Get momentum values
        q_array = self.k_array
        
        # Directory for distributions data
        data_directory = 'data/snmd/kvnn_%d/%s' % (self.kvnn, edf)
        
        # Create file name
        file_name = '%s_%s_channels' % (nucleus, nucleon)
        # Add each channel to file name
        for channel in self.channels:
            file_name += '_%s' % channel
        if self.generator == 'Block-diag':
            file_name += '_LambdaBD_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        else:
            file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = replace_periods(file_name) + '.dat'
        
        # Load R values and nucleonic densities (the R_array's are the same)
        R_array, rho_p_array = load_density(nucleus, 'proton', Z, N, edf)
        if edf == 'AV18' and Z == N: # e.g., AV18 He4 densities
            rho_n_array = rho_p_array
        else: # e.g., AV18 He8 densities
            R_array, rho_n_array = load_density(nucleus, 'neutron', Z, N, edf)
        dR = R_array[2] - R_array[1] # Assuming linear spacing

        # Calculate n_\lambda^\tau(q) for each q in q_array
        if nucleon == 'proton':
            n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
                                q_array, R_array, dR, rho_p_array, rho_n_array)
        elif nucleon == 'neutron':
            n_I_array, n_delU_array, n_delU2_array = self.n_contributions(
                                q_array, R_array, dR, rho_n_array, rho_p_array)
            
        # Total momentum distribution
        n_total_array = n_I_array + n_delU_array + n_delU2_array
    
        # Open file and write header where we allocate roughly 18 centered
        # spaces for each label
        f = open(data_directory + '/' + file_name, 'w')
        header = '#' + '{:^17s}{:^18s}{:^18s}{:^18s}{:^18s}'.format('q',
                  'total', '1', '\delta U', '\delta U^2')
        f.write(header + '\n')
    
        # Loop over momenta q
        for iq, q in enumerate(q_array):

            # Write to data file following the format from the header
            line = '{:^18.6f}{:^18.6e}{:^18.6e}{:^18.6e}{:^18.6e}'.format( q,
                            n_total_array[iq], n_I_array[iq], n_delU_array[iq],
                            n_delU2_array[iq] )
            f.write('\n' + line)

        # Close file
        f.close()
        
        
    def n_lambda_interp(self, nucleus, nucleon, Z, N, edf='SLY4'):
        """
        Interpolate the single-nucleon momentum distribution for the specified
        file.

        Parameters
        ----------
        nucleus : str
            Specify nucleus (e.g., 'O16', 'Ca40', etc.)
        nucleon : str
            Specify 'proton' or 'neutron'.
        Z : int
            Proton number of the nucleus.
        N : int
            Neutron number of the nucleus.
        edf : str, optional
            Name of EDF (e.g., 'SLY4').
            
        Returns
        -------
        output : tuple
            Tuple of functions that depend only on momentum q [fm^-1] where
            each function corresponds to contributions to n_\lambda(q): total,
            1, \delta U, and \delta U^2.

        """
        
        # Directory for distributions data
        data_directory = 'data/snmd/kvnn_%d/%s' % (self.kvnn, edf)
        
        # Get file name
        file_name = '%s_%s_channels' % (nucleus, nucleon)
        # Add each channel to file name
        for channel in self.channels:
            file_name += '_%s' % channel
        if self.generator == 'Block-diag':
            file_name += '_LambdaBD_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        else:
            file_name += '_lamb_%.2f_kmax_%.1f' % (self.lamb, self.kmax)
        file_name = replace_periods(file_name) + '.dat'
        
        # Load data which includes all contributions to n_\lambda(q)
        data = np.loadtxt(data_directory + '/' + file_name)
        
        # Split data into 1-D arrays for each column
        q_array = data[:, 0] # Momentum in fm^-1
        n_total_array = data[:, 1] # Total distribution
        n_I_array = data[:, 2] # 1 term
        n_delU_array = data[:, 3] # \delta U term
        n_delU2_array = data[:, 4] # \delta U^2 term
        
        # Interpolate each array (UnivariateSpline is for smoothing whereas
        # interp1d gives closer value to the actual calculation)
        n_total_func = interp1d(q_array, n_total_array, bounds_error=False,
                                kind='linear', fill_value='extrapolate')
        n_I_func = interp1d(q_array, n_I_array, bounds_error=False,
                            kind='linear', fill_value='extrapolate')
        n_delU_func = interp1d(q_array, n_delU_array, bounds_error=False,
                               kind='linear', fill_value='extrapolate')
        n_delU2_func = interp1d(q_array, n_delU2_array, bounds_error=False,
                                kind='linear', fill_value='extrapolate')
        
        # Return all contributions with total first
        # Note, these are functions of q
        return n_total_func, n_I_func, n_delU_func, n_delU2_func