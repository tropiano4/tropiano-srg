#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: snmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 25, 2021
# 
# Calculates SRG-evolved single-nucleon momentum distributions for nuclei
# assuming the evolved wave function is given by a free Fermi gas. (Note,
# 'snmd' stands for single-nucleon momentum distribution.) Combine these
# functions with LDA.py for nuclear-averaged momentum distributions. This
# script is based off single_particle_momentum_dist.py in Old_codes.
#
# Revision history:
#   XX/XX/XX --- ...
#
#------------------------------------------------------------------------------


import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
from integration import gaussian_quadrature_mesh
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


# --- Tests and bug checks --- #
# 1. How much do things change with different K and x meshes?


class single_nucleon_momentum_distributions(object):
    
    
    def __init__(self, kvnn, channels, lamb, kmax=0.0, kmid=0.0, ntot=0):
        """
        Saves momentum arrays, grids, and matrix elements of \delta U and
        \delta U^{\dagger} given the input potential and SRG \lambda.
        Evaluates and saves the pp and pn matrix elements.
        
        Parameters
        ----------
        kvnn : int
            This number specifies the potential.
        channels : tuple
            Partial wave channels to include in the calculation.
        lamb : float
            SRG evolution parameter lambda [fm^-1].
        kmax : float, optional
            Maximum value in the momentum mesh [fm^-1].
        kmid : float, optional
            Mid-point value in the momentum mesh [fm^-1].
        ntot : int, optional
            Number of momentum points in mesh.
            
        """
        
        # --- Set up --- #
        
        # Save highest allowed L based on input channels (assuming channels
        # is ordered correctly!)
        highest_L = 0
        for channel in channels:
            next_L = self.channel_L_value(channel)
            if next_L > highest_L:
                highest_L = next_L
        
        # Load and save momentum and angular arrays for integration
        
        # Relative momentum k [fm^-1]
        k_array, k_weights = vnn.load_momentum(kvnn, '1S0')
        ntot = len(k_array)
        self.k_array, self.k_weights, self.ntot = k_array, k_weights, ntot
        
        # Total momentum K [fm^-1] where we put more points toward K=0 fm^-1
        Kmax = 3.0 # Max momentum
        Kmid = 1.5 # Mid-point
        Ntot = 20 # Total number of points
        Nmod = 10 # Total number of points in the low-K region
        K_array, K_weights = gaussian_quadrature_mesh(Kmax, Ntot, xmid=Kmid,
                                                      nmod=Nmod)
        self.K_array, self.K_weights, self.Ntot = K_array, K_weights, Ntot
        
        # cos(\theta) angles for averaging
        xtot = 9
        x_array, x_weights = leggauss(xtot)
        self.x_array, self.x_weights, self.xtot = x_array, x_weights, xtot
        
        # For dividing out momenta/weights
        factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
        # For coupled-channel matrices
        factor_array_cc = np.concatenate( (factor_array, factor_array) )
        
        
        # --- Evaluate matrix elements for each channel --- #
        
        # Initialize pp and pn matrix elements
        deltaU_pp = np.zeros( (ntot, ntot) ) # \delta U linear term
        deltaU_pn = np.zeros( (ntot, ntot) )
        deltaU2_pp = np.zeros( (ntot, ntot) ) # \delta U \delta U^\dagger term
        deltaU2_pn = np.zeros( (ntot, ntot) )
        
        # Allowed channels for pp (and nn) up through the D-waves
        pp_channels = ['1S0', '3P0', '3P1', '3P2', '1D2']
        
        # Loop over channels and evaluate matrix elements
        for channel in channels:

            # Load SRG transformation
            H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
            H_evolved = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 
                                             method='srg', generator='Wegner', 
                                             lamb=lamb)
            # Load U(k, k') [unitless]
            U_matrix_unitless = SRG_unitary_transformation(H_initial,
                                                           H_evolved)

            # Isolate 2-body term and convert to fm^3
            if vnn.coupled_channel(channel):
                I_matrix_unitless = np.eye( 2*ntot, 2*ntot )
                row, col = np.meshgrid(factor_array_cc, factor_array_cc)
            else:
                I_matrix_unitless = np.eye(ntot, ntot)
                row, col = np.meshgrid(factor_array, factor_array)
            delta_U_matrix_unitless = U_matrix_unitless - I_matrix_unitless
            delta_U_matrix = delta_U_matrix_unitless / row / col
            
            # 2J+1 factor
            J = int( channel[-1] )
            
            # Add to the pp and pn terms
            # Coupled-channel
            if vnn.coupled_channel(channel):
                    
                # First L of coupled-channel
                # Isospin CG's=1/\sqrt(2) for pn
                deltaU_pn += (2*J+1)/2 * delta_U_matrix[:ntot, :ntot]
                deltaU2_pn += (2*J+1)/4 * ( delta_U_matrix[:ntot, :ntot]**2 + \
                                            delta_U_matrix[:ntot, ntot:]**2 )

                # Isospin CG's=1 for pp
                if channel in pp_channels:
                    deltaU_pp += (2*J+1) * delta_U_matrix[:ntot, :ntot]
                    deltaU2_pp += (2*J+1) * ( delta_U_matrix[:ntot, :ntot]**2 \
                                            + delta_U_matrix[:ntot, ntot:]**2 )
                    
                # Decide whether to add second L based on highest allowed
                # L value (e.g., 0 + 2 <= 2 meaning we include the 3D1-3D1
                # part of the coupled 3S1-3D1 channel if we input D-waves
                # in channels)
                if self.channel_L_value(channel) + 2 <= highest_L:
                    deltaU_pn += (2*J+1)/2 * delta_U_matrix[ntot:, ntot:]
                    deltaU2_pn += (2*J+1)/4 * ( \
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
                deltaU2_pn += (2*J+1)/4 * delta_U_matrix**2
                
                # Isospin CG's=1 for pp
                if channel in pp_channels:
                    deltaU_pp += (2*J+1) * delta_U_matrix
                    deltaU2_pp += (2*J+1) * delta_U_matrix**2

        # Evaluate pp and pn < k | \delta U | k > terms
        # These are (ntot, 1) arrays
        self.deltaU_pp_k = np.diag(deltaU_pp)
        self.deltaU_pn_k = np.diag(deltaU_pn)
        
        # Interpolate pp and pn < k | \delta U | k' > and
        # < k | \delta U^{\dagger} | k' > matrix elements to evaluate at
        # |q_vec-K_vec/2| for \deltaU \deltaU^{\dagger} term
        # These are functions of two momentum variables
        self.deltaU2_pp_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pp)
        self.deltaU2_pn_func = RectBivariateSpline(k_array, k_array,
                                                   deltaU2_pn)
        
        
        # --- Set up grids for integration over K and x --- #
        
        # Make these grids of k, K, and x for evaluating \theta functions
        self.k_grid_2d, self.x_grid_2d = np.meshgrid(k_array, x_array,
                                                     indexing='ij')
        self.k_grid_3d, self.K_grid_3d, self.x_grid_3d = np.meshgrid(k_array,
                                               K_array, x_array, indexing='ij')
        
        # The next few are grids for attaching integration weights
        _, self.x_weights_2d = np.meshgrid(k_array, x_weights, indexing='ij')
        _, _, self.x_weights_3d = np.meshgrid(k_array, K_array, x_weights,
                                              indexing='ij')
        _, self.K_weights_2d = np.meshgrid(k_array, K_weights, indexing='ij')

        
    def channel_L_value(self, channel):
        """
        Returns the L value associated with the given partial wave channel.

        Parameters
        ----------
        channel : str
            The partial wave channel (e.g. '1S0').

        Returns
        -------
        L : int
            Total orbital angular momentum associated with partial wave
            channel [unitless].

        """
        
        # This gives 'S', 'P', etc.
        channel_letter = channel[1]
        
        if channel_letter == 'S':
            return 0
        elif channel_letter == 'P':
            return 1
        elif channel_letter == 'D':
            return 2
        elif channel_letter == 'F':
            return 3
        elif channel_letter == 'G':
            return 4
        elif channel_letter == 'H':
            return 5
        else:
            print('Input channel is outside the range of this function.')
            return None
        
        
    # Grids for \theta functions and angle-averaging
    # Return n(q, kF_1, kF_2) total or 1, \delta U, \delta U\delta U^{\dagger}
    # or pp, nn, pn, and np contributions
    
    
    def theta_q_2k(self, kF, q):
        """
        Evaluates \theta( k_F - \abs(q_vec - 2k_vec) ) for every momentum k
        and angle x where \abs(q_vec - 2k_vec) = \sqrt(q^2 + 4k^2 + 4kqx).

        Parameters
        ----------
        kF : float
            Fermi momentum [fm^-1].
        q : float
            Single-nucleon momentum [fm^-1].

        Returns
        -------
        output : 2-D ndarray
            \theta function evaluated at each point in k and x. This is a
            (ntot, xtot), unitless array.

        """
    
        q_2k_magnitude = q**2 + 4 * self.k_grid_2d**2 - \
                         4 * q * self.k_grid_2d * self.x_grid_2d
        
        # This returns a (ntot, xtot) array of boolean values at every point
        # converted to 1's or 0's by multiplying 1
        theta_q_2k = ( q_2k_magnitude < kF ) * 1

        # Return the weighted matrix for integration purposes
        return theta_q_2k * self.x_weights_2d
    
    
    def theta_K_k(self, kF, sign=1):
        """
        Evaluates \theta( k_F - \abs(K_vec/2 + sign*k_vec) ) for every
        momentum k and angle x where \abs(K_vec/2 + sign * k_vec) = 
        \sqrt(K^2/4 + k^2 + sign*Kkx).

        Parameters
        ----------
        kF : float
            Fermi momentum [fm^-1].
        sign : int, optional
            Sign of k_vec in the magnitude. Only other option is -1.

        Returns
        -------
        output : 3-D ndarray
            \theta function evaluated at each point in k, K, and x. This is a
            (ntot, Ntot, xtot), unitless array.

        """
    
        K_k_magnitude = self.K_grid_3d**2/4 + self.k_grid_3d**2 + \
                        sign * self.k_grid_3d * self.K_grid_3d * self.x_grid_3d
        
        # This returns a (ntot, Ntot, xtot) array of boolean values at every
        # point converted to 1's or 0's by multiplying 1
        theta_K_k = ( K_k_magnitude < kF ) * 1

        # Return the weighted matrix for integration purposes
        return theta_K_k * self.x_weights_3d
    

    def n_lambda(self, q, kF_1, kF_2):
        """
        """
        
        return None
    
    
    
    
if __name__ == '__main__':
    
    # --- Testing --- #
    kvnn = 6
    kmax = 10.0
    kmid = 2.0
    ntot = 120
    lamb = 1.35
    channels = ['1S0', '3S1', '3P0', '3P1', '3P2', '1P1']
    
    test = single_nucleon_momentum_distributions(kvnn, channels, lamb, kmax,
                                                 kmid, ntot)
    a = test.theta_q_2k(1.35, 1.0)
    b = test.theta_K_k(1.35)
    print(a)
    print(a[0])
    print(len(b), len(b[0]), len(b[0][0]) )