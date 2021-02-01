#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a seperate file with an
# extension _testv#.py where v# corresponds to the version number. For example,
# momentum_projection_operator_testv1.py. Use the revision history below to
# document when and why these files are created.
#
# Revision history:
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py in Old_codes.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#   01/26/21 --- Renamed to test_script.py.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Calculate n_{\lambda}(q, Q=0) for different nuclei using LDA. Split into
#   pp, pn, nn, np contributions as in Ryckebusch_2019oya.


from os import getcwd, chdir
import numpy as np
# Scripts made by A.T.
from operators import find_q_index
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


def n_lambda_pair(pair, q, kvnn, lamb, k_F):
    """
    Gives the N,Np contribution to the SRG-evolved pair momentum distribution
    under the following approximations:
        1. Q_vector = 0.
        2. q >> k_F.
        3. Taking only S-waves: 1S0, 3S1-3S1, 3S1-3D1.
    
    Parameters
    ----------
    pair : str
        Type of pair: 'pp', 'pn', 'nn', or 'np'.
    q : float
        High momentum value [fm^-1].
    kvnn : int
        This number specifies the potential.
    lamb : float
        Evolution parameter lambda [fm^-1].
    k_F : float
        Fermi momentum [fm^-1].

    Returns
    -------
    XXXX FINISH THIS PART XXXXX
    
    Notes
    -----
    1. Make functions for each factor with arguments T, M_T, etc.
    2. Should lamb be replaced by k_F as an argument instead?
    3. Should it be (U_matrix - I) / row / col or U_matrix / row / col - I?
    4. Check the normalization of this function. Integrating over d3q and d3Q
       should give 1, correct?
    
    """
    
    # Set T, M_T, S, etc. based on the entered pair ('pp', 'pn', etc.)'
    if pair == 'pp' or pair == 'nn':
        
        # T = 1
        # M_T = 1 or M_T = -1 (doesn't matter!)
        # Make function that takes T, M_T as argument?
        # Clebsch-Gordan coefficients for isospin are 1 here
        # < 1/2 1/2 | 1 1 > = 1 or < -1/2 -1/2 | 1 -1 > = 1
        isospin_factor = 1
        
        # Select the allowed value of S and M_S
        # S = 0
        # M_S = 0
        # Make function that takes S, M_S as argument?
        
        # We sum over m_s, m_s', m_s_1, m_s_2 restricted to the following
        # values:
        #   m_s = up, m_s' = down, m_s_1 = up, m_s_2 = down
        #   m_s = up, m_s' = down, m_s_1 = down, m_s_2 = up
        #   m_s = down, m_s' = up, m_s_1 = up, m_s_2 = down
        #   m_s = down, m_s' = up, m_s_1 = down, m_s_2 = up
        # Each time we get factors of 1/\sqrt(2) or -1/\sqrt(2) but
        # appearing in even powers - thus, it's always a factor of 1/4.
        # But there are 4 possible combinations and we are summing over them!
        spin_factor = 1
        
        # Select the associated partial wave channel
        channel = '1S0'
        # This means L=0, M_L=0, J=0, M_J=0
        # < M_L M_S | J M_J> factors are always 1 in this case
        total_ang_momentum_factor = 1
        # Is this ever not 1?
        
        # Expansion of k_vector and q_vector factors
        k_expansion_factor = 1 / (2*np.pi**2) # ( \sqrt(2/\pi)*Y00 )^2
        q_expansion_factor = 1 / (2*np.pi**2) # ( \sqrt(2/\pi)*Y00 )^2
        
        # Let's combine all the relevant factors here to be concise
        factor = isospin_factor * spin_factor * total_ang_momentum_factor * \
                 k_expansion_factor * q_expansion_factor
        
        # Load initial and evolved Hamiltonians with momentum mesh
        k_array, k_weights = vnn.load_momentum(kvnn, channel)
        # Length of k_array
        ntot = len(k_array)
        factor_array = k_array * np.sqrt( 2*k_weights / np.pi )
        row, col = np.meshgrid(factor_array, factor_array)
        
        H_initial = vnn.load_hamiltonian(kvnn, channel) # MeV
        # Use band-diagonal evolution as default
        H_evolved = vnn.load_hamiltonian(kvnn, channel, method='srg', 
                                         generator='Wegner', lamb=lamb) # MeV
        
        # Calculate SRG transformation
        U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
        
        # Subtract out identity matrix and divide out momentum and weights
        I = np.eye(ntot, ntot)
        delta_U_matrix = (U_matrix - I) / row / col
        
        # Now find the matrix element of delta_U(k, q)
        q_index = find_q_index(q, k_array)
        delta_U_vector = delta_U_matrix[:ntot, q_index]
        # Do the same for delta_U^{\dagger}
        delta_U_dag_vector = delta_U_matrix.T[q_index, :ntot]
        
        # The overall contribution to the evolved pair momentum distribution
        # as a function of k (which is summed to k_F later)
        pair_contribution_k = delta_U_vector * delta_U_dag_vector * factor
        # 1/(2*\pi^2)^2 * \delta U_1S0(k,q) * \delta U_1S0^{\dagger}(q,k)
                              
        # Sum over k up to k_F (this should probably be an intergral?)
        # \sum_k k_i^2 w_i -> \int dk k^2
        # I have \sum_k... But the weights and k factors should have come from
        # U O U^t (two integrals over k and k')
        # One of the k' is done by \delta function
        # There is one more. So yes, do 2 / \pi \sum_k k_i^2 w_i
        integrand = factor_array * pair_contribution_k
        
        pair_contribution = 1/8 * np.sum(integrand)

    # Worry about this later   
    elif pair == 'pn' or pair == 'np':
        
        T_values = [0, 1]
        
        # Loop over T_values and calculate isospin coefficients
    
        # Sum over each spin projection
    
    
        # Select the associated partial wave channel
        
        pair_contribution_k = 0

        # should have a vector that depends on k values
        
    else:
        
        print("Enter a valid pair: 'pp', 'pn', 'nn', or 'np'.")
        return None

    return pair_contribution


# Perhaps this function should go in a script within Densities/HFBRAD_SLY4
def load_density(nucleus, nucleon, Z, N):
    """
    Loads a nucleonic density for the given nucleus.
    
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
        
    Returns
    -------
    r_array : 1-D ndarray
        Coordinates array [fm].
    rho_array : 1-D ndarray
        Nucleonic density as a function of r [# of protons or neutrons / vol].
    
    """


    cwd = getcwd()
    
    # Go to directory corresponding to specified nucleus
    densities_directory = 'Densities/HFBRAD_SLY4/%s' % nucleus
    chdir(densities_directory)

    # Load .dens file
    file_name = '%s_%d_%d.dens' % (nucleon, N, Z)
    table = np.loadtxt(file_name)
    
    chdir(cwd) # Go back to current working directory
    
    r_array = table[:, 0]
    rho_array = table[:, 1]
    
    return r_array, rho_array


# Check normalization of \int dr r^2 \rho_N(r) = ???
# 4*\pi \int dr r^2 \rho_Z(r) = Z (likewise for N)


def local_density_approximation():
    
    # k_F = something dependent on \rho(r)
    # Derive this
    # Evaluate k_F at each r_i plugging this into n_lambda_pair as well
    #   So you have an array n_lambda_pair_i(k_F_i) \rho_i
    #   Then sum the array and weight with 4 \pi r^2 dr (and divide by
    #   4 \pi \int r^2 dr \rho(r))
    
    return None
