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
#   Testing an automatic way of selecting relevant quantum numbers for SRC
#    work (e.g., given m_s, m_s', m_t, m_t' what are the allowed S, T, etc.?)


import numpy as np


def total_isospin(m_t, m_tp):
    """
    Returns the possible values of T and M_T given m_t and m_tp.
    
    Parameters
    ----------
    m_t : float
        Isospin projection of particle 1. Must be 1/2 or -1/2.
    m_tp : float
        Isospin projection of particle 2. Must be 1/2 or -1/2.
        
    Returns
    -------
    T_array : 1-D ndarray
        Possible values of total isospin.
    M_T : int
        Total isospin projection.
        
    """
    
    M_T = int( m_t + m_tp )
    
    if abs(M_T) == 1:
        T_array = np.array( [1] )
    else:
        T_array = np.array( [0, 1] )
    
    return T_array, M_T


def total_spin(m_s, m_sp, T):
    """
    Returns the possible values of S and M_S given m_s, m_sp, and T.
    
    Parameters
    ----------
    m_s : float
        Spin projection of particle 1. Must be 1/2 or -1/2.
    m_sp : float
        Spin projection of particle 2. Must be 1/2 or -1/2.
    T : int
        Total isospin projection.
        
    Returns
    -------
    S : int
        Total spin.
    M_S : int
        Total spin projection.
        
    """
    
    M_S = int( m_s + m_sp )
    
    if T == 0:
        S = 1
    else:
        S = 0
        
    if abs(M_S) == 1 and T == 0:
        S = 1
    elif abs(M_S) == 1 and T == 1:
        S = None # Not allowed
    elif M_S == 0 and T == 0:
        S = 1
    else:
        S = 0
        
    return S, M_S


def orbital_angular_momentum(S):
    """
    Returns the possible orbital angular momentum values excluding 3+.

    Parameters
    ----------
    S : int
        Total spin.

    Returns
    -------
    L_array : 1-D ndarray
        Possible values of orbital angular momentum.

    """
    
    if S == 0:
        L_array = np.array( [0] ) # Only 1S0
    else:
        L_array = np.array( [0, 2] ) # S=1 -> 3S1 - 3D1 coupled channel
        
    return L_array
        

def orbital_angular_momentum_projection(L):
    """
    Returns possible orbital angular momentum projections.

    Parameters
    ----------
    L : int
        Orbital angular momentum.

    Returns
    -------
    M_L_array : 1-D ndarray
        Possible values of orbital angular momentum projection.

    """
    
    if L == 0:
        M_L_array = np.array( [0] )
    else: # L=2
        M_L_array = np.array( [-2, -1, 0, 1, 2] )
        
    return M_L_array


def total_angular_momentum(M_L, M_S):
    """
    Returns the possible values for 

    Parameters
    ----------
    M_L : TYPE
        DESCRIPTION.
    M_S : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    return None



if __name__ == '__main__':
    
    particle_spins = np.array( [-0.5, 0.5] )
    
    for m_t in particle_spins:
        for m_tp in particle_spins:
            
            print('_'*50)
            print( 'm_t = %.1f, m_tp = %.1f' % (m_t, m_tp) )
            T_array, M_T = total_isospin(m_t, m_tp)
            
            for T in T_array:
                
                print( 'T = %d, M_T = %d' % (T, M_T) )
                print('-'*50)
                
                for m_s in particle_spins:
                    for m_sp in particle_spins:
                        
                        print( 'm_s = %.1f, m_sp = %.1f' % (m_s, m_sp) )
                        S, M_S = total_spin(m_s, m_sp, T)
                        if S == 1 or S == 0:
                            print( 'S = %d, M_S = %d\n' % (S, M_S) )
                        else:
                            print('Not allowed.\n')
                            
                print('-'*50)
                