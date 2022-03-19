#!/usr/bin/env python3

"""
File: srg_transformation.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 3, 2019

Function that returns an SRG unitary transformation given an initial and
SRG-evolved Hamiltonian.

Last update: March 17, 2022

"""

# To-do: Update all instances of SRG_unitary_transformation calls

# Python imports
import numpy as np
import numpy.linalg as la


def get_transformation(H_initial, H_evolved):
    """
    SRG unitary transformation built out of eigenvectors of the initial and 
    evolved Hamiltonians.
    
    Parameters
    ----------
    H_initial : 2-D ndarray
        Initial Hamiltonian matrix [MeV].
    H_evolved : 2-D ndarray
        Evolved Hamiltonian matrix [MeV].
        
    Returns
    -------
    U_matrix : 2-D ndarray
        SRG unitary transformation matrix.
        
    """
    
    # Length of the matrices
    ntot = len(H_initial)

    # Diagonalize the initial Hamiltonian
    _, vecs_initial = la.eigh(H_initial)
    # Diagonalize the evolved Hamiltonian
    _, vecs_evolved = la.eigh(H_evolved)

    # Initialize unitary transformation U with same dimension as Hamiltonians
    U_matrix = np.zeros( (ntot, ntot) )
    
    # Transformation is given by summing over the outer product of evolved and
    # initial eigenvectors
    for alpha in range(ntot):
        
        # Eigenvectors (these are already sorted correctly from eigh)
        psi_alpha_initial = vecs_initial[:, alpha]
        psi_alpha_evolved = vecs_evolved[:, alpha]
        
        # Make sure the phases match
        if psi_alpha_initial.T @ psi_alpha_evolved < 0:
            psi_alpha_evolved = -psi_alpha_evolved
        
        # Outer product of eigenvectors
        U_matrix += np.outer(psi_alpha_evolved, psi_alpha_initial)
        
    return U_matrix