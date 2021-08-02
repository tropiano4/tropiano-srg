#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: srg_unitary_transformation.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 3, 2019
# 
# SRG unitary transformation function returns an SRG unitary transformation 
# given an initial and SRG evolved Hamiltonian in units MeV.
#
# Revision history:
#   04/02/21 --- Switched la.eig to la.eigh for hermitian matrices after
#                issues with eig returning complex eigenvalues with 0j values.
#
#------------------------------------------------------------------------------


import numpy as np
import numpy.linalg as la


def SRG_unitary_transformation(H_initial, H_evolved):
    """
    SRG unitary transformation built out of eigenvectors of the initial and 
    evolved Hamiltonian matrices.
    
    Parameters
    ----------
    H_initial : 2-D ndarray
        Initial Hamiltonian matrix in units MeV.
    H_evolved : 2-D ndarray
        Evolved Hamiltonian matrix in units MeV.
        
    Returns
    -------
    U : 2-D ndarray
        SRG unitary transformation matrix (unitless).
        
    """
    
    # Length of the matrices
    N = len(H_initial)

    # Diagonalize the initial Hamiltonian
    _, vecs_initial = la.eigh(H_initial)
    # Diagonalize the evolved Hamiltonian
    _, vecs_evolved = la.eigh(H_evolved)

    # Initialize unitary transformation U
    U = np.zeros( (N, N) )
    
    # The unitary transformation is given by summing over the outer product of 
    # evolved and initial eigenvectors
    for alpha in range(N):
        
        # Eigenvectors (these are already sorted correctly from eigh)
        psi_alpha_initial = vecs_initial[:, alpha]
        psi_alpha_evolved = vecs_evolved[:, alpha]
        
        # Make sure the phases match using dot product
        if psi_alpha_initial.T @ psi_alpha_evolved < 0:
            psi_alpha_evolved = -psi_alpha_evolved
        
        # Outer product of eigenvectors
        U += np.outer(psi_alpha_evolved, psi_alpha_initial)
        
    return U