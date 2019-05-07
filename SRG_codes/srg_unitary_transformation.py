# Created 05/03/19 by A.T. (tropiano.4@osu.edu)

# SRG unitary transformation function: this function returns an SRG unitary
# transformation given an initial and SRG evolved Hamiltonian in units MeV.


import numpy as np
import numpy.linalg as la


def SRG_unitary_transformation(H0_matrix, Hs_matrix):
    '''Returns an SRG unitary transformation built out of eigenvectors of the
    initial and evolved matrices.'''
    
    # Arguments
    
    # H0_matrix (2-D NumPy array): Initial Hamiltonian in units MeV
    # Hs_matrix (2-D NumPy array): Evolved Hamiltonian in units MeV
    
    # Dimension of the matrices
    N = len(H0_matrix)
    
    # Diagonalize unevolved Hamiltonian
    eig_unevolved, vecs_unevolved = la.eig(H0_matrix)
    # Diagonalize evolved Hamiltonian
    eig_evolved, vecs_evolved = la.eig(Hs_matrix)
    
    # Store eigenvalue and eigenvector pairs in a dictionary which will sort by 
    # lowest eigenvalue naturally
    d_unevolved = {}
    d_evolved = {}
    for i in range(N):
        d_unevolved[eig_unevolved[i]] = vecs_unevolved[:,i]
        d_evolved[eig_evolved[i]] = vecs_evolved[:,i]
        
    # Sort eigenvalues from lowest to highest energy
    eig_unevolved = np.sort(eig_unevolved)
    eig_evolved = np.sort(eig_evolved)

    # Initialize unitary transformation U
    U = np.zeros((N,N))
    # Unitary transformation is given by summing over outer product of evolved 
    # and unevolved eigenvectors
    for alpha in range(N):
        
        # Eigenvectors
        psi_alpha_unevolved = d_unevolved[eig_unevolved[alpha]]
        psi_alpha_evolved = d_evolved[eig_evolved[alpha]]
        
        # Make sure the phases match using dot-product
        if psi_alpha_unevolved.T @ psi_alpha_evolved < 0:
            psi_alpha_evolved = -psi_alpha_evolved
        
        # Outer product of eigenvectors
        U += np.outer(psi_alpha_evolved,psi_alpha_unevolved)
        
    return U