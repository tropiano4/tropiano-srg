# Created 05/03/19 by A.T. (tropiano.4@osu.edu)

# SRG unitary transformation function: this function returns an SRG unitary
# transformation given an initial and SRG evolved Hamiltonian in units MeV


import numpy as np
import numpy.linalg as la


def load_omega_matrix(kvnn, channel, kmax, kmid, ntot, generator, lamb, 
                      lambda_bd=0.00, k_magnus=6, ds=1e-5):
    '''Description.'''
    
    return None


def Magnus_unitary_transformation(kvnn, channel, kmax, kmid, ntot, generator, \
                                  lamb, lambda_bd=0.00, k_magnus=6, ds=1e-5):
    '''Returns a Magnus unitary transformation by loading an evolved Omega
    matrix and taking U(s) = exp^Omega(s).'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    # generator (string): SRG generator ('Wegner', 'T', 'Block-diag')
    # lamb (float): Evolution parameter lambda in units fm^-1
    # lambda_bd (float): Lambda value for block-diagonal decoupling (e.g. 2.00 
    # fm^-1)
    # k_magnus (integer): ...
    # ds (float): ...
    
    # Load Omega(s)
    Os_matrix = 0.0
    
    # Exponentiate matrix for unitary transformation
    U = expm(Os_matrix)
    
    return U