#!/usr/bin/env python3

"""
File: wave_function.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 3, 2019

Function for calculation of wave functions given a Hamiltonian.

Last update: March 17, 2022

"""

# Python imports
import numpy as np
import numpy.linalg as la

# Imports from A.T. codes
from .tools import find_index


def wave_function(
        H_matrix, eps=-2.22, U_matrix=None, print_normalization=False):
    """
    Diagonalizes the Hamiltonian and returns the wave function of the state
    nearest energy = eps. The wave function is unitless, that is, the momenta 
    and weights are factored in such that \sum_i { |psi(k_i)|^2 } = 1. For an 
    evolved wave function, enter an SRG transformation U_matrix.
    
    Parameters
    ----------
    H_matrix : 2-D ndarray
        Hamiltonian matrix [MeV].
    eps : float, optional
        Energy of the desired state [MeV]. Default is deuteron.
    U_matrix : 2-D ndarray, optional
        SRG transformation matrix [unitless]. If no transformation is provided,
        the function will not evolve the operator.
    print_normalization : bool, optional
        Option to print the normalization of the wave function as a check.
    
    Returns
    -------
    psi : 1-D ndarray
        Wave function for the specified state [unitless].
        
    Notes
    -----
    We use the following completeness relation:
        
        1 = 2/\pi \int_0^{\infty} dk k^2 |k><k|,
        
    which means the unitless wave functions include the integration measure.
        
    """

    # Diagonalize Hamiltonian
    eigenvalues, eigenvectors = la.eig(H_matrix)

    # Index of the wave function
    eps_index = find_index(eps, eigenvalues)

    # Full wave function (unitless)
    psi = eigenvectors[:, eps_index]

    # Evolve wave function?
    if U_matrix is not None:
        psi = U_matrix @ psi

    # Check normalization?
    if print_normalization:
        normalization = np.sum(psi ** 2)
        print(f"Normalization = {normalization:.5f}.")

    return psi
