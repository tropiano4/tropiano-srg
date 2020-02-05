#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: sdt_functions.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     January 28, 2020
# 
# Functions for comparing two Hamiltonians or potentials based off spectral
# distribution theory (SDT). See https://arxiv.org/pdf/1308.5963.pdf and
# https://arxiv.org/pdf/nucl-th/0703076.pdf for more details.
#
#------------------------------------------------------------------------------


import numpy as np


def expectation_value(A):
    """
    Calculates the expectation value of A with dimensionality N. Related to
    the trace of the operator.
    
    Parameters
    ----------
    A : 2-D ndarray
        Input matrix.
    
    Returns
    -------
    output : float
        Expectation value of A.
    
    """
    
    return np.trace(A) / len(A)


def sigma(A):
    """
    Calculates the positive square root of the variance of A.
    
    Parameters
    ----------
    A : 2-D ndarray
        Input matrix.
    
    Returns
    -------
    output : float
        Square root of the variance of A.
    
    """
    
    # A^2
    A2 = A @ A
    
    return np.sqrt( expectation_value(A2) - expectation_value(A)**2 )


def inner_product(A, B):
    """
    Calculates the inner product of two matrices in SDT.
    
    Parameters
    ----------
    A : 2-D ndarray
        First input matrix.
    B : 2-D ndarray
        Second input matrix.
    
    Returns
    -------
    output : float
        Inner product of A and B.
    
    """
    
    # C is the matrix product of A^{\dagger} B (we're assuming A is real)
    C = A.T @ B
    
    return expectation_value(C) - expectation_value(A.T) * expectation_value(B)


def correlation_coefficient(A, B):
    """
    Calculates the correlation coefficient between two matrices in SDT.
    
    Parameters
    ----------
    A : 2-D ndarray
        First input matrix.
    B : 2-D ndarray
        Second input matrix.
    
    Returns
    -------
    output: float
        Correlation coefficient.
    
    """
    
    return inner_product(A, B) / ( sigma(A) * sigma(B) )