#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: sdt_functions.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     January 28, 2019
# 
# Description of program.
# https://arxiv.org/pdf/1308.5963.pdf
# https://arxiv.org/pdf/nucl-th/0703076.pdf
#
#------------------------------------------------------------------------------


import numpy as np


def expectation_value(A):
    """
    Description of function.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    
    return np.trace(A) / len(A)


def sigma(A):
    """
    Description of function.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    
    A2 = A @ A
    
    return np.sqrt( expectation_value(A2) - expectation_value(A)**2 )


def inner_product(A, B):
    """
    Description of function.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    
    C = A.T @ B
    
    
    return expectation_value(C) - expectation_value(A.T) * expectation_value(B)


def correlation_coefficient(A, B):
    """
    Description of function.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    
    return inner_product(A, B) / ( sigma(A) * sigma(B) )