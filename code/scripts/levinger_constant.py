#!/usr/bin/env python3

"""
File: levinger_constant.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: July 27, 2022

Function that computes the Levinger constant L. This quantity relates nuclear
photo-absorption to deuteron photo-disintegration
    
    \sigma_A(E_\gamma) = L N*Z/A \sigma_d(E_\gamma),
    
at high energies E_\gamma. We use momentum distributions at high relative
momentum to extract the Levinger constant L.

Last update: July 27, 2022

"""

# Python imports
import numpy as np


def compute_L(A, Z, N, n_pn_array, n_d_array):
    """
    Levinger constant L evaluated by averaging over the ratio of the in-
    medium pn momentum distribution (with C.o.M. momentum Q integrated out)
    over the deuteron momentum distribution.
    
        L \approx A/(N*Z) * n_pn^A(q) / n^d(q),
        
    for q >> SRG \lambda.

    Parameters
    ----------
    A : int
        Nuclear mass number (A = Z + N).
    Z : int
        Proton number of nucleus.
    N : int
        Neutron number of nucleus.
    n_pn_array : 1-D ndarray
        In-medium pn pair momentum distribution [fm^3] with respect to a
        range of high relative momentum values, where the C.o.M. momentum Q
        is integrated out.
    n_d_array : 1-D ndarray
        Deuteron momentum distribution [fm^3] with respect to the same range
        of high relative momentum values.
    
    Returns
    -------
    output : float
        Levinger constant [unitless].
        
    """
    
    # Check that the numbers A, Z, and N make sense
    if A != Z+N:
        raise RuntimeError("Incomptatible values of A, Z, and N.")
        
    ratio_array = A/(N*Z) * n_pn_array / n_d_array
    
    return np.mean(ratio_array)