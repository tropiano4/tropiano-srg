#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: gauss_legendre_quad.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     June 15, 2020
# 
# Calculates nodes and weights for Gauss-Legendre quadrature integration. Use
# this for generating momentum space nodes and weights.
#
# Revision history:
#   xx/xx/xx --- ...
#
#------------------------------------------------------------------------------


import numpy as np
import numpy.linalg as la
from numpy.polynomial.legendre import leggauss
from observables import phase_shifts


# Define function and interval
k_min = 0.0
k_max = 10.0
ntot = 15
x_array, k_weights = leggauss(ntot)
k_array = 0.5 * (x_array + 1) * (k_max - k_min) + k_min

# Set up kinetic energy
# h-bar^2 / M [MeV fm^2]
hbar_sq_over_M = 41.47
T_matrix = hbar_sq_over_M * np.diag( k_array**2 ) # MeV

# Set up potential
lamb = 4.0 # units fm^-1
V_0 = 10.0 # units MeV fm^3
k_row, k_col = np.meshgrid(k_array, k_array)
# Units are MeV fm^3
V_matrix = V_0 * np.exp( -( k_row**2 / lamb**2 + k_col**2 / lamb**2 ) )

# Convert to MeV
factor_array = np.sqrt( (2*k_weights) / np.pi ) * k_array
row, col = np.meshgrid(factor_array, factor_array)

# Calculate Hamiltonian in MeV and eigenenergies
H_matrix_MeV = T_matrix + V_matrix * row * col
eigenvalues, eigenvectors = la.eig(H_matrix_MeV)


delta_array = phase_shifts(2 * np.sort(eigenvalues), V_matrix / hbar_sq_over_M,
                           k_array, k_weights)

for i in range(len(eigenvalues)):
    print('E=%.2f MeV, \delta=%.2f deg' % (np.sort(eigenvalues)[i], delta_array[i]))