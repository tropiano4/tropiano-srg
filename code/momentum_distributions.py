#!/usr/bin/env python3

"""
File: name.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: Month Day, Year

1-3 summary of the script.

Last update: March 17, 2022

"""

# To-do: ...

# Python imports

# Imports from A.T. codes


# code


# Work flow for momentum distributions code:
# Don't worry about file organization yet

# [class Momentum_distributions]
# __init__ shouldn't do much.

# Get \delta U [class Momentum_distributions]
# Set-up \delta U_{\tau,\tau'}(k,k') matrix elements and interpolate.
# This relies on kvnn, channels (except for deuteron), \lambda, generator, and
# possibly additional arguments for an initial Hamiltonian that's already
# softened (same kvnn but \lambda_initial), or an initial Hamiltonian that is 
# inverse-transformed by a harder potential (hard kvnn and \delta \lambda).
# Will call Potential(kvnn, channel, kmax, kmid, ntot).
# * Should momentum distributions modules take in deltaU_matrices, deltaU_func,
#   self.deltaU_matrices, self.deltaU_func, or momentum distributions object?

# Compute I term [sub-classes of Momentum_distributions]
# snmd.py and dmd.py are similar up to a factor of 4*\pi (suggests using 
# super()). pmd.py has different arguments. Each function relies on input mesh-
# grids of dimension corresponding to the number of integrals. Calls functions
# to weight with \theta(\kF-q) or evaluate angular average of two \theta's.

# Compute \delta U term [sub-classes of Momentum_distributions]
# All scripts are fairly different in this part. Again relies on meshgrids.
# Also calls functions to evaluate angular average of \theta function.

# Compute \delta U \delta U^\dagger term [sub-classes of Momentum_distributions]
# Same as above.

# [sub-classes of Momentum_distributions] Not sure about this yet!
# There should probably be a distinction of multiplying by \theta(...) and
# multiplying by \int dx/2 \theta(...). Where should these functions be?
# 1. Separate .py file in code/modules as sub-class of momentum distributions?

# [sub-classes of Momentum_distributions]
# Do only n_contributions(meshgrids, kF_tau_array, kF_taup_array) which
# calculates the contributions to the momentum distribution (for each script:
# snmd.py, pmd.py, dmd.py).

# [sub-classes of Momentum_distributions]
# Add method that calculates n(q, Q=0) for pmd.py.

# [Main script]
# SRG currently does the calculation using run_srg() which calls
# save_potential(). In the same spirit, we could do 
# run_momentum_distribution() and save_momentum_distribution() where the
# latter (maybe) relies on np.savetxt(). Include option for n(q, Q=0).
# Will rely on densities.py. Do not break pn and np convention.
# run_momentum_distribution() will call set_matrix_elements() which saves
# \delta U matrix elements (called using super()). Alternatively, can make the
# matrix elements an argument of n_contributions().

# [Main script]
# Add load_momentum_distribution() function which gives the
# momentum distribution given distribution_type = 'pn', 'pp', 'nn', 'p', 'n',
# or 'd', and optional arguments: nucleus, Z, N, density_type,
# interpolate=False, contributions=False, zero_total_momentum=False.
# * This could extend to three functions with an overheaf function if you want.
# * output : Either q_array, (and Q_array?), n_array or n_func.

# [Main script]
# Add method to integrate out Q-dependence of pmd's. See quasideuteron.ipynb
# for example.

# Copy the last function of dmd.py and add to src.ipynb for now. Pretty sure
# it's only used in one figure suggesting that it should be part of the
# function that plots that figure.

# Same idea for partial wave decomposition (see above).

# Add a note to try and understand the K integration in deuteron (does that
# even make sense?)