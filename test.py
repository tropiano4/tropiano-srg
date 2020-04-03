#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a seperate file with an
# extension _testv#.py where v# corresponds to the version number. For example,
# momentum_projection_operator_testv1.py. Use the revision history below to
# document when and why these files are created.
#
# Revision history:
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/19 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Generating a few figures relevant to r^2 operator evolution to understand
#   how r^2 evolves for different potentials and SRG generators.


from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from Figures import figures_functions as ff
import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Set up --- #

# Specify potential and SRG-evolution here


# Plotting specifications that are dependent on the settings above


# --- Main calculations --- #

# Specify coordinates array
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)


# --- Plot momentum distribution --- #


# --- Plot r^2 operator --- #


# --- Plot matrix elements of expectation value --- #