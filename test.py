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
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py in Old_codes.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Test convert_ticks_to_labels function here


import numpy as np
from Figures import figures_functions as ff

a = 5

x_max = 4.0
x_stepsize = 1.0 # Step-size in labeling tick marks
x_ticks = np.arange(0.0, x_max + x_stepsize, x_stepsize)
x_ticks_strings = ff.convert_ticks_to_labels(x_ticks)
print(x_ticks_strings)

ylim = (-4.5, 2.5)
y_stepsize = 1.0
y_ticks = np.arange(ylim[0] + 0.5, ylim[1] + 0.5, y_stepsize)
y_ticks_strings = ff.convert_ticks_to_labels(y_ticks)
print(y_ticks_strings)

ylim = (-4.0, 2.0)
y_stepsize = 1.0
y_ticks = np.arange(ylim[0], ylim[1] + y_stepsize, y_stepsize)
y_ticks_strings = ff.convert_ticks_to_labels(y_ticks)
print(y_ticks_strings)

colorbar_limits = (-1.0, 1.0)
mn = colorbar_limits[0]
mx = colorbar_limits[1]
levels_ticks = np.linspace(mn, mx, 9)
levels_ticks_strings = ff.convert_ticks_to_labels(levels_ticks)
print(levels_ticks_strings)

colorbar_limits = (-0.4, 0.4)
mn = colorbar_limits[0]
mx = colorbar_limits[1]
levels_ticks = np.linspace(mn, mx, 9)
levels_ticks_strings = ff.convert_ticks_to_labels(levels_ticks)
print(levels_ticks_strings)