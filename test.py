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
#   Looking at mesh-dependence in spherical bessel functions.


#from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from Figures import figures_functions as ff
#import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
#from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Set up --- #

# Specify potential and SRG-evolution here
kvnn = 10
channel = '3S1'
# Test two different meshes
kmax_1 = 10.0
kmid_1 = 2.0
kmax_2 = 30.0
kmid_2 = 4.0
ntot = 120

# Legend label
legend_label = r'$k_{\rm max}=%.1f$' + ' fm' + r'$^{-1}$'

# Load momentum and weights, and define factor_array
k_array_1, k_weights_1 = lsp.load_momentum(kvnn, channel, kmax=kmax_1, 
                                           kmid=kmid_1, ntot=ntot)
factor_array_1 = np.sqrt(2/np.pi) * k_array_1 * np.sqrt(k_weights_1)
row_1, col_1 = np.meshgrid(factor_array_1, factor_array_1)
k_array_2, k_weights_2 = lsp.load_momentum(kvnn, channel, kmax=kmax_2, 
                                           kmid=kmid_2, ntot=ntot)
factor_array_2 = np.sqrt(2/np.pi) * k_array_2 * np.sqrt(k_weights_2)
row_2, col_2 = np.meshgrid(factor_array_2, factor_array_2)

# Specify coordinates array
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

hank_trans_1 = op.hankel_transformation('3S1', k_array_1, r_array, dr)
hank_trans_2 = op.hankel_transformation('3S1', k_array_2, r_array, dr)

# Fix r value
#r = 0.0
#r = 2.0
r = 10.0
r_index = op.find_q_index(r, r_array)

# Take off diagonal elements
hank_trans_od_1 = hank_trans_1[:,r_index]
hank_trans_od_2 = hank_trans_2[:,r_index]

# Divide by factor_array_1? -> THIS IS NOT MESH-INDEPENDENT
#opt = True
opt = False
if opt:
    hank_trans_od_1 /= factor_array_1
    hank_trans_od_2 /= factor_array_2


# Plot momentum distributions
plt.close('all')
f, ax = plt.subplots(figsize=(4, 4))
    
# This will raise an error if you aren't plotting with respect to k
ax.plot(k_array_1, hank_trans_od_1, color=ff.xkcd_colors(0),
        label=legend_label % kmax_1, linestyle='solid', linewidth=2.0)
ax.plot(k_array_2, hank_trans_od_2, color=ff.xkcd_colors(1),
        label=legend_label % kmax_2, linestyle='dashdot', linewidth=2.0)

# Specify axes limits
ax.set_xlim([0.0, 4.0])
if opt:
    ax.set_ylim([-5, 5])
else:
    ax.set_ylim([-0.2, 0.8])
ax.legend(loc='upper right', fontsize=13, frameon=False)

# Set axes labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=18)
ax.set_ylabel(r'$\sqrt{dr} r j_0(kr)$', fontsize=20)

anchored_text = AnchoredText(r'$r=%.2f$' % r_array[r_index] + ' fm',
                             loc='lower right', prop=dict(size=14), 
                             frameon=False)
ax.add_artist(anchored_text)

# Enlarge axes tick marks
ax.tick_params(labelsize=14)

plt.show()