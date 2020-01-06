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
#   08/28/19 --- Testing the deuteron momentum distribution using projection
#                operators from operators.py
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   09/10/19 --- Using this script to run SRG evolution on several potentials.
#   09/24/19 --- Comparing wave functions from different SRG-evolved potentials
#                by looking at momentum distribution functions.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   10/22/19 --- Comparing the initial and SRG block-diagonal evolved deuteron
#                wave functions squared. Possible connection between V_low-k
#                and block-diagonal SRG.
#   10/29/19 --- Testing r^2 operator and RMS half-radius of deuteron.
#   01/03/20 --- Looking at block-diagonal unitary transformations.
#
#------------------------------------------------------------------------------


from matplotlib.offsetbox import AnchoredText
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp
#import operators as op
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


kvnn = 111
#channel = '1P1'
#channel = '3S1'
channel = '1S0'
kmax = 8.0
kmid = 2.0
ntot = 120

k_array, _ = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)

generator = 'Block-diag'
lamb = 1.0
lambda_bd = 2.0

H_evolved = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, 'srg', 
                                generator, lamb, lambda_bd)
U_matrix = SRG_unitary_transformation(H_initial, H_evolved)

if channel == '1P1' or channel == '1S0':
    I = np.eye(ntot, ntot)
elif channel == '3S1':
    I = np.eye(2*ntot, 2*ntot)

delta_U = U_matrix - I

x_label = "k' [fm" + r'$^{-1}$' + ']'
y_label = 'k [fm' + r'$^{-1}$' + ']'
mx = 1e-2
mn = -1e-2
color_style = 'jet'

f, ax = plt.subplots(1, 1, figsize=(4, 3.5))
if channel == '1P1' or channel == '1S0':
    c = ax.pcolormesh(k_array, k_array, delta_U, cmap=color_style, vmin=mn, 
                      vmax=mx)
elif channel == '3S1':
    c = ax.pcolormesh(k_array, k_array, delta_U[:ntot,:ntot], cmap=color_style, 
                      vmin=mn, vmax=mx)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.tick_params(labeltop=True)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
channel_anchored_text = AnchoredText(ff.channel_label_conversion(channel),
                                     loc='lower right')
ax.add_artist(channel_anchored_text)
plt.gca().invert_yaxis()

f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
cbar = f.colorbar(c, cax=cbar_ax)
# Set colorbar label
cbar.ax.set_title(r"$\delta U(k,k')$")

plt.show()

# Do log-scale

mxx = 1e-2
mnn = 1e-5
color_style_log = 'Blues'

f, ax = plt.subplots(1, 1, figsize=(4, 3.5))
if channel == '1P1' or channel == '1S0':
    c = ax.pcolormesh(k_array, k_array, abs(delta_U), cmap=color_style_log, 
                      norm=colors.LogNorm(vmin=mnn, vmax=mxx))
elif channel == '3S1':
    c = ax.pcolormesh(k_array, k_array, abs(delta_U[:ntot,:ntot]), cmap=color_style_log, 
                      norm=colors.LogNorm(vmin=mnn, vmax=mxx))
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.tick_params(labeltop=True)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
channel_anchored_text = AnchoredText(ff.channel_label_conversion(channel),
                                     loc='lower right')
ax.add_artist(channel_anchored_text)
plt.gca().invert_yaxis()

f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
cbar = f.colorbar(c, cax=cbar_ax)
# Set colorbar label
cbar.ax.set_title(r"$\delta U(k,k')$")

plt.show()