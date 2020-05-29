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
#   Testing mesh-dependence in V, ataq, and r^2


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from Figures import figures_functions as ff
#import observables as ob
import operators as op
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# --- Set up --- #

# Potential and channel
kvnn = 10
channel = '3S1'

# SRG details
generator = 'Wegner'
lamb = 1.5

# Mesh details
ntot = 120
mesh_1 = (10.0, 2.0, ntot)
mesh_2 = (30.0, 4.0, ntot)
meshes = [mesh_1, mesh_2]

# Operator details
q = 3.0
r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

# Initialize dictionary to store operators
d = {}

# Store unitary transformations in dictionary
d['U'] = {}
for mesh in meshes:
    
    H_initial = lsp.load_hamiltonian(kvnn, channel, kmax=mesh[0],
                                     kmid=mesh[1], ntot=ntot)
    H_evolved = lsp.load_hamiltonian(kvnn, channel, kmax=mesh[0],
                                     kmid=mesh[1], ntot=ntot, method='srg',
                                     generator=generator, lamb=1.5)
    d['U'][mesh] = SRG_unitary_transformation(H_initial, H_evolved)


# --- Potentials --- #

d['V'] = {}

# Load potentials
d['V'][kvnn] = {}
d['V'][901] = {}
for mesh in meshes:
    d['V'][kvnn][mesh] = lsp.load_potential(kvnn, channel, kmax=mesh[0],
                                            kmid=mesh[1], ntot=ntot)
    d['V'][901][mesh] = lsp.load_potential(901, channel, kmax=mesh[0],
                                           kmid=mesh[1], ntot=ntot,
                                           method='srg', generator=generator,
                                           lamb=1.2)
    
# Plot slices of initial EM N3LO potentials
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))
for i, mesh in enumerate(meshes):
    
    k_array, _ = lsp.load_momentum(kvnn, channel, kmax=mesh[0],
                                   kmid=mesh[1], ntot=ntot)
    V_diag = np.diag( d['V'][kvnn][mesh] )[:ntot]
    if i == 0:
        curve_label = 'Momentum mesh 1'
    else:
        curve_label = 'Momentum mesh 2'
    ax.plot(k_array, V_diag, color=ff.xkcd_colors(i), label=curve_label,
            linestyle=ff.line_styles(i), linewidth=2.0)

# Specify axes limits
ax.set_xlim( [0.0, 4.0] )
ax.set_ylim( [-2.0, 1.0] )
# Set x and y labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=16)
ax.set_ylabel('V(k,k) [fm]', fontsize=18)
# Set legend
ax.legend(loc='upper right', frameon=False, fontsize=14)
# Add kvnn label
kvnn_label = ff.kvnn_label_conversion(kvnn)
anchored_text = AnchoredText(kvnn_label, loc='lower right', frameon=False, 
                             prop=dict(size=18))
ax.add_artist(anchored_text)
# Save figure
f.savefig('V_EM_N3LO_test.pdf', bbox_inches='tight')

# Plot slices of evolved \Lambda=9 fm^-1 potentials
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))
for i, mesh in enumerate(meshes):
    
    k_array, _ = lsp.load_momentum(901, channel, kmax=mesh[0],
                                   kmid=mesh[1], ntot=ntot)
    V_diag = np.diag( d['V'][901][mesh] )[:ntot]
    if i == 0:
        curve_label = 'Momentum mesh 1'
    else:
        curve_label = 'Momentum mesh 2'
    ax.plot(k_array, V_diag, color=ff.xkcd_colors(i), label=curve_label,
            linestyle=ff.line_styles(i), linewidth=2.0)

# Specify axes limits
ax.set_xlim( [0.0, 4.0] )
ax.set_ylim( [-4.0, 1.5] )
# Set x and y labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=16)
ax.set_ylabel('V(k,k) [fm]', fontsize=18)
# Set legend
ax.legend(loc='upper right', frameon=False, fontsize=14)
# Add kvnn label
kvnn_label = ff.kvnn_label_conversion(901)
anchored_text = AnchoredText(kvnn_label, loc='lower right', frameon=False,
                             prop=dict(size=18))
ax.add_artist(anchored_text)
# Save figure
f.savefig('V_L9_test.pdf', bbox_inches='tight')


# --- a^{\dagger}_q a_q --- #

d['ataq'] = {}

# Load ataq operators
for mesh in meshes:
    
    k_array, k_weights = lsp.load_momentum(kvnn, channel, kmax=mesh[0],
                                           kmid=mesh[1], ntot=ntot)
    factor_array = k_array * np.sqrt(k_weights) * np.sqrt(2/np.pi)
    factor_array = np.concatenate( (factor_array, factor_array) )
    row, col = np.meshgrid(factor_array, factor_array)
    d['ataq'][mesh] = op.momentum_projection_operator(q, k_array, k_weights,
                      channel, U=d['U'][mesh]) / row / col
                
# Plot slices of evolved ataq operators
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))
for i, mesh in enumerate(meshes):
    
    k_array, _ = lsp.load_momentum(kvnn, channel, kmax=mesh[0],
                                   kmid=mesh[1], ntot=ntot)
    ataq_diag = np.diag( d['ataq'][mesh] )[:ntot]
    if i == 0:
        curve_label = 'Momentum mesh 1'
    else:
        curve_label = 'Momentum mesh 2'
    ax.plot(k_array, ataq_diag, color=ff.xkcd_colors(i), label=curve_label,
            linestyle=ff.line_styles(i), linewidth=2.0)

# Specify axes limits
ax.set_xlim( [0.0, 4.0] )
ax.set_ylim( [-0.003, 0.012] )
# Set x and y labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=16)
ax.set_ylabel(r'$a^{\dagger}_q a_q$' + '(k,k) [fm' + r'$^6$' + ']', fontsize=18)
# Set legend
ax.legend(loc='upper left', frameon=False, fontsize=11)
# Add kvnn label
kvnn_label = ff.kvnn_label_conversion(kvnn)
anchored_text = AnchoredText(kvnn_label, loc='lower left', frameon=False,
                             prop=dict(size=18))
ax.add_artist(anchored_text)
# Save figure
f.savefig('ataq_test.pdf', bbox_inches='tight')


# --- r^2 --- #

d['r2'] = {}

# Load ataq operators
for mesh in meshes:
    
    k_array, k_weights = lsp.load_momentum(kvnn, channel, kmax=mesh[0],
                                           kmid=mesh[1], ntot=ntot)
    factor_array = k_array * np.sqrt(k_weights) * np.sqrt(2/np.pi)
    factor_array = np.concatenate( (factor_array, factor_array) )
    row, col = np.meshgrid(factor_array, factor_array)
    initial_operator = op.r2_operator(k_array, k_weights, r_array, dr)
    d['r2'][mesh] = ( op.r2_operator(k_array, k_weights, r_array, dr,
                      U=d['U'][mesh]) - initial_operator ) / row / col

# Plot slices of r^2 for kvnn=10 - kmid dependence
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))
for i, mesh in enumerate(meshes):
    
    k_array, _ = lsp.load_momentum(kvnn, channel, kmax=mesh[0],
                                   kmid=mesh[1], ntot=ntot)
    r2_diag = np.diag( d['r2'][mesh] )[:ntot]
    if i == 0:
        curve_label = 'Momentum mesh 1'
    else:
        curve_label = 'Momentum mesh 2'
    ax.plot(k_array, r2_diag, color=ff.xkcd_colors(i), label=curve_label,
            linestyle=ff.line_styles(i), linewidth=2.0)

# Specify axes limits
ax.set_xlim( [0.0, 0.4] )
ax.set_ylim( [-1000, 3000] )
# Set x and y labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=16)
ax.set_ylabel(r'$\Delta r^2_{\lambda}$' + '(k,k) [fm' + r'$^5$' + ']', fontsize=18)
# Set legend
ax.legend(loc='upper right', frameon=False, fontsize=14)
# Add kvnn label
kvnn_label = ff.kvnn_label_conversion(kvnn)
anchored_text = AnchoredText(kvnn_label, loc='center right', frameon=False,
                             prop=dict(size=18))
ax.add_artist(anchored_text)
# Save figure
f.savefig('r2_test.pdf', bbox_inches='tight')


# Show r_max dependence
r_max_alt = 25.0
r_array_alt = np.arange(r_min, r_max_alt + dr, dr)
k_array, k_weights = lsp.load_momentum(kvnn, channel, kmax=mesh_1[0], 
                                       kmid=mesh_1[1], ntot=ntot)
factor_array = k_array * np.sqrt(k_weights) * np.sqrt(2/np.pi)
factor_array = np.concatenate( (factor_array, factor_array) )
row, col = np.meshgrid(factor_array, factor_array)
initial_operator = op.r2_operator(k_array, k_weights, r_array_alt, dr)
d['r2']['alt'] = ( op.r2_operator(k_array, k_weights, r_array_alt, dr,
                   U=d['U'][mesh_1]) - initial_operator ) / row / col

# Plot slices of evolved r^2 for kvnn=10 - IR cutoff
plt.close('all')
f, ax = plt.subplots(figsize=(4,4))
k_array, _ = lsp.load_momentum(kvnn, channel, kmax=mesh_1[0], kmid=mesh_1[1], ntot=ntot)
r2_diag = np.diag( d['r2'][mesh_1] )[:ntot]
r2_diag_alt = np.diag( d['r2']['alt'] )[:ntot]
ax.plot(k_array, r2_diag, color=ff.xkcd_colors(0), 
        label=r'$r_{max}=%.1f$ fm'%r_max, linestyle=ff.line_styles(0), linewidth=2.0)
ax.plot(k_array, r2_diag_alt, color=ff.xkcd_colors(1),
        label=r'$r_{max}=%.1f$ fm'%r_max_alt, linestyle=ff.line_styles(1), linewidth=2.0)

# Specify axes limits
ax.set_xlim( [0.0, 0.2] )
ax.set_ylim( [-1000, 3000] )
# Set x and y labels
ax.set_xlabel('k [fm' + r'$^{-1}$' + ']', fontsize=16)
ax.set_ylabel(r'$\Delta r^2_{\lambda}$' + '(k,k) [fm' + r'$^5$' + ']', fontsize=18)
# Set legend
ax.legend(loc='upper right', frameon=False, fontsize=14)
# Add kvnn label
kvnn_label = ff.kvnn_label_conversion(kvnn)
anchored_text = AnchoredText(kvnn_label, loc='center right', frameon=False,
                             prop=dict(size=18))
ax.add_artist(anchored_text)
# Save figure
f.savefig('r2_r_max_test.pdf', bbox_inches='tight')