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
#
#------------------------------------------------------------------------------


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
import operators as op
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lsp
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# Potential specifications
kvnn_list = [79, 111, 222]
channel = '3S1'
ntot = 120

generator = 'Wegner'
lamb = 1.5

r_min = 0.005
r_max = 30.2
dr = 0.005
r_array = np.arange(r_min, r_max + dr, dr)

# Labels
x_label = 'k [fm' + r'$^{-1}$' + ']'
y_label = r'$\Delta r^2_{\lambda}$' + '(k,0)'
legend_size = 11
legend_location = 'upper right'
anchored_text_label = r'$G=H_D$' + ', ' + r'$\lambda=%.1f$' % lamb + ' fm' + r'$^{-1}$'
anchored_text_size = 14
anchored_text_location = 'lower left'

x_max = 5.0
y_min = -3000
y_max = 10000
    
# Fontsize for labels and tick marks
axes_label_size = 18
    
# Size of figure
row_number = 1
col_number = 1
figure_size = (4*col_number, 4*row_number) # (width, height)
    
d = {}
# --- Load data and plot --- #
for kvnn in kvnn_list:
    
    d[kvnn] = {}
    
    H_initial = lsp.load_hamiltonian(kvnn, channel)
    k_array, k_weights = lsp.load_momentum(kvnn, channel)
    factor_array = 1
    factor_array = np.sqrt(k_weights)
    #factor_array = k_array * np.sqrt(k_weights) * 2/np.pi
    #r2_matrix_init = op.r2_operator(k_array, k_weights, r_array, dr)[:ntot, :ntot]
    H_evolved = lsp.load_hamiltonian(kvnn, channel, method='srg',
                                     generator=generator, lamb=lamb)
    # Calculate the unitary transformation
    U_matrix = SRG_unitary_transformation(H_initial, H_evolved)
    # Evolved r^2 operator (spurious factor of pi?)
    r2_matrix = op.r2_operator(k_array, k_weights, r_array, dr, U_matrix)
    
    d[kvnn] = r2_matrix[:ntot, 0] / factor_array
    
# Take differences
d['EMN-RKE'] = d[79] - d[111]
d['Gez-RKE'] = d[222] - d[111]
d['EMN-Gez'] = d[79] - d[222]
    
# Initialize figure
plt.close('all')
f, ax = plt.subplots(row_number, col_number, figsize=figure_size)
    
ax.plot(k_array, d['EMN-RKE'], label='EMN-RKE')
ax.plot(k_array, d['Gez-RKE'], label='Gez.-RKE')
ax.plot(k_array, d['EMN-Gez'], label='EMN-Gez.')
# Specify axes limits
ax.set_xlim( (0, x_max) )
ax.set_ylim( (y_min, y_max) )
# Set axes labels
ax.legend(loc=legend_location, fontsize=legend_size, frameon=False)
ax.set_xlabel(x_label, fontsize=axes_label_size)
ax.set_ylabel(y_label, fontsize=axes_label_size)
anchored_text = AnchoredText(anchored_text_label, 
                                       prop=dict(size=anchored_text_size),
                                       loc=anchored_text_location,
                                       frameon=False)
ax.add_artist(anchored_text)
plt.show()