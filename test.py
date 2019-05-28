#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes.
#
# Last thing tested:
#   Using plots to verify correct SRG evolution of various potentials with
#   updates to SRG code.
#
#------------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from Potentials.vsrg_macos import load_save_potentials as lp

kvnn = 6
#kvnn = 10
#kvnn = 105
#kvnn = 222
#kvnn = 900

channel = '3S1'

kmax = 30.0
#kmax = 10.0
#kmax = 8.0

kmid = 4.0
#kmid = 2.0

ntot = 120

#gen = 'Wegner'
gen = 'T'
#gen = 'Block-diag'
lambda_bd = 2.00
#lambda_bd = 3.00

k_array = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)[0]
    
# Lambda values
lambda_list = [3.0, 2.0, 1.5]

# Limit of momentum in fm^-1
k_max = 3.0
    
# Store evolved potentials in dictionary
d = {}
    
for lamb in lambda_list:

    # Load potential
    V_full = lp.load_potential(kvnn, channel, kmax, kmid, ntot, method='srg', 
                               generator=gen, lamb=lamb, lambda_bd=lambda_bd)
    # Resize to S-S block
    V = V_full[:120, :120]
    
    d[lamb] = V

# Potential label
p_lbl = 'kvnn = %d'%kvnn           
            
# Generator label
if gen == 'Wegner':
    G_lbl = r'$G=H_{D}$'
elif gen == 'T':
    G_lbl = r'$G=T_{rel}$'
elif gen == 'Block-diag':
    G_lbl = r'$G=H_{BD}$' + '\n' + r'$\Lambda_{BD} = %.2f \/ fm^{-1}$' \
            % lambda_bd
    
# Limits of color bar in fm
mx = 1.0
mn = -1.0
    
# Plot 1x3 figure     
plt.close('all')
    
f, (ax1, ax2, ax3)= plt.subplots(1, 3, sharex=True, sharey=True, 
                                 figsize=(12, 4) )
    
# lambda = 2.8 fm^-1
ax1.pcolormesh(k_array, k_array, d[ lambda_list[0] ], vmin=mn, vmax=mx, 
               cmap='jet')
ax1.set_xlim( (0, k_max) )
ax1.set_ylim( (0, k_max) )
ax1.set_xlabel(r"$k \/ \/ [fm^{-1}]$", fontsize=16)
ax1.set_ylabel(r"$k' \/ \/ [fm^{-1}]$", fontsize=16)
ax1.xaxis.set_label_position('top')
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=True)
stepsize = 1.0
ax1.xaxis.set_ticks( np.arange(0.0, k_max + stepsize, stepsize) )
ax1.yaxis.set_ticks( np.arange(0.0, k_max + stepsize, stepsize) )
anchored_text_1 = AnchoredText( p_lbl, prop=dict(size=22), loc=3, 
                               frameon=False )
ax1.add_artist(anchored_text_1)
ax1.set_aspect('equal')
    
# lambda = 2.0 fm^-1
ax2.pcolormesh(k_array, k_array, d[ lambda_list[1] ], vmin=mn, vmax=mx, 
               cmap='jet')
ax2.set_xlim( (0, k_max) )
ax2.set_ylim( (0, k_max) )
ax2.set_xlabel(r"$k \/ \/ [fm^{-1}]$", fontsize=16)
ax2.xaxis.set_label_position('top')
ax2.xaxis.tick_top()
ax2.tick_params(labeltop=True)
ax2.set_aspect('equal')
    
# lambda = 1.2 fm^-1
c = ax3.pcolormesh(k_array, k_array,d[ lambda_list[2] ], vmin=mn, vmax=mx, 
                   cmap='jet')
ax3.set_xlim((0, k_max))
ax3.set_ylim((0, k_max))
ax3.set_xlabel(r"$k \/ \/ [fm^{-1}]$", fontsize=16)
ax3.xaxis.set_label_position('top')
ax3.xaxis.tick_top()
ax3.tick_params(labeltop=True)
anchored_text_2 = AnchoredText( G_lbl, prop=dict(size=22), loc=3, 
                               frameon=False)
ax3.add_artist(anchored_text_2)
ax3.set_aspect('equal')

plt.gca().invert_yaxis()
f.subplots_adjust(hspace=0.1, wspace=0.1)
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes( (0.85, 0.15, 0.05, 0.7) )
cbar = f.colorbar(c, cax=cbar_ax)
cbar.set_label(r'$[fm]$', rotation=0, labelpad=25, fontsize=22)

plt.show()