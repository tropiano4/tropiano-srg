# Scratch Work

# Test codes in progress in this script.
# Last thing tested: Testing SciPy's ode solver for srg_kinetic_enegy.py and
# srg_wegner.py


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from Potentials.vsrg_macos import load_save_potentials as lp
  

kvnn = 901
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120

#gen = 'T'
gen = 'Wegner'

k_array = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)[0]
    
# Lambda values
lambda_list = [2.8,2.0,1.2]

# Limit of momentum in fm^-1
k_max = 3.0
    
# Store evolved potentials in dictionary
d = {}
    
for lamb in lambda_list:

    # Load S-S block
    V = lp.load_potential(kvnn, channel, kmax, kmid, ntot, method='srg', 
                          generator=gen, lamb=lamb)[:120,:120]
    
    d[lamb] = V

# Potential label
p_lbl = r'$\Lambda=9 \/ \/ fm^{-1}$'           
            
# Generator label
#G_lbl = r'$G=T_{rel}$'
G_lbl = r'$G=H_{D}$'
    
# Limits of color bar in fm
mx = 1.0
mn = -1.0
    
# Plot 1x3 figure     
plt.close('all')
    
f, (ax1,ax2,ax3)= plt.subplots(1,3,sharex=True,sharey=True,figsize=(12,4))
    
# lambda = 2.8 fm^-1
ax1.pcolormesh(k_array,k_array,d[lambda_list[0]],vmin=mn,vmax=mx,cmap='jet')
ax1.set_xlim([0,k_max])
ax1.set_ylim([0,k_max])
ax1.set_xlabel(r"$k \/ \/ [fm^{-1}]$",fontsize=16)
ax1.set_ylabel(r"$k' \/ \/ [fm^{-1}]$",fontsize=16)
ax1.xaxis.set_label_position('top')
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=True)
stepsize = 1.0
ax1.xaxis.set_ticks(np.arange(0.0,k_max+stepsize,stepsize))
ax1.yaxis.set_ticks(np.arange(0.0,k_max+stepsize,stepsize))
anchored_text_1 = AnchoredText(p_lbl,prop=dict(size=22),loc=3,frameon=False)
ax1.add_artist(anchored_text_1)
ax1.set_aspect('equal')
    
# lambda = 2.0 fm^-1
ax2.pcolormesh(k_array,k_array,d[lambda_list[1]],vmin=mn,vmax=mx,cmap='jet')
ax2.set_xlim([0,k_max])
ax2.set_ylim([0,k_max])
ax2.set_xlabel(r"$k \/ \/ [fm^{-1}]$",fontsize=16)
ax2.xaxis.set_label_position('top')
ax2.xaxis.tick_top()
ax2.tick_params(labeltop=True)
ax2.set_aspect('equal')
    
# lambda = 1.2 fm^-1
c = ax3.pcolormesh(k_array,k_array,d[lambda_list[2]],vmin=mn,vmax=mx,cmap='jet')
ax3.set_xlim([0,k_max])
ax3.set_ylim([0,k_max])
ax3.set_xlabel(r"$k \/ \/ [fm^{-1}]$",fontsize=16)
ax3.xaxis.set_label_position('top')
ax3.xaxis.tick_top()
ax3.tick_params(labeltop=True)
anchored_text_2 = AnchoredText(G_lbl,prop=dict(size=22),loc=3,frameon=False)
ax3.add_artist(anchored_text_2)
ax3.set_aspect('equal')

plt.gca().invert_yaxis()
f.subplots_adjust(hspace=0.1,wspace=0.1)
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = f.colorbar(c,cax=cbar_ax)
cbar.set_label(r'$[fm]$',rotation=0,labelpad=25,fontsize=22)

plt.show()