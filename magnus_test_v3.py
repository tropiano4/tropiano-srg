#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: magnus_test_v3.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     January 24, 2020
# 
# Test how the Magnus expansion converges looking at the norm of \eta(s) and
# \Omega(s) for various potentials.
#
#------------------------------------------------------------------------------


from os import chdir, getcwd
from math import factorial
#from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
from scipy.linalg import expm
from sympy import bernoulli
import time
# Scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp


class Magnus_test(object):
    
    
    def __init__(self, H_initial, k_magnus, ds):
        
        # Convert Hamiltonian and relative kinetic energy to scattering units
        hbar_sq_over_M = 41.47
        H_initial /= hbar_sq_over_M
        self.H_initial = H_initial
        self.N = len(H_initial)
        self.k_magnus = k_magnus
        self.ds = ds
        
        # Initialize factorial and Bernoulli number arrays for summations
        factorial_array = np.zeros(k_magnus+1)
        magnus_factors = np.zeros(k_magnus+1)
        for i in range(k_magnus+1):
            factorial_array[i] = factorial(i)
            magnus_factors[i] = bernoulli(i) / factorial(i)
        self.factorial_array = factorial_array
        self.magnus_factors = magnus_factors


    def commutator(self, A, B):

        return A @ B - B @ A

    
    def derivative(self, O_evolved):
        
        # Compute the evolving Hamiltonian with the BCH formula
        #H_evolved = self.bch_formula(self.H_initial, O_evolved, 25)
        # Use expm instead
        H_evolved = expm(O_evolved) @ self.H_initial @ expm(-O_evolved)

        # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
        # evolving Hamiltonian)
        G = np.diag( np.diag(H_evolved))
        eta = self.commutator( G, H_evolved)
        eta_norm = la.norm(eta)
        
        # Initial nested commutator ad_Omega^0(eta)
        ad = eta
        
        # k = 0 term
        dO_matrix = self.magnus_factors[0] * ad
        
        # Sum from k = 1 to k = k_magnus
        for k in range(1, self.k_magnus + 1):
            
            ad = self.commutator(O_evolved, ad)
            dO_matrix += self.magnus_factors[k] * ad
    
        return dO_matrix, eta_norm


    def euler_method(self, lamb):
        
        # Set initial s value and initial value of omega
        s = 0.0
        s_final = 1.0 / lamb**4.0
        ds = self.ds
        s_array = np.arange(0.0, s_final, ds)
        n = len(s_array)
        eta_array = np.zeros(n+1)
        omega_array = np.zeros(n+1)
        
        # Initial Magnus omega matrix
        O_evolved = np.zeros( (self.N, self.N) )

        # Start time
        print('Starting...')
        t0 = time.time() # Start time
        
        # Step in s until s_final is reached
        i = 0
        for s in s_array:

            # Next step in s
            try:
                dO, eta_norm = self.derivative(O_evolved)
                O_evolved += dO * ds
                eta_array[i] = eta_norm
                omega_array[i] = la.norm(O_evolved)
            # Check for OverflowError (this occurs when omega grows to infinity)
            except OverflowError:
                print('_'*85)
                print('Infinities or NaNs encountered in omega matrix.')
                print('Try using a smaller step-size')
                # Store in dictionary
                d = {}
                d['s_array'] = s_array[:i]
                d['eta_array'] = eta_array[:i]
                d['omega_array'] = omega_array[:i]
                return d
            except ValueError:
                print('_'*85)
                print('Infinities or NaNs encountered in omega matrix.')
                print('Try using a smaller step-size')
                # Store in dictionary
                d = {}
                d['s_array'] = s_array[:i]
                d['eta_array'] = eta_array[:i]
                d['omega_array'] = omega_array[:i]
                return d
            
            i += 1
                
        # Step to last s value
        ds_exact = s_final - s_array[-1]

        dO, eta_norm = self.derivative(O_evolved)
        O_evolved += dO * ds_exact
        eta_array[i] = eta_norm
        omega_array[i] = la.norm(O_evolved)
        s_array = np.append(s_array, s_final)
        # Evolve Hamiltonian
        H_evolved = expm(O_evolved) @ self.H_initial @ expm(-O_evolved)

        t1 = time.time() # End time
        mins = round( (t1 - t0) / 60.0, 2)
        print( 'H(s) done evolving to final lambda = %.1f fm^-1 after %f mins'
              % (lamb, mins) )
        
        # Store in dictionary
        d = {}
        d['H'] = H_evolved
        d['s_array'] = s_array
        d['eta_array'] = eta_array
        d['omega_array'] = omega_array
        
        return d


# Make executable
if __name__ == '__main__':
    
    
    # --- Set-up --- #
    
    cwd = getcwd()
    hbar_sq_over_M = 41.47
    kvnn_list = [900, 901, 902]
    #kvnn_list = [79, 111, 222]
    channel = '3S1'
    ntot = 120
    generator = 'Wegner'
    lamb = 1.5
    k_magnus = 6
    
    # Mesh and evolution specifications
    kmax = 30.0
    #kmax = 10.0
    kmid = 4.0
    #kmid = 2.0
    
    # --- Evolve and plot relevant figures --- #
    
    plt.close('all')
    row = 1
    col = 2
    f, (ax1, ax2) = plt.subplots(row, col, figsize=(4*col, 4*row))
    
    # Labels, fontsize, and location
    x_label = 's' + ' [fm' + r'$^4$' + ']'
    x_label_size = 18
    y_label_1 = r'$||\eta(s)||$' + ' [fm' + r'$^{-4}$' + ']'
    y_label_2 = r'$||\Omega(s)||$'
    y_label_size = 20
    legend_label_size = 12
    
    i = 0
    for kvnn in kvnn_list:
    
        H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)

        if kvnn == 900:
            magnus = Magnus_test(H_initial, k_magnus, 1e-5)
            curve_color = ff.xkcd_colors(0)
            curve_style = ff.line_styles(2)
        elif kvnn == 901:
            magnus = Magnus_test(H_initial, k_magnus, 1e-6)
            curve_color = ff.xkcd_colors(1)
            curve_style = ff.line_styles(1)
        elif kvnn == 902:
            magnus = Magnus_test(H_initial, k_magnus, 1e-5)
            curve_color = ff.xkcd_colors(2)
            curve_style = ff.line_styles(0)
        d = magnus.euler_method(lamb)

        kvnn_label = ff.kvnn_label_conversion(kvnn)
                      
        # Plot

        ax1.loglog(d['s_array'], d['eta_array'], color=curve_color, 
                   linestyle=curve_style, label=kvnn_label)
        ax2.loglog(d['s_array'], d['omega_array'], color=curve_color, 
                   linestyle=curve_style, label=kvnn_label)
            
        i += 1
    
    # Axes limits
    ax1.set_xlim([1e-5, 1e0])
    ax2.set_xlim([1e-5, 1e0])
    ax1.set_ylim([1e0, 1e7])
    ax2.set_ylim([1e-3, 1e4])
    # Enlarge axes tick marks
    ax1.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)
    # Axes labels
    ax1.set_xlabel(x_label, fontsize=x_label_size)
    ax1.set_ylabel(y_label_1, fontsize=y_label_size)
    ax2.set_xlabel(x_label, fontsize=x_label_size)
    ax2.set_ylabel(y_label_2, fontsize=y_label_size)
    # Add legend
    ax1.legend(loc='upper right', frameon=False, fontsize=legend_label_size)
    f.tight_layout()
    ax2.legend(loc='upper right', frameon=False, fontsize=legend_label_size)
    f.tight_layout()

    file_name = 'eta_omega_norms_kvnns_%d_%d_%d_%s_%s_k_magnus_%d_ds%.1e' % \
                 (kvnn_list[0], kvnn_list[1], kvnn_list[2], channel, 
                  generator, k_magnus, 1e-5)
    # Replace '.' with ',' in file name since LaTeX doesn't like periods
    file_name = ff.replace_periods_with_commas(file_name)
                
    # Save figure
    chdir('Figures/Operator_evolution')
    f.savefig(file_name+'.pdf', bbox_inches='tight')
    chdir(cwd)