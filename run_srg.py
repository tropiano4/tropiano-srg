#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: run_srg.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# SRG evolves a given Hamiltonian in units MeV. The run_srg function returns a 
# dictionary of the evolved Hamiltonian (in scattering units [fm^-2]) at each 
# point in lambda (which serves as the key).
#
# Revision history:
#   May 28, 2019 --- Updated along with various updates to SRG codes. Split 
#                    original script evolve_hamiltonian.py into run_srg.py and 
#                    run_magnus.py.
#
#------------------------------------------------------------------------------


import numpy as np
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
from SRG_codes import srg_wegner
from SRG_codes import srg_kinetic_energy
from SRG_codes import srg_block_diagonal


def run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
            lambda_bd=0.00, save=True):
    """
    SRG evolves a specified Hamiltonian to several values of lambda [fm^-1] and 
    has the option to save the evolved potentials.
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel ('1S0', '3S1', etc.)
    kmax : float
        Maximum value in the momentum mesh.
    kmid : float
        Mid-point value in the momentum mesh.
    ntot : int
        Number of momentum points in mesh.
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lambda_array : 1-D ndarray
        Lambda evolution values in units fm^-1.
    lambda_bd : float, optional
        Cutoff for block-diagonal decoupling.
    save : bool, optional
        If true, saves the evolved potentials within Potentials/vsrg_macos.
        
    Returns
    -------
    d : dict
        Dictionary storing each evolved Hamiltonian with keys (floats)
        corresponding to each lambda value (e.g. d[1.5] returns the evolved
        Hamiltonian at lambda = 1.5 fm^-1).
        
    """

    # Load initial Hamiltonian, kinetic energy and weights
    H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    T_rel = lp.load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)     
    k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # h-bar^2 / M [MeV fm^2] for conversion from MeV to scattering units
    hbar_sq_over_M = 41.47
    
    # Initial value of lambda in units fm^-1 (most potentials in 
    # Potentials/vsrg_macos are generated at lambda = 12 fm^-1 but check 
    # run_generate_vsrg_vlowk.pl to be sure)
    lambda_initial = 12.0

    # Initialize SRG class
    if generator == 'Wegner':
        
        evolve = srg_wegner.SRG(H_initial)
      
    elif generator == 'T':

        evolve = srg_kinetic_energy.SRG(H_initial, T_rel)
        
    elif generator == 'Block-diag':
        
        # Whether Hamiltonian is coupled channel or not (cc = True or False)
        cc = lp.coupled_channel(channel) 
        evolve = srg_block_diagonal.SRG(H_initial, lambda_bd, k_array, cc)
        
    # Time the evolution and return dictionary d of evolved Hamiltonians where
    # the keys are lambda values
    t0 = time.time() # Start time
    d = evolve.evolve_hamiltonian(lambda_initial, lambda_array)
    t1 = time.time() # End time
    
    # Print details
    mins = round((t1-t0)/60.0,2) # Minutes elapsed evolving H(s)
    print('_'*85)
    print( 'H(s) done evolving to final lambda = %.1f fm^-1 after %f minutes'
          % (lambda_array[-1], mins) )
    print('_'*85)
    print('\nSpecifications:\n')
    print( 'kvnn = %d, channel = %s' % (kvnn, channel) )
    print( 'kmax = %.1f, kmid = %.1f, ntot = %d' % (kmax, kmid, ntot) )
    print( 'method = srg, generator = %s' % generator )
    if generator == 'Block-diag':
        print( 'block-diagonal lambda = %.2f fm^-1' % lambda_bd )
    
    # Writes evolved potentials to files if save = True
    if save:

        # Save evolved potential for each lambda value
        for lamb in lambda_array:

            H_evolved = d[lamb] # Scattering units here [fm^-2]
            # Subtract off kinetic energy (need to convert T_rel from MeV to 
            # fm^-2)
            V_evolved = H_evolved - T_rel / hbar_sq_over_M
            # Save evolved potential
            lp.save_potential(k_array, k_weights, V_evolved, kvnn, channel, 
                              kmax, kmid, ntot, 'srg', generator, lamb, 
                              lambda_bd)
                
    # Otherwise, only return the dictionary d
    return d


# Make run_srg.py executable
if __name__ == '__main__':
    
    
    # Specify potential
    
    kvnn = 6 # AV18
    #kvnn = 10 # EM N3LO
    #kvnn = 105 # RKE N3LO at Lambda = 400 MeV
    #kvnn = 106 # RKE N3LO at Lambda = 450 MeV
    #kvnn = 107 # RKE N3LO at Lambda = 500 MeV
    #kvnn = 111 # RKE N4LO at Lambda = 450 MeV
    #kvnn = 112 # RKE N4LO at Lambda = 500 MeV
    #kvnn = 900 # Wendt at Lambda = 4 fm^-1
    #kvnn = 901 # Wendt at Lambda = 9 fm^-1
    #kvnn = 902 # Wendt at Lambda = 20 fm^-1
    #kvnn = 222  # Gezerlis et al N2LO local potential at R_0 = 1.0 fm cutoff
    #kvnn = 224  # Gezerlis et al N2LO local potential at R_0 = 1.2 fm cutoff
    
    channel = '3S1'
    
    kmax = 30.0
    #kmax = 8.0
    #kmax = 10.0
    
    kmid = 4.0
    #kmid = 2.0
    
    ntot = 120
    
    # Specify evolution
    
    #generator = 'Wegner'
    generator = 'T'
    #generator = 'Block-diag'
    
    #lambda_bd = 1.00 # This won't affect the band-diagonal generators
    lambda_bd = 2.00 
    #lambda_bd = 3.00
    #lambda_bd = 4.00
    
    #lambda_array = np.array( (10.0, 2.8, 2.0, 1.2) )
    lambda_array = np.array( (6.0, 3.0, 2.0, 1.5) )
    #lambda_array = np.array( (2.8) )
    #lambda_array = np.array( (1.2) )
    
    # Save the evolved Hamiltonian?
    save = True
    #save = False
    
    # Evolve Hamiltonian
    d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
                lambda_bd, save)