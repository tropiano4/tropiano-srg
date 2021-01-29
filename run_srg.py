#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: run_srg.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# SRG evolves a given potential. The run_srg function returns a dictionary of
# the evolved Hamiltonian (in scattering units [fm^-2]) at each point in
# lambda (which serves as the key). Includes an option to save the evolved
# potentials.
#
# Revision history:
#   05/28/19 --- Updated along with various updates to SRG codes. Split 
#                original script evolve_hamiltonian.py into run_srg.py and 
#                run_magnus.py.
#
#------------------------------------------------------------------------------


import numpy as np
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import vnn
from SRG import srg_wegner
from SRG import srg_kinetic_energy
from SRG import srg_block_diagonal


def run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
            lambda_bd=0.00, save=True):
    """
    SRG evolves a specified Hamiltonian to several values of lambda [fm^-1]
    and has the option to save the evolved potentials.
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0').
    kmax : float
        Maximum value in the momentum mesh [fm^-1].
    kmid : float
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int
        Number of momentum points in mesh.
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lambda_array : 1-D ndarray
        Lambda evolution values [fm^-1].
    lambda_bd : float, optional
        Cutoff for block-diagonal decoupling [fm^-1].
    save : bool, optional
        If true, saves the evolved potentials within Potentials/vsrg_macos.
        
    Returns
    -------
    d : dict
        Dictionary storing each evolved Hamiltonian with keys (floats)
        corresponding to each lambda value (e.g. d[1.5] returns the evolved
        Hamiltonian at lambda = 1.5 fm^-1).
        
    """

    # Load initial Hamiltonian, kinetic energy, momentum, and weights
    H_initial = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    T_rel = vnn.load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)     
    k_array, k_weights = vnn.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # h-bar^2 / M [MeV fm^2] for conversion from MeV to scattering units
    hbar_sq_over_m = 41.47
    
    # Initial value of lambda in units fm^-1
    # Technically, this value should be infinity, but we can take 20 fm^-1
    # which is reasonably large
    lambda_initial = 20.0

    # Initialize SRG class
    if generator == 'Wegner':
        
        evolve = srg_wegner.SRG(H_initial)
      
    elif generator == 'T':

        evolve = srg_kinetic_energy.SRG(H_initial, T_rel)
        
    elif generator == 'Block-diag':
        
        # Whether Hamiltonian is coupled channel or not (cc = True or False)
        cc = vnn.coupled_channel(channel) 
        evolve = srg_block_diagonal.SRG(H_initial, lambda_bd, k_array, cc)
        
    # Time the evolution and return dictionary d of evolved Hamiltonians where
    # the keys are lambda values
    t0 = time.time() # Start time
    d = evolve.evolve_hamiltonian(lambda_initial, lambda_array)
    t1 = time.time() # End time
    
    # Print details
    mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
    print('_'*85)
    print( 'H(s) done evolving to final lambda = %.2f fm^-1 after %f minutes'
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
            V_evolved = H_evolved - T_rel / hbar_sq_over_m
            # Save evolved potential
            vnn.save_potential(k_array, k_weights, V_evolved, kvnn, channel, 
                              kmax, kmid, ntot, 'srg', generator, lamb, 
                              lambda_bd)
                
    # Otherwise, only return the dictionary d
    return d


# Make run_srg.py executable
if __name__ == '__main__':
    
    
    # Specifications
    
    # Potentials
    kvnns = (6, 79, 111, 222)
    
    # Channels
    channels = ('3S1', '1S0')
    
    # Momentum mesh
    kmax = 10.0
    kmid = 2.0
    ntot = 120
    
    # Generators
    generators = ('Wegner', 'Block-diag')
    
    # \lambda or \Lambda_BD value in fm^-1
    lamb = 1.35 # Approximately the Fermi momentum

    # Loop over each evolution
    for kvnn in kvnns:
        for channel in channels:
            for generator in generators:
                if generator == 'Block-diag':
                    d = run_srg(kvnn, channel, kmax, kmid, ntot, generator,
                                np.array( [1.0] ), lambda_bd=lamb, save=True)
                else:
                    d = run_srg(kvnn, channel, kmax, kmid, ntot, generator,
                                np.array( [lamb] ), save=True)