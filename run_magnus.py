#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: run_magnus.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     June 3, 2019
# 
# Magnus evolves a given Hamiltonian in units MeV. The run_magnus function 
# returns a dictionary of the evolved Hamiltonian (in scattering units [fm^-2])
# and omega matrix at each point in lambda (which serves as the key).
#
#------------------------------------------------------------------------------


import numpy as np
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
from Magnus_codes import magnus_wegner
from Magnus_codes import magnus_kinetic_energy


def run_magnus(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
               k_magnus=6, ds=1e-5, save=True):
    """
    Magnus evolves a specified Hamiltonian to several values of lambda [fm^-1]
    and has the option to save the evolved potentials and omega matrices.
    
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
        SRG generator 'Wegner' or 'T'.
    lambda_array : 1-D ndarray
        Lambda evolution values in units fm^-1.
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... )
    ds : float, optional
        Step-size in the flow parameter s.
    save : bool, optional
        If true, saves the evolved potentials within Potentials/vsrg_macos.
        
    Returns
    -------
    d : dict
        Dictionary storing each evolved Hamiltonian and omega matrix with keys
        (floats) corresponding to each lambda value (e.g. d['hamiltonian'][1.5]
        and d['omega'][1.5] returns the evolved Hamiltonian and omega matrix 
        at lambda = 1.5 fm^-1, respectively).
        
    """

    # Load initial Hamiltonian, kinetic energy and weights
    H_initial = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    T_rel = lp.load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)     
    k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # h-bar^2 / M [MeV fm^2] for conversion from MeV to scattering units
    hbar_sq_over_M = 41.47

    # Initialize SRG class
    if generator == 'Wegner':
        
        evolve = magnus_wegner.Magnus(H_initial, k_magnus, ds)
      
    elif generator == 'T':

        evolve = magnus_kinetic_energy.Magnus(H_initial, T_rel, k_magnus, ds)
        
    # Time the evolution and return dictionary d of evolved Hamiltonians and
    # omega matrices where the keys are lambda values
    t0 = time.time() # Start time
    # Check for OverflowError (this occurs when omega grows to infinity)
    try:
        d = evolve.evolve_hamiltonian(lambda_array)
    except OverflowError:
        print('_'*85)
        print('Infinities NaNs encountered in omega matrix.')
        print('Try using a smaller step-size or running magnus_split.py.')
        return None
    except ValueError:
        print('_'*85)
        print('Infinities NaNs encountered in omega matrix.')
        print('Try using a smaller step-size or running magnus_split.py.')
        return None
    t1 = time.time() # End time
    
    # Print details
    mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
    print('_'*85)
    print( 'H(s) done evolving to final lambda = %.1f fm^-1 after %f minutes'
          % (lambda_array[-1], mins) )
    print('_'*85)
    print('\nSpecifications:\n')
    print( 'kvnn = %d, channel = %s' % (kvnn, channel) )
    print( 'kmax = %.1f, kmid = %.1f, ntot = %d' % (kmax, kmid, ntot) )
    print( 'method = magnus, generator = %s' % generator )
    print( 'k_magnus = %d, ds = %.1e' % (k_magnus, ds) )
    
    # Writes evolved potentials and omega matrices to files if save = True
    if save:

        # Save evolved potential and omega matrix for each lambda value
        for lamb in lambda_array:

            # Scattering units here [fm^-2]
            H_evolved = d['hamiltonian'][lamb]
            # Subtract off kinetic energy (need to convert T_rel from MeV to
            # fm^-2)
            V_evolved = H_evolved - T_rel / hbar_sq_over_M
            
            # Save evolved potential
            lp.save_potential(k_array, k_weights, V_evolved, kvnn, channel, 
                              kmax, kmid, ntot, 'magnus', generator, lamb,
                              k_magnus=k_magnus, ds=ds)
            
            O_evolved = d['omega'][lamb]
            # Save evolved omega matrix
            lp.save_omega(k_array, O_evolved, kvnn, channel, kmax, kmid, ntot,
                          generator, lamb, k_magnus, ds=ds)
                
    # Otherwise, only return the dictionary d
    return d


# Make run_magnus.py executable
if __name__ == '__main__':
    
    
    # Specify potential
    
    #kvnn = 10 # EM N3LO
    #kvnn = 112 # RKE N4LO at Lambda = 500 MeV
    #kvnn = 222  # Gezerlis et al N2LO local potential at R_0 = 1.0 fm cutoff
    kvnn = 900 # Wendt at Lambda = 4 fm^-1
    #kvnn = 901 # Wendt at Lambda = 9 fm^-1
    #kvnn = 902 # Wendt at Lambda = 20 fm^-1
    
    channel = '3S1'
    
    kmax = 30.0
    #kmax = 8.0
    #kmax = 10.0
    
    kmid = 4.0
    #kmid = 2.0
    
    ntot = 120
    
    # Specify evolution
    
    generator = 'Wegner'
    #generator = 'T'
    
    # Evolve to lambda = 10 fm^-1 separately at a step-size 10^-6
    lambda_array = np.array( [2.8, 2.0, 1.2] )
    #lambda_array = np.array( [10.0] )
    #lambda_array = np.array( [6.0, 3.0, 2.0, 1.5] )
    #lambda_array = np.array( [2.8] )
    #lambda_array = np.array( [1.2] )
    
    k_magnus = 2
    #k_magnus = 4
    #k_magnus = 6
    #k_magnus = 10
    #k_magnus = 14
    
    ds = 1e-5
    #ds = 1e-6
    
    # Save the evolved Hamiltonian and omega?
    save = True
    #save = False
    
    # Evolve Hamiltonian
    d = run_magnus(kvnn, channel, kmax, kmid, ntot, generator, lambda_array,
                   k_magnus, ds, save)