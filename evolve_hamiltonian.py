# Created 05/01/19 by A.T. (tropiano.4@osu.edu)
# Updated 05/07/19 by A.T. to run Magnus implementation

# Run SRG or Magnus. This script runs an SRG or Magnus evolution for a given
# Hamiltonian in units MeV. SRG code return a dictionary of the evolved 
# Hamiltonian (in units fm^-2) at each point in lambda (which serves as the 
# key). Magnus code returns the evolved Hamiltonian and Omega matrix at a given
# value of lambda.


import numpy as np
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import load_save_potentials as lp
from SRG_codes import srg_wegner
from SRG_codes import srg_kinetic_energy
from SRG_codes import srg_block_diagonal
from Magnus_codes import magnus_wegner
from Magnus_codes import magnus_kinetic_energy


def run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array, \
            lambda_bd=0.00, save=True):
    '''This function SRG evolves a specified Hamiltonian to several values of
    lambda [fm^-1] and has the option to save the evolved potentials.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    # generator (string): SRG generator ('Wegner', 'T', 'Block-diag')
    # lambda_array (1-D NumPy array): Array of lambda evolution values 
    # lambda_bd (float): Lambda value for block-diagonal decoupling (e.g. 2.00 
    # fm^-1)
    # save (Boolean): Option to save evolved potentials. If true, writes the
    # evolved potential formatted similar to the initial potentials.
    
    # Load Hamiltonian, kinetic energy and weights
    H0_matrix = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    T0_matrix = lp.load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)     
    k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # h-bar^2 / M [MeV fm^2] for conversion of MeV to fm^-2
    hbar_sq_over_M = 41.47

    # Initialize SRG class
    if generator == 'Wegner':
        
        evolve = srg_wegner.SRG(H0_matrix)
      
    elif generator == 'T':

        evolve = srg_kinetic_energy.SRG(H0_matrix, T0_matrix)
        
    elif generator == 'Block-diag':
        
        cc = lp.coupled_channel(channel)
        evolve = srg_block_diagonal.SRG(H0_matrix, lambda_bd, k_array, cc)
        
    # Time the evolution and return dictionary d of evolved Hamiltonians
    t0 = time.time() # Start time
    d = evolve.evolve_hamiltonian(lambda_array)
    t1 = time.time() # End time
    
    mins = round((t1-t0)/60.0,2) # Minutes elapsed evolving H(s)
    print('_'*85)
    print('H(s) done evolving to final lambda = %.1f fm^-1 after %f minutes'%( \
          lambda_array[-1], mins))
    print('_'*85)
    print('\nSpecifications:\n')
    print('kvnn = %d, channel = %s'%(kvnn, channel))
    print('kmax = %.1f, kmid = %.1f, ntot = %d'%(kmax, kmid, ntot))
    print('method = srg, generator = %s'%generator)
    if generator == 'Block-diag':
        print('block-diagonal lambda = %.2f fm^-1'%lambda_bd)
    
    # Writes evolved potentials to files if save = True
    if save:

        # Save evolved potential for each lambda value
        for lamb in lambda_array:

            Hs_matrix = d[lamb] # Units are fm^-2 here
            # Subtract off kinetic energy (need to convert from MeV to fm^-2)
            Vs_matrix = Hs_matrix - T0_matrix/hbar_sq_over_M
            # Save evolved potential
            lp.save_potential(k_array, k_weights, Vs_matrix, kvnn, channel, kmax, \
                              kmid, ntot, 'srg', generator, lamb, lambda_bd)
                
    # Otherwise, return the dictionary d
    else:
        
        return d
    
    
def run_magnus(kvnn, channel, kmax, kmid, ntot, method, generator, lambda_array, \
               k_magnus=6, ds=1e-5, save=True):
    '''This function evolves a specified Hamiltonian to several values of
    lambda [fm^-1] and has the option to save the evolved potentials.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    # method (string): The evolution method ('srg' or 'magnus')
    # generator (string): SRG generator ('Wegner', 'T', 'Block-diag')
    # lambda_array (1-D NumPy array): Array of lambda evolution values 
    # k_magnus (integer): Number of terms to include in Magnus sum
    # ds (float): Step-size in the flow parameter s
    # save (Boolean): Option to save evolved potentials. If true, writes the
    # evolved potential formatted similar to the initial potentials.
    
    # Load Hamiltonian, kinetic energy and weights
    H0_matrix = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    T0_matrix = lp.load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)     
    k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # h-bar^2 / M [MeV fm^2] for conversion of MeV to fm^-2
    hbar_sq_over_M = 41.47

    # Initialize Magnus class
    if generator == 'Wegner':
            
        evolve = magnus_wegner.Magnus(H0_matrix, k_magnus, ds)
            
    elif generator == 'T':
            
        evolve = magnus_kinetic_energy.Magnus(H0_matrix, T0_matrix, k_magnus, ds)
        
    # Time the evolution and return dictionary d of evolved Hamiltonians and 
    # omega matrices
    d = {}
    d['hamiltonian'] = {}
    d['omega'] = {}
    
    t0 = time.time() # Start time
    for lamb in lambda_array:
        Hs_matrix, Os_matrix = evolve.evolve_hamiltonian(lamb)
        d['hamiltonian'][lamb] = Hs_matrix
        d['omega'][lamb] = Os_matrix
    t1 = time.time() # End time
    
    mins = round((t1-t0)/60.0,2) # Minutes elapsed evolving H(s)
    print('_'*85)
    print('H(s) done evolving to final lambda = %.1f fm^-1 after %f minutes'%( \
          lambda_array[-1], mins))
    print('_'*85)
    print('\nSpecifications:\n')
    print('kvnn = %d, channel = %s'%(kvnn, channel))
    print('kmax = %.1f, kmid = %.1f, ntot = %d'%(kmax, kmid, ntot))
    print('method = magnus, generator = %s'%generator)
    print('k_magnus = %d, ds = %.1e'%(k_magnus, ds))
    
    # Writes evolved potentials to files if save = True
    # Otherwise, just return the dictionary d
    if save:
                
        # Save evolved potential and Omega for each lambda value
        for lamb in lambda_array:
            
            # Need a smaller step-size for large values of lambda
            if lamb >= 10.0 and ds > 1e-6:
                step_size = 1e-6
            else:
                step_size = ds

            Hs_matrix = d['hamiltonian'][lamb] # Units are fm^-2 here
            # Subtract off kinetic energy (need to convert from MeV to fm^-2)
            Vs_matrix = Hs_matrix - T0_matrix/hbar_sq_over_M
            # Save evolved potential
            lp.save_potential(k_array, k_weights, Vs_matrix, kvnn, channel, \
                              kmax, kmid, ntot, 'magnus', generator, lamb, \
                              k_magnus=k_magnus, ds=step_size)
                
            Os_matrix = d['omega'][lamb]
            lp.save_omega(k_array, Os_matrix, kvnn, channel, kmax, kmid, ntot, \
                          generator, lamb, k_magnus=k_magnus, ds=step_size)

    return d
    

if __name__ == '__main__':
    
    
    # Specify potential
    
    #kvnn = 6 # AV18
    #kvnn = 10 # EM N3LO
    #kvnn = 105 # RKE N3LO at Lambda = 400 MeV
    #kvnn = 106 # RKE N3LO at Lambda = 450 MeV
    #kvnn = 107 # RKE N3LO at Lambda = 500 MeV
    #kvnn = 111 # RKE N4LO at Lambda = 450 MeV
    #kvnn = 112 # RKE N4LO at Lambda = 500 MeV
    #kvnn = 119 # Wendt at Lambda = 4 fm^-1
    #kvnn = 120 # Wendt at Lambda = 9 fm^-1
    #kvnn = 121 # Wendt at Lambda = 20 fm^-1
    #kvnn = 222  # Gezerlis et al N2LO local potential at R_0 = 1.0 fm cutoff
    kvnn = 224  # Gezerlis et al 2LO local potential at R_0 = 1.2 fm cutoff
    
    channel = '3S1'
    
    #kmax = 30.0
    #kmax = 8.0
    kmax = 10.0
    #kmid = 4.0
    kmid = 2.0
    ntot = 120
    
    
    # Specify evolution
    
    method = 'srg'
    #method = 'magnus'
    k_magnus = 6 # This won't affect SRG computations
    ds = 1e-5 # This won't affect SRG computations
    #ds = 1e-6
    
    generator = 'Wegner'
    #generator = 'T'
    #generator = 'Block-diag'
    #lambda_bd = 1.00 # This won't affect the band-diagonal generators
    lambda_bd = 2.00 
    #lambda_bd = 3.00
    #lambda_bd = 4.00
    
    #lambda_array = np.array([10.0,2.8,2.0,1.2])
    #lambda_array = np.array([10.0,2.8,2.0,1.5])
    #lambda_array = np.array([10.0,2.8])
    #lambda_array = np.array([10.0])
    #lambda_array = np.array([1.2])
    
    # Save the evolved Hamiltonian?
    save = True
    #save = False
    
    
    # Evolve Hamiltonian

    if method == 'srg':
        
        d = run_srg(kvnn, channel, kmax, kmid, ntot, generator, lambda_array, 
                    lambda_bd=lambda_bd, save=save)
        
    elif method == 'magnus':
        
        d = run_magnus(kvnn, channel, kmax, kmid, ntot, method, generator, 
                       lambda_array, k_magnus=k_magnus, ds=ds, save=save)