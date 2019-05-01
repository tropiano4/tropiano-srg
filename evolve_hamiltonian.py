# Created 05/01/19 by A.T. (tropiano.4@osu.edu)

# Run SRG or Magnus. This script runs an SRG or Magnus evolution for a given
# Hamiltonian in units MeV.


import time
import numpy as np
# Scripts made by A.T.
import load_save_potentials as L
import srg_wegner
import srg_kinetic_energy
import srg_block_diagonal
#import magnus_wegner
#import magnus_kinetic_energy


def main(kvnn, channel, kmax, kmid, ntot, method, generator, lambda_array, \
         lambda_bd=0.00, save=True):
    '''This function evolves a specified Hamiltonian to several values of
    lambda [fm^-1] and has the option to save evolved potentials.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    # method (string): The evolution method ('srg' or 'magnus')
    # generator (string): SRG generator ('Wegner', 'T', 'Block-diag')
    # lambda_array (NumPy array): Array of lambda values to be evolved to
    # lambda_bd (float): Lambda value for block-diagonal decoupling (e.g. 2.00 
    # fm^-1)
    # save (Boolean): Option to save evolved potentials. If true, writes the
    # data formatted as initial potentials.
    
    # Load Hamiltonian
    H0_matrix = L.load_H(kvnn, channel, kmax, kmid, ntot)
    
    # Initialize SRG class
    if generator == 'Wegner':
        
        srg = srg_wegner.SRG(H0_matrix)
      
    elif generator == 'T':

        T0_matrix = L.load_T(kvnn, channel, kmax, kmid, ntot)     
        srg = srg_kinetic_energy.SRG(H0_matrix, T0_matrix)
        
    elif generator == 'Block-diag':
        
        k_array = L.load_mesh(kvnn, channel, kmax, kmid, ntot)[0]
        cc = L.coupled_channel(channel)
        srg = srg_block_diagonal.SRG(H0_matrix, lambda_bd, k_array, cc)
        
    # Time the evolution and return dictionary d
    t0 = time.time() # Start time
    d = srg.evolve_hamiltonian(lambda_array)
    t1 = time.time() # End time
    
    mins = round((t1-t0)/60.0,2) # Minutes elapsed evolving H(s)
    print('_'*80)
    print('H(s) done evolving to lambda = %.1f fm^-1 after %f minutes'%( \
          lambda_array[-1], mins))
    print('_'*80)
    print('Specifications:\n')
    print('kvnn = %d, channel = %s'%(kvnn, channel))
    print('kmax = %.1f, kmid = %.1f, ntot = %d'%(kmax, kmid, ntot))
    print('method = %s, generator = %s'%(method, generator))
    if generator == 'Block-diag':
        print('block-diagonal lambda = %.1f'%lambda_bd)
    
    # Writes evolved potentials to files if save = True
    if save:
        
        # Use a function from load_save_potentials here
        a = 1
        
    # Otherwise, return the dictionary d
    else:
        
        return d
    

if __name__ == '__main__':
    
    # Example: Evolve EM N3LO
    
    # Specify potential
    kvnn = 10
    channel = '3S1'
    kmax = 30.0
    kmid = 4.0
    ntot = 120
    
    # Specify evolution
    
    method = 'srg'
    #method = 'magnus'
    
    generator = 'Wegner'
    #generator = 'T'
    #generator = 'Block-diag'
    lambda_bd = 2.00
    
    #lambda_array = np.array([10.0,2.8,2.0,1.2])
    #lambda_array = np.array([10.0,2.8])
    lambda_array = np.array([10.0])

    # Run evolution with saving
    #main(kvnn, channel, kmax, kmid, ntot, method, generator, lambda_array, \
    #     lambda_bd)
    
    # Run evolution without saving
    d = main(kvnn, channel, kmax, kmid, ntot, method, generator, \
             lambda_array, lambda_bd, save=False)