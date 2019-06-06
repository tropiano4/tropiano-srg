#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: load_save_potentials.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 22, 2019
# 
# Loads potentials from Potentials/vsrg_macos directory.
# Potentials are organized by a kvnn number and the details of their momentum
# mesh (kmax, kmid, and ntot). Below are several examples of kvnn numbers:
#
# kvnn   description        
#                         
#   6     Argonne v18                                
#  10     Entem/Machleidt N3LO (500 MeV cutoff)      
#  12     Entem/Machleidt N3LO (600 MeV cutoff)
#  13     Entem/Machleidt N3LOW (400 MeV cutoff)
#  32     Epelbaum et al N3LO (550/600 MeV cutoff)
#  90     Reinert-Krebs-Epelbaum LO (400 MeV cutoff)
# 222     Gezerlis et al N2LO local potential at R_0 = 1.0 fm cutoff
# 224     Gezerlis et al N2LO local potential at R_0 = 1.2 fm cutoff
# 900     Wendt LO non-local potential at 4 fm^-1 cutoff
# 901     Wendt LO non-local potential at 9 fm^-1 cutoff
# 902     Wendt LO non-local potential at 20 fm^-1 cutoff
#
# Also saves SRG or Magnus evolved potentials in the same directory.
#
# Revision history:
#   May 7, 2019 --- Save and load Magnus evolved potentials and omega matrices.
#
#------------------------------------------------------------------------------


from os import getcwd, chdir
import numpy as np


def load_momentum(kvnn, channel, kmax, kmid, ntot):
    """
    Load momentum and weights for given potential

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
        
    Returns
    -------
    k_array : 1-D ndarray
        Momentum array.
    k_weights : 1-D ndarray
        Momentum weights.
        
    """
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/'+ \
    'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d'%(kvnn_string, kmax, kmid, ntot)
    
    chdir(potential_directory)
    
    # Initialize arrays
    k_array = np.zeros(ntot)
    k_weights = np.zeros(ntot)
    
    # Open mesh file and load
    mesh_file = 'vsrg_%s_kvnn_%s_lam12.0_reg_0_3_0_mesh.out'%(channel,
                                                              kvnn_string)
    f = open(mesh_file, 'r')
    
    i = 0
    for line in f:
        unit = line.strip().split()
        k_array[i] = unit[0]
        k_weights[i] = unit[1]
        i += 1
    
    # Close file
    f.close()
      
    # Go back to original directory
    chdir(cwd)
        
    return k_array, k_weights
    

def load_potential(kvnn, channel, kmax, kmid, ntot, method='initial', 
                   generator='Wegner', lamb=1.5, lambda_bd=0.00, k_magnus=6, 
                   ds=1e-5):
    """
    Loads an NN potential in units fm.
    
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
    method : str, optional
        The evolution method 'srg' or 'magnus'. Choose 'initial' if unevolved.
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda in units fm^-1.
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling (e.g. 2.00 fm^-1).
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only).
        
    Returns
    -------
    V : 2-D ndarray
        NN potential matrix in units fm.
        
    """
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
    
    potential_directory = 'Potentials/vsrg_macos/'+ \
    'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d'%(kvnn_string, kmax, kmid, ntot)
    
    # Name of potential file (split on cases of initial, SRG or Magnus evolved)
    if method == 'initial':
        
        vnn_file = 'vnn_%s_kvnn_%s_lam12.0_reg_0_3_0.out'%(channel, 
                                                           kvnn_string)
        
    elif method == 'srg': 
    
        if generator == 'Block-diag':
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f.out'%(channel,
                        kvnn_string, method, generator, lambda_bd, lamb)
        else: 
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f.out'%(channel,
                        kvnn_string, method, generator, lamb)
            
    elif method == 'magnus': 
        
        if generator == 'Block-diag':
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f_k%d_ds%.1e.out'%(
                        channel, kvnn_string, method, generator, lambda_bd,
                        lamb, k_magnus, ds)
        else: 
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f_k%d_ds%.1e.out'%(
                        channel, kvnn_string, method, generator, lamb,
                        k_magnus, ds)
    
    chdir(potential_directory)
        
    # Load output file
    data = np.loadtxt(vnn_file)

    chdir(cwd)
    
    # Coupled channel potential
    if coupled_channel(channel):
        v11 = np.reshape(data[:,2], (ntot,ntot))
        v12 = np.reshape(data[:,3], (ntot,ntot))
        v21 = np.reshape(data[:,4], (ntot,ntot))
        v22 = np.reshape(data[:,5], (ntot,ntot))
        V = np.vstack( ( np.hstack( (v11, v12) ), np.hstack( (v21, v22) ) ) )
    else:
        V = np.reshape(data[:,2], (ntot, ntot))

    return V


def load_kinetic_energy(kvnn, channel, kmax, kmid, ntot):
    """
    Loads relative kinetic energy in units MeV.

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
        
    Returns
    -------
    T : 2-D ndarray
        Relative kinetic energy matrix in units MeV.
        
    """
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # Load momentum array
    k_array = load_momentum(kvnn, channel, kmax, kmid, ntot)[0]
    # Matrix of (h-bar*k)^2 / M along diagonal (n x n)
    T = np.diag(k_array**2)*hbar_sq_over_M
    
    # Coupled channel potential
    if coupled_channel(channel):
        
        # Length of k_array
        n = len(k_array)
        # Matrix of zeros (n x n)
        o = np.zeros((n, n))
        T = np.vstack( ( np.hstack( (T, o) ), np.hstack( (o, T) ) ) )
        
    return T
    

def load_hamiltonian(kvnn, channel, kmax, kmid, ntot, method='initial', 
                     generator='Wegner', lamb=1.2, lambda_bd=0.00, k_magnus=6, 
                     ds=1e-5):
    """
    Loads Hamiltonian for given potential in units MeV.
    
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
    method : str, optional
        The evolution method 'srg' or 'magnus'. Choose 'initial' if unevolved.
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda in units fm^-1.
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling (e.g. 2.00 fm^-1).
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only).
        
    Returns
    -------
    H : 2-D ndarray
        Hamiltonian matrix in units MeV.
        
    """
    
    # Load relative kinetic energy in units MeV
    T = load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)
    
    # Load potential in units fm
    V = load_potential(kvnn, channel, kmax, kmid, ntot, method, generator, \
                       lamb, lambda_bd, k_magnus, ds)
    
    # Load momentum and weights
    k_array, k_weights = load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # Convert to units MeV
    V = convert2MeV(k_array, k_weights, V, coupled_channel(channel))
    
    # Add T and V for Hamiltonian
    H = T+V
    
    return H


def load_omega(kvnn, channel, kmax, kmid, ntot, generator, lamb, k_magnus=6,
               ds=1e-5):
    """
    Loads a Magnus evolved omega matrix.
    
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
    lamb : float
        Evolution parameter lambda in units fm^-1.
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... )
    ds : float, optional
        Step-size in the flow parameter s.
        
    Returns
    -------
    O : 2-D ndarray
        Magnus evolved omega matrix.
        
    """
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
    
    potential_directory = 'Potentials/vsrg_macos/'+ \
    'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d'%(kvnn_string, kmax, kmid, ntot)
    
    # Name of omega file
    omega_file = 'omega_%s_kvnn_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
                  (channel, kvnn_string, generator, lamb, k_magnus, ds)
    
    chdir(potential_directory)
        
    # Load output file
    data = np.loadtxt(omega_file)

    chdir(cwd)
    
    # Coupled channel potential
    if coupled_channel(channel):
        o11 = np.reshape(data[:,2], (ntot,ntot))
        o12 = np.reshape(data[:,3], (ntot,ntot))
        o21 = np.reshape(data[:,4], (ntot,ntot))
        o22 = np.reshape(data[:,5], (ntot,ntot))
        O = np.vstack( ( np.hstack( (o11, o12) ), np.hstack( (o21, o22) ) ) )
    else:
        O = np.reshape(data[:,2], (ntot, ntot))

    return O


def save_potential(k_array, k_weights, V, kvnn, channel, kmax, kmid, ntot, 
                   method, generator, lamb, lambda_bd=0.00, k_magnus=6, 
                   ds=1e-5):
    """
    Saves an SRG or Magnus evolved potential in units fm.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array.
    k_weights : 1-D ndarray
        Momentum weights.
    V : 2-D ndarray
        Potential matrix in units fm^-2.
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
    method : str
        The evolution method 'srg' or 'magnus'.
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float
        Evolution parameter lambda in units fm^-1.
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling (e.g. 2.00 fm^-1).
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only).
        
    """

    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/'+ \
    'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d'%(kvnn_string, kmax, kmid, ntot)
    
    chdir(potential_directory)
    
    # Name potential file and save
    if method == 'srg': 
    
        if generator == 'Block-diag':
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f.out'%(channel,
                        kvnn_string, method, generator, lambda_bd, lamb)
        else: 
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f.out'%(channel, 
                        kvnn_string, method, generator, lamb)
            
    elif method == 'magnus': 

        vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
                    (channel, kvnn_string, method, generator, lamb, k_magnus,
                     ds)
    
    f = open(vnn_file,'w')
    
    # Length of momentum array
    n = len(k_array)
    
    # Coupled-channel potentials are written differently - write each sub-
    # block as a column
    if coupled_channel(channel):
        
        header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format('k',
                  'kp', 'V11', 'V12', 'V21', 'V22')
        f.write('#'+header+'\n')
        
        for i in range(n):
            
            k = k_array[i]
            w = k_weights[i]
            
            for j in range(n):
                
                kp = k_array[j]
                wp = k_weights[j]
                # For conversion to fm multiply by this factor
                factor = np.pi/(2.0*k*kp*np.sqrt(w*wp))
                
                v11 = V[i,j]*factor
                v12 = V[i,j+n]*factor
                v21 = V[i+n,j]*factor
                v22 = V[i+n,j+n]*factor
                
                line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
                        k, kp, v11, v12, v21, v22)
                
                f.write(line+'\n')
                
    else:
        
        header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'V')
        f.write('#'+header+'\n')
        
        for i in range(n):
            
            k = k_array[i]
            w = k_weights[i]
            
            for j in range(n):
                
                kp = k_array[j]
                wp = k_weights[j]
                # For conversion to fm multiply by this factor
                factor = np.pi/(2.0*k*kp*np.sqrt(w*wp))
                
                v = V[i,j]*factor
                
                line = '{:^15.6f}{:^15.6f}{:^23e}'.format(k, kp, v)
                
                f.write(line+'\n')
                
    f.close()

    chdir(cwd)
    
    
def save_omega(k_array, O, kvnn, channel, kmax, kmid, ntot, generator, lamb,
               k_magnus, ds=1e-5):
    """
    Saves a Magnus evolved omega matrix.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array.
    O : 2-D ndarray
        Magnus evolved omega matrix.
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
    lamb : float
        Evolution parameter lambda in units fm^-1.
    k_magnus : int
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... )
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only).
        
    """

    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/'+ \
                          'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
                           (kvnn_string, kmax, kmid, ntot)
    
    chdir(potential_directory)
    
    # Name omega file and save
    omega_file = 'omega_%s_kvnn_%s_%s_lambda%.1f_k%d_ds%.1e.out'%(channel,
                  kvnn_string, generator, lamb, k_magnus, ds)
    
    f = open(omega_file,'w')
    
    # Length of momentum array
    n = len(k_array)
    
    if coupled_channel(channel):
        
        header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format('k', 'kp', 
                  'O11', 'O12', 'O21', 'O22')
        f.write('#'+header+'\n')
        
        for i in range(n):
            
            k = k_array[i]
            
            for j in range(n):
                
                kp = k_array[j]
                
                o11 = O[i,j]
                o12 = O[i,j+n]
                o21 = O[i+n,j]
                o22 = O[i+n,j+n]
                
                line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
                        k, kp, o11, o12, o21, o22)
                
                f.write(line+'\n')
                
    else:
        
        header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'O')
        f.write('#'+header+'\n')
        
        for i in range(n):
            
            k = k_array[i]
            
            for j in range(n):
                
                kp = k_array[j]
                
                o = O[i,j]
                
                line = '{:^15.6f}{:^15.6f}{:^23e}'.format(k, kp, o)
                
                f.write(line+'\n')
                
    f.close()

    chdir(cwd)


def coupled_channel(channel):
    """
    Truth value on whether the given channel is a coupled channel.
    
    Parameters
    ----------
    channel : str
        The partial wave channel ('1S0', '3S1', etc.)
    
    Returns
    -------
    boolean_value : bool
        True if the channel is coupled channel and false otherwise.
        
    """
    
    # List of coupled channels
    coupled_channels = ['3S1','3P2','3D3','3F4','3G5','3H6','3I7','3K8','3L9',
                        '3M10','3N11','3O12','3Q13']

    # This is true or false
    boolean_value = channel in coupled_channels
    
    return boolean_value


def convert2MeV(k_array, k_weights, V_fm, coupled_channel=False):
    """
    Converts a potential from units fm to MeV with momentum array and 
    weights.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array.
    k_weights : 1-D ndarray
        Momentum weights.
    V_fm : 2-D ndarray
        Potential matrix in units fm.
    coupled_channel : bool, optional
        True if the channel is coupled channel and false otherwise.

    Returns
    -------
    V_MeV : 2-D ndarray
        Potential matrix in units MeV.

    """
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # If coupled channel, double the length of the momentum and weights arrays
    if coupled_channel:
        gp = np.concatenate( (k_array, k_array) )
        gw = np.concatenate( (k_weights, k_weights) )
    else:
        gp = k_array
        gw = k_weights
    
    # Build grids of momentum and weights factor
    row, col = np.meshgrid(gp*np.sqrt(gw), gp*np.sqrt(gw))
    
    # Multiply the potential by 2/pi*row*col gives fm^-2 conversion and 
    # multiplying by hbar_sq_over_M gives the MeV conversion
    V_MeV = V_fm * 2/np.pi * row * col * hbar_sq_over_M
    
    return V_MeV