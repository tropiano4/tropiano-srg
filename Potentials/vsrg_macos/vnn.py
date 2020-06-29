#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: vnn.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 22, 2019
# 
# Loads potentials from Potentials/vsrg_macos directory. Potentials are 
# organized by a kvnn number and the details of their momentum mesh (kmax, 
# kmid, and ntot). Below are several examples of kvnn numbers:
#
# kvnn   description        
#                         
#   6     Argonne v18                                
#  10     Entem/Machleidt N3LO (500 MeV cutoff)      
#  12     Entem/Machleidt N3LO (600 MeV cutoff)
#  13     Entem/Machleidt N3LOW (400 MeV cutoff)
#  32     Epelbaum et al N3LO (550/600 MeV cutoff)
#  79     Entem-Machleidt-Nosyk N4LO (500 MeV cutoff)
#  90     Reinert-Krebs-Epelbaum LO (400 MeV cutoff)
# 222     Gezerlis et al N2LO local potential (1.0 fm cutoff)
# 224     Gezerlis et al N2LO local potential (1.2 fm cutoff)
# 900     Wendt LO non-local potential (4 fm^-1 cutoff)
# 901     Wendt LO non-local potential (9 fm^-1 cutoff)
# 902     Wendt LO non-local potential (20 fm^-1 cutoff)
#
# Also saves SRG or Magnus evolved potentials in the same directory. This
# script includes several other useful potential-oriented functions.
#
# Revision history:
#   05/07/19 --- Save and load Magnus evolved potentials and omega matrices.
#   02/26/20 --- Updated to include default momentum mesh specifications.
#   05/05/20 --- Updated to include function that converts potential to a new
#                momentum mesh.
#   06/29/20 --- Renamed to 'vnn.py'.
#
#------------------------------------------------------------------------------


from os import getcwd, chdir
import numpy as np
from scipy.interpolate import interp2d


def load_momentum(kvnn, channel, kmax=0.0, kmid=0.0, ntot=0):
    """
    Load momentum and weights for given potential.

    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
        
    Returns
    -------
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    k_weights : 1-D ndarray
        Momentum weights [fm^-1].
        
    """
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential

    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/vsrg_kvnn_' + \
                          '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
                          (kvnn_string, kmax, kmid, ntot)
    
    chdir(potential_directory)
    
    # Initialize arrays
    k_array = np.zeros(ntot)
    k_weights = np.zeros(ntot)
    
    # Open mesh file and load
    mesh_file = 'vsrg_%s_kvnn_%s_lam12.0_reg_0_3_0_mesh.out' % \
                (channel, kvnn_string)
    f = open(mesh_file, 'r')
    
    i = 0
    for line in f:
        unit = line.strip().split() # This is a list
        k_array[i] = unit[0] # First element of list
        k_weights[i] = unit[1] # Second element of list
        i += 1
    
    # Close file
    f.close()
      
    # Go back to original directory
    chdir(cwd)
        
    return k_array, k_weights
    

def load_potential(kvnn, channel, kmax=0.0, kmid=0.0, ntot=0, method='initial', 
                   generator='Wegner', lamb=1.5, lambda_bd=0.00, k_magnus=6, 
                   ds=1e-5):
    """
    Loads an NN potential.
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
    method : str, optional
        The evolution method 'srg' or 'magnus'. Choose 'initial' if unevolved.
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda [fm^-1].
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling [fm^-1].
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only) [fm^4].
        
    Returns
    -------
    V : 2-D ndarray
        NN potential matrix [fm].
        
    """
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
    
    potential_directory = 'Potentials/vsrg_macos/vsrg_kvnn_' + \
                          '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
                          (kvnn_string, kmax, kmid, ntot)
    
    # Name of potential file (split on cases of initial, SRG or Magnus evolved)
    if method == 'initial':
        
        vnn_file = 'vnn_%s_kvnn_%s_lam12.0_reg_0_3_0.out' % \
                   (channel, kvnn_string)
        
    elif method == 'srg': 
    
        if generator == 'Block-diag':
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f.out' % (channel, 
                        kvnn_string, method, generator, lambda_bd, lamb)
        else: 
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f.out' % (channel,
                        kvnn_string, method, generator, lamb)
            
    elif method == 'magnus': 
        
        if generator == 'Block-diag':
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f_k%d_ds%.1e.out' % \
                       (channel, kvnn_string, method, generator, lambda_bd,
                        lamb, k_magnus, ds)
        else: 
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
                       (channel, kvnn_string, method, generator, lamb,
                        k_magnus, ds)
    
    chdir(potential_directory)
        
    # Load output file
    data = np.loadtxt(vnn_file)

    chdir(cwd)
    
    # Coupled channel potential
    if coupled_channel(channel):
        v11 = np.reshape(data[:,2], (ntot, ntot))
        v12 = np.reshape(data[:,3], (ntot, ntot))
        v21 = np.reshape(data[:,4], (ntot, ntot))
        v22 = np.reshape(data[:,5], (ntot, ntot))
        V = np.vstack( ( np.hstack( (v11, v12) ), np.hstack( (v21, v22) ) ) )
    else:
        V = np.reshape(data[:,2], (ntot, ntot))

    return V


def load_kinetic_energy(kvnn, channel, kmax=0.0, kmid=0.0, ntot=0):
    """
    Loads relative kinetic energy.

    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
        
    Returns
    -------
    T : 2-D ndarray
        Relative kinetic energy matrix [MeV].
        
    """
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
    
    # Load momentum array
    k_array = load_momentum(kvnn, channel, kmax, kmid, ntot)[0]
    
    # Matrix of (h-bar*k)^2 / M along diagonal (n x n)
    T = np.diag(k_array**2) * hbar_sq_over_M
    
    # Coupled channel potential
    if coupled_channel(channel):
        
        # Length of k_array
        n = len(k_array)
        # Matrix of zeros (n x n)
        O = np.zeros( (n, n) )
        T = np.vstack( ( np.hstack( (T, O) ), np.hstack( (O, T) ) ) )
        
    return T
    

def load_hamiltonian(kvnn, channel, kmax=0.0, kmid=0.0, ntot=0, 
                     method='initial', generator='Wegner', lamb=1.5, 
                     lambda_bd=0.00, k_magnus=6, ds=1e-5):
    """
    Loads Hamiltonian for a given NN potential.
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
    method : str, optional
        The evolution method 'srg' or 'magnus'. Choose 'initial' if unevolved.
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda [fm^-1].
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling [fm^-1].
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only) [fm^4].
        
    Returns
    -------
    H : 2-D ndarray
        Hamiltonian matrix [MeV].
        
    """
    
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
    
    # Load relative kinetic energy in units MeV
    T = load_kinetic_energy(kvnn, channel, kmax, kmid, ntot)
    
    # Load potential in units fm
    V = load_potential(kvnn, channel, kmax, kmid, ntot, method, generator,
                       lamb, lambda_bd, k_magnus, ds)
    
    # Load momentum and weights
    k_array, k_weights = load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # Convert to units MeV
    V = convert2MeV(k_array, k_weights, V, coupled_channel(channel))
    
    # Add T and V for Hamiltonian
    H = T + V
    
    return H


def load_omega(kvnn, channel, kmax=0.0, kmid=0.0, ntot=0, generator='Wegner', 
               lamb=1.5, lambda_bd=0.00, k_magnus=6, ds=1e-5):
    """
    Loads a Magnus evolved omega matrix.
    
    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda [fm^-1].
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling [fm^-1].
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only) [fm^4].
        
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
        
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
    
    potential_directory = 'Potentials/vsrg_macos/vsrg_kvnn_' + \
                          '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
                          (kvnn_string, kmax, kmid, ntot)
    
    # Name of omega file
    omega_file = 'omega_%s_kvnn_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
                 (channel, kvnn_string, generator, lamb, k_magnus, ds)
    
    chdir(potential_directory)
        
    # Load output file
    data = np.loadtxt(omega_file)

    chdir(cwd)
    
    # Coupled channel potential
    if coupled_channel(channel):
        o11 = np.reshape(data[:,2], (ntot, ntot))
        o12 = np.reshape(data[:,3], (ntot, ntot))
        o21 = np.reshape(data[:,4], (ntot, ntot))
        o22 = np.reshape(data[:,5], (ntot, ntot))
        O = np.vstack( ( np.hstack( (o11, o12) ), np.hstack( (o21, o22) ) ) )
    else:
        O = np.reshape(data[:,2], (ntot, ntot))

    return O


def save_potential(k_array, k_weights, V, kvnn, channel, kmax=0.0, kmid=0.0, 
                   ntot=0, method='srg', generator='Wegner', lamb=1.5, 
                   lambda_bd=0.00, k_magnus=6, ds=1e-5):
    """
    Saves an SRG or Magnus evolved potential.
    
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
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
    method : str, optional
        The evolution method 'srg' or 'magnus'. Choose 'initial' if unevolved.
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda [fm^-1].
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling [fm^-1].
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only) [fm^4].
        
    """

    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/vsrg_kvnn_' + \
                          '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
                          (kvnn_string, kmax, kmid, ntot)
    
    chdir(potential_directory)
    
    # Name potential file and save
    if method == 'srg': 
    
        if generator == 'Block-diag':
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f.out' % (channel,
                        kvnn_string, method, generator, lambda_bd, lamb)
        else: 
            
            vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f.out' % (channel, 
                        kvnn_string, method, generator, lamb)
            
    elif method == 'magnus': 

        vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f_k%d_ds%.1e.out' % (channel, 
                    kvnn_string, method, generator, lamb, k_magnus, ds)
    
    else:
        
        vnn_file = 'vnn_%s_kvnn_%s_lam12.0_reg_0_3_0.out' % (channel, 
                                                             kvnn_string)
    
    f = open(vnn_file,'w')
    
    # Length of momentum array
    n = len(k_array)
    
    # Coupled-channel potentials are written differently - write each sub-
    # block as a column
    if coupled_channel(channel):
        
        header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format('k',
                  'kp', 'V11', 'V12', 'V21', 'V22')
        f.write('#' + header + '\n')
        
        for i in range(n):
            
            k = k_array[i]
            w = k_weights[i]
            
            for j in range(n):
                
                kp = k_array[j]
                wp = k_weights[j]
                # For conversion to fm multiply by this factor
                factor = np.pi / ( 2.0 * k * kp * np.sqrt(w*wp) )
                
                v11 = V[i, j] * factor
                v12 = V[i, j+n] * factor
                v21 = V[i+n, j] * factor
                v22 = V[i+n, j+n] * factor
                
                line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
                        k, kp, v11, v12, v21, v22)
                
                f.write(line + '\n')
                
    else:
        
        header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'V')
        f.write('#' + header + '\n')
        
        for i in range(n):
            
            k = k_array[i]
            w = k_weights[i]
            
            for j in range(n):
                
                kp = k_array[j]
                wp = k_weights[j]
                # For conversion to fm multiply by this factor
                factor = np.pi / ( 2.0 * k * kp * np.sqrt(w*wp) )
                
                v = V[i,j]*factor
                
                line = '{:^15.6f}{:^15.6f}{:^23e}'.format(k, kp, v)
                
                f.write(line+'\n')
                
    f.close()

    chdir(cwd)
    
    
def save_omega(k_array, O, kvnn, channel, kmax=0.0, kmid=0.0, ntot=0, 
               generator='Wegner', lamb=1.5, lambda_bd=0.00, k_magnus=6, 
               ds=1e-5):
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
        The partial wave channel (e.g. '1S0')
    kmax : float, optional
        Maximum value in the momentum mesh [fm^-1].
    kmid : float, optional
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int, optional
        Number of momentum points in mesh [fm^-1].
    generator : str, optional
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float, optional
        Evolution parameter lambda [fm^-1].
    lambda_bd : float, optional
        Lambda value for block-diagonal decoupling [fm^-1].
    k_magnus : int, optional
        Number of terms to include in Magnus sum (that is,
        dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
    ds : float, optional
        Step-size in the flow parameter s (for Magnus only) [fm^4].
        
    """

    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    # Set default mesh specifications if necessary
    if kmax == 0.0:
        kmax, kmid, ntot = mesh_specifications(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/vsrg_kvnn_' + \
                          '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
                          (kvnn_string, kmax, kmid, ntot)
    
    chdir(potential_directory)
    
    # Name omega file and save
    omega_file = 'omega_%s_kvnn_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
                 (channel, kvnn_string, generator, lamb, k_magnus, ds)
    
    f = open(omega_file,'w')
    
    # Length of momentum array
    n = len(k_array)
    
    if coupled_channel(channel):
        
        header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format('k', 'kp', 
                  'O11', 'O12', 'O21', 'O22')
        f.write('#' + header + '\n')
        
        for i in range(n):
            
            k = k_array[i]
            
            for j in range(n):
                
                kp = k_array[j]
                
                o11 = O[i, j]
                o12 = O[i, j+n]
                o21 = O[i+n, j]
                o22 = O[i+n, j+n]
                
                line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
                        k, kp, o11, o12, o21, o22)
                
                f.write(line + '\n')
                
    else:
        
        header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'O')
        f.write('#' + header + '\n')
        
        for i in range(n):
            
            k = k_array[i]
            
            for j in range(n):
                
                kp = k_array[j]
                
                o = O[i,j]
                
                line = '{:^15.6f}{:^15.6f}{:^23e}'.format(k, kp, o)
                
                f.write(line + '\n')
                
    f.close()

    chdir(cwd)


def coupled_channel(channel):
    """
    Truth value on whether the given channel is a coupled channel.
    
    Parameters
    ----------
    channel : str
        The partial wave channel (e.g. '1S0')
    
    Returns
    -------
    boolean_value : bool
        True if the channel is coupled channel and false otherwise.
        
    """
    
    # List of coupled channels
    coupled_channels = ['3S1', '3P2', '3D3', '3F4', '3G5', '3H6', '3I7',
                        '3K8', '3L9', '3M10', '3N11', '3O12', '3Q13']

    # This is true or false
    boolean_value = channel in coupled_channels
    
    return boolean_value


def convert2MeV(k_array, k_weights, V_fm, coupled_channel=False):
    """
    Converts an NN potential from units fm to MeV with momentum array and 
    weights.
    
    Parameters
    ----------
    k_array : 1-D ndarray
        Momentum array [fm^-1.
    k_weights : 1-D ndarray
        Momentum weights [fm^-1.
    V_fm : 2-D ndarray
        Potential matrix [fm].
    coupled_channel : bool, optional
        True if the channel is a coupled channel and false otherwise.

    Returns
    -------
    V_MeV : 2-D ndarray
        Potential matrix [MeV].

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
    row, col = np.meshgrid(gp * np.sqrt(gw), gp * np.sqrt(gw))
    
    # Multiply the potential by 2/pi*row*col gives fm^-2 conversion and 
    # multiplying by hbar_sq_over_M gives the MeV conversion
    V_MeV = V_fm * 2/np.pi * row * col * hbar_sq_over_M
    
    return V_MeV


def mesh_specifications(kvnn):
    """
    Returns the default mesh specifications for a particular potential.

    Parameters
    ----------
    kvnn : int
        This number specifies the potential.

    Returns
    -------
    kmax : float
        Maximum value in the momentum mesh [fm^-1].
    kmid : float
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int
        Number of momentum points in mesh [fm^-1].

    """
    
    # In all cases so far, ntot = 120
    ntot = 120
    
    # Argonne v18 and Wendt LO (4, 9, 20 fm^-1)
    if kvnn in [10, 900, 901, 902]:
        kmax = 30.0
        kmid = 4.0
    # EM N3LO (500 MeV), EMN N4LO (450, 500, 550 MeV), RKE N3LO (400, 450, 500
    # MeV), RKE N4LO (400, 450, 500 MeV), Gezerlis N2LO (1, 1.2 fm)
    elif kvnn in [6, 74, 79, 84, 105, 106, 107, 110, 111, 112, 113, 222, 
                  224]:
        kmax = 10.0
        kmid = 2.0
    else:
        kmax = 10.0
        kmid = 2.0

    return kmax, kmid, ntot


def convert_potential_to_new_mesh(kvnn, channel, old_mesh, new_mesh):
    """
    Change a potential's momentum mesh and saves to new folder.

    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0')
    old_mesh : tuple
        Values which specify the old momentum mesh (kmax, kmid, ntot) where
        kmax and kmid are floats and ntot is an int.
    new_mesh : tuple
        Values which specify the new momentum mesh (kmax, kmid, ntot) where
        kmax and kmid are floats and ntot is an int.

    """
    
    # Load old mesh here
    kmax_old, kmid_old, ntot_old = old_mesh
    k_array_old, _ = load_momentum(kvnn, channel, kmax=kmax_old,
                                   kmid=kmid_old, ntot=ntot_old)
    
    # Interpolate V to old mesh
    # Start by loading V_matrix_old in units fm
    V_matrix_old = load_potential(kvnn, channel, kmax=kmax_old, kmid=kmid_old,
                                  ntot=ntot_old)
    # If coupled-channel, do one sub-block at a time
    if coupled_channel(channel):
        V_11_func = interp2d( k_array_old, k_array_old, 
                              V_matrix_old[:ntot_old, :ntot_old] )
        V_12_func = interp2d( k_array_old, k_array_old, 
                              V_matrix_old[:ntot_old, ntot_old:] )
        V_21_func = interp2d( k_array_old, k_array_old, 
                              V_matrix_old[ntot_old:, :ntot_old] )
        V_22_func = interp2d( k_array_old, k_array_old, 
                              V_matrix_old[ntot_old:, ntot_old:] )
    else:
        V_func = interp2d(k_array_old, k_array_old, V_matrix_old)
    
    # Load new mesh here
    kmax_new, kmid_new, ntot_new = new_mesh
    k_array_new, k_weights = load_momentum(kvnn, channel, kmax=kmax_new,
                                           kmid=kmid_new, ntot=ntot_new)
    
    # Generate new potential V_matrix_new (still units fm)
    # If coupled-channel, do one-sub-block at a time
    if coupled_channel(channel):
        V11 = V_11_func(k_array_new, k_array_new)
        V12 = V_12_func(k_array_new, k_array_new)
        V21 = V_21_func(k_array_new, k_array_new)
        V22 = V_22_func(k_array_new, k_array_new)
        V_matrix_new = np.vstack( ( np.hstack( (V11, V12) ), 
                                    np.hstack( (V21, V22) ) ) )
    else:
        V_matrix_new = V_func(k_array_new, k_array_new)
    
    # Convert to units fm^-2
    factor_array = np.sqrt(2 * k_weights / np.pi) * k_array_new
    # Double the length of factor_array if coupled-channel
    if coupled_channel(channel):
        factor_array = np.concatenate( (factor_array, factor_array) )
    row, col = np.meshgrid(factor_array, factor_array)
    V_matrix_new *= row * col
    
    # Save to file (this reads V in unit fm^-2!)
    save_potential(k_array_new, k_weights, V_matrix_new, kvnn, channel,
                   kmax=kmax_new, kmid=kmid_new, ntot=ntot_new,
                   method='initial')