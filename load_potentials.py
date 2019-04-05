# Created 03/22/19 by A.T. (tropiano.4@osu.edu)

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

from os import getcwd,chdir
import numpy as np


def load_mesh(kvnn,channel,kmax,kmid,ntot):
    '''Load momentum and weights for given potential.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
        
    potential_directory = 'Potentials/vsrg_macos/'+ \
    'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d'%(kvnn_string,kmax,kmid,ntot)
    
    chdir(potential_directory)
    
    # Initialize arrays
    k_array = np.zeros(ntot)
    k_weights = np.zeros(ntot)
    
    # Open mesh file and load
    mesh_file = 'vsrg_%s_kvnn_%s_lam12.0_reg_0_3_0_mesh.out'%(channel,kvnn_string)
    f = open(mesh_file,'r')
    
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
        
    return k_array,k_weights
    

def load_V(kvnn,channel,kmax,kmid,ntot):
    '''Loads an NN potential in units fm.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    
    # NEXT UPDATE: Generalize to SRG/Magnus evolved files
    
    # Get current working directory
    cwd = getcwd()
    
    # Go to directory of specified potential
    # Convert kvnn to string
    if kvnn < 10:
        kvnn_string = '0'+str(kvnn)
    else:
        kvnn_string = str(kvnn)
    
    potential_directory = 'Potentials/vsrg_macos/'+ \
    'vsrg_kvnn_%s_lam12.0_kmax%d_kmid%d_ntot%d'%(kvnn_string,kmax,kmid,ntot)
    
    chdir(potential_directory)
    
    # Open potential file and load
    vnn_file = 'vnn_%s_kvnn_%s_lam12.0_reg_0_3_0.out'%(channel,kvnn_string)
        
    # Load output file
    data = np.loadtxt(vnn_file)

    chdir(cwd)
    
    # Coupled channel potential
    if coupled_channel(channel):
        v11 = np.reshape(data[:,2],(ntot,ntot))
        v12 = np.reshape(data[:,3],(ntot,ntot))
        v21 = np.reshape(data[:,4],(ntot,ntot))
        v22 = np.reshape(data[:,5],(ntot,ntot))
        V = np.vstack((np.hstack((v11,v12)),np.hstack((v21,v22))))
    else:
        V = np.reshape(data[:,2],(ntot,ntot))

    return V


def load_T(kvnn,channel,kmax,kmid,ntot):
    '''Loads relative kinetic energy in units MeV.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # Load momentum array
    k_array = load_mesh(kvnn,channel,kmax,kmid,ntot)[0]
    # Matrix of (h-bar*k)^2 / M along diagonal (n x n)
    T = np.diag(k_array**2)*hbar_sq_over_M
    
    # Coupled channel potential
    if coupled_channel(channel):
        
        # Length of k_array
        n = len(k_array)
        # Matrix of zeros (n x n)
        o = np.zeros((n,n))
        T = np.vstack((np.hstack((T,o)),np.hstack((o,T))))
        
    return T
    

def load_H(kvnn,channel,kmax,kmid,ntot):
    '''Loads Hamiltonian for given potential in units MeV.'''
    
    # Arguments
    
    # kvnn (integer): This number specifies the potential
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    # kmax (float): Maximum value in the momentum mesh
    # kmid (float): Mid-point value in the momentum mesh
    # ntot (integer): Number of momentum points in mesh
    
    # Load relative kinetic energy in units MeV
    T = load_T(kvnn,channel,kmax,kmid,ntot)
    
    # Load potential in units fm
    V = load_V(kvnn,channel,kmax,kmid,ntot)
    
    # Load momentum and weights
    k_array,k_weights = load_mesh(kvnn,channel,kmax,kmid,ntot)
    
    # Convert to units MeV
    V = convert2MeV(k_array,k_weights,V,coupled_channel(channel))
    
    # Add T and V for Hamiltonian
    return T+V


def coupled_channel(channel):
    '''Returns the truth value on whether the given channel is a coupled
    channel.'''
    
    # Arguments
    
    # channel (string): The partial wave channel ('1S0', '3S1', etc.)
    
    # List of coupled channels
    coupled_channels = ['3S1','3P2','3D3','3F4','3G5','3H6','3I7','3K8','3L9',\
                        '3M10','3N11','3O12','3Q13']

    # This is true or false
    boolean_value = channel in coupled_channels
    return boolean_value


def convert2MeV(k_array,k_weights,V,coupled_channel=False):
    '''Converts a potential from units fm to MeV with momentum array and 
    weights.'''
    
    # Arguments
    
    # k_array (NumPy array): Momentum array
    # k_weights (NumPy array): Momentum weights
    # V (NumPy array): Potential matrix in units fm
    # coupled_channel (Boolean): Value corresponding to whether the potential
    # is coupled channel or not
    
    # h-bar^2 / M [MeV fm^2]
    hbar_sq_over_M = 41.47
    
    # If coupled channel, double the length of the momentum and weights arrays
    if coupled_channel:
        gp = np.concatenate((k_array,k_array))
        gw = np.concatenate((k_weights,k_weights))
    else:
        gp = k_array
        gw = k_weights
    
    # Build grids of momentum and weights factor
    row,col = np.meshgrid(gp*np.sqrt(gw),gp*np.sqrt(gw))
    
    # Multiply the potential by 2/pi*row*col gives fm^-2 conversion
    V *= 2/np.pi*row*col
    
    # Multiply by h-bar^2 / M [MeV fm^2] to obtain MeV units
    return V*hbar_sq_over_M


if __name__ == '__main__':
    
    # Example for loading momentum mesh and weights below
    V = load_V(6,'3S1',30.0,4.0,120)