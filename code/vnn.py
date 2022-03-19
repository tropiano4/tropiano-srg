#!/usr/bin/env python3

"""
File: vnn.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: March 22, 2019

Handles potentials from data/vnn directory. Potentials are organized by a kvnn
number, a partial wave channel, and the details of their momentum mesh. Below
are several examples of kvnn numbers:

    kvnn    description        
    
    5       Nijmegen II              
    6       Argonne v18   
    7       CD-Bonn                      
    10      Entem-Machleidt N3LO (500 MeV cutoff)      
    12      Entem-Machleidt N3LO (600 MeV cutoff)
    13      Entem-Machleidt N3LOW (400 MeV cutoff)
    32      Epelbaum et al N3LO (550/600 MeV cutoff)
    79      Entem-Machleidt-Nosyk N4LO (500 MeV cutoff)
    90      Reinert-Krebs-Epelbaum LO (400 MeV cutoff)
    111	    Reinert-Krebs-Epelbaum N4LO (450 MeV cutoff)
    222     GT+ N2LO (1.0 fm cutoff)
    224     GT+ N2LO (1.2 fm cutoff)
    900     Wendt LO non-local potential (4 fm^-1 cutoff)
    901     Wendt LO non-local potential (9 fm^-1 cutoff)
    902     Wendt LO non-local potential (20 fm^-1 cutoff)

Also runs SRG and Magnus scripts, saving evolved potentials in the same
directory as above. This class includes several other useful potential-
oriented functions.

Last update: March 18, 2022

"""

# To-do: Update potential_directory. How do you format long strings correctly?
# To-do: Change file naming convention to be more consistent across scripts.
# To-do: Use np.savetxt() to save potentials?
# To-do: Fix save_potential() - why is 'initial' even a possibility?
# To-do: Should probably incorporate dU/ds = \eta U method?

# Python imports
import numpy as np
import time

# Imports from A.T. codes
from modules.srg import SRG
import modules.tools as tl

class Potential(object):
    
    def __init__(self, kvnn, channel, kmax, kmid, ntot):
        """
        Save the inputs of the potential. Also records the relevant directory.

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

        """

        self.kvnn = kvnn
        self.channel = channel
        self.kmax = kmax
        self.kmid = kmid
        self.ntot = ntot
        
        # Save h-bar^2 / M for units conversion 
        self.hbar_sq_over_m = 41.47  # [MeV fm^2]
        
        # Need kvnn as string (can cause error if kvnn < 10)
        if kvnn < 10:
            kvnn_string = '0'+str(kvnn)
        else:
            kvnn_string = str(kvnn)
        self.kvnn_string = kvnn_string
        
        # Save whether the potential is coupled-channel or not
        self.coupled_channel_bool = tl.coupled_channel(self.channel)
            
        # Get potential directory
        potential_directory = '../potentials/vsrg_macos/vsrg_kvnn_' + \
                              '%s_lam12.0_kmax%d_kmid%d_ntot%d/' % \
                              (kvnn_string, kmax, kmid, ntot)
        self.potential_directory = potential_directory
        
    def load_mesh(self):
        """
        Loads the relative momentum and weights associated with the potential.

        Returns
        -------
        k_array : 1-D ndarray
            Momentum array [fm^-1].
        k_weights : 1-D ndarray
            Momentum weights [fm^-1].

        """

        # Load mesh file
        mesh_file = 'vsrg_%s_kvnn_%s_lam12.0_reg_0_3_0_mesh.out' % \
                    (self.channel, self.kvnn_string)    
        data = np.loadtxt(self.potential_directory + mesh_file)
    
        # Momentum is the first column and weights are the second
        k_array = data[:, 0]  # fm^-1
        k_weights = data[:, 1]  # fm^-1

        return k_array, k_weights
    
    def load_potential(
            self, method='initial', generator='Wegner', lamb=0.0,
            lambda_bd=0.0, k_magnus=0, ds=1e-5):
        """
        Loads the potential. This function is capable of loading SRG- or
        Magnus-evolved potentials as well.
    
        Parameters
        ----------
        method : str, optional
            The evolution method 'srg' or 'magnus'. Choose 'initial' if you
            want the unevolved potential.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float, optional
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        k_magnus : int, optional
            Number of terms to include in Magnus sum (that is,
            dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
        ds : float, optional
            Step-size in the flow parameter s (for Magnus only) [fm^4].
        
        Returns
        -------
        V_matrix : 2-D ndarray
            Potential matrix [fm].
        
        """
        
        # Initialize file name
        file_name = f'vnn_{self.channel}_kvnn_{self.kvnn_string}'
        
        # Name of potential file (split on cases of initial, SRG, or Magnus)
        if method == 'initial':
        
            file_name += '_lam12.0_reg_0_3_0'
        
        # File name is different for SRG- or Magnus-evolved potentials
        else:
    
            if generator == 'Block-diag':
                
                file_name += \
                    f'_{method}_{generator}{lambda_bd:.2f}_lambda{lamb:.1f}'

            else:
            
                # Added this bit to load potentials that were evolved to \lambda
                # values that specify two decimal places (e.g. lamb = 1.35)
                if lamb == round(lamb, 1): # Standard is one decimal place
                
                    file_name += f'_{method}_{generator}_lambda{lamb:.1f}'
                
                else:
            
                    file_name += f'_{method}_{generator}_lambda{lamb:.2f}'
            
            if method == 'magnus':
                
                file_name += f'_k{k_magnus:d}_ds{ds:.1e}'

        file_name += '.out'

        # Load output file
        data = np.loadtxt(self.potential_directory + file_name)

        # Coupled-channel potential?
        if self.coupled_channel_bool:
        
            v11 = np.reshape(data[:, 2], (self.ntot, self.ntot))
            v12 = np.reshape(data[:, 3], (self.ntot, self.ntot))
            v21 = np.reshape(data[:, 4], (self.ntot, self.ntot))
            v22 = np.reshape(data[:, 5], (self.ntot, self.ntot))
            V_matrix = tl.build_coupled_channel_matrix(v11, v12, v21, v22)
        
        else:
        
            V_matrix = np.reshape(data[:, 2], (self.ntot, self.ntot))

        return V_matrix
    
    def load_kinetic_energy(self):
        """
        Loads relative kinetic energy.
        
        Returns
        -------
        T_matrix : 2-D ndarray
            Relative kinetic energy matrix [MeV].
        
        """
    
        # Load momentum array [fm^-1]
        k_array, _ = self.load_mesh()
    
        # Matrix of (h-bar*k)^2 / M along diagonal [MeV]
        T_matrix = np.diag(k_array**2) * self.hbar_sq_over_m
    
        # Coupled-channel operator?
        if self.coupled_channel_bool:
        
            # Matrix of zeros (n x n)
            zeros = np.zeros((self.ntot, self.ntot))
        
            # Build coupled-channel T_rel matrix
            T_matrix = tl.build_coupled_channel_matrix(T_matrix, zeros, zeros,
                                                       T_matrix)
        
        return T_matrix
    
    def load_hamiltonian(
            self, method='initial', generator='Wegner', lamb=0.0,
            lambda_bd=0.0, k_magnus=0, ds=1e-5):
        """
        Loads the NN Hamiltonian. This function is capable of loading SRG-
        or Magnus-evolved Hamiltonians as well.
    
        Parameters
        ----------
        method : str, optional
            The evolution method 'srg' or 'magnus'. Choose 'initial' if you
            want the unevolved potential.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float, optional
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        k_magnus : int, optional
            Number of terms to include in Magnus sum (that is,
            dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
        ds : float, optional
            Step-size in the flow parameter s (for Magnus only) [fm^4].
        
        Returns
        -------
        H_matrix : 2-D ndarray
            Hamiltonian matrix [MeV].
        
        """
        
        # Load relative kinetic energy [MeV]
        T_matrix = self.load_kinetic_energy()
    
        # Load potential [fm]
        V_matrix_fm = self.load_potential(method, generator, lamb, lambda_bd,
                                          k_magnus, ds)

        # Convert to units MeV
        V_matrix = self.convert_V_to_MeV(V_matrix_fm)
    
        # Add to obtain Hamiltonian [MeV]
        H = T_matrix + V_matrix
    
        return H
    
    def save_potential(
            self, V_matrix, method='srg', generator='Wegner', lamb=0.0,
            lambda_bd=0.0, k_magnus=0, ds=1e-5):
        """
        Saves the potential. This function is used for saving SRG- or Magnus-
        evolved potentials.
    
        Parameters
        ----------
        V_matrix : 2-D ndarray
            Potential matrix [fm].
        method : str, optional
            The evolution method 'srg' or 'magnus'. Choose 'initial' if you
            want the unevolved potential.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float, optional
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        k_magnus : int, optional
            Number of terms to include in Magnus sum (that is,
            dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
        ds : float, optional
            Step-size in the flow parameter s (for Magnus only) [fm^4].
        
        """
        
        # Initialize file name
        file_name = f'vnn_{self.channel}_kvnn_{self.kvnn_string}'
        
        # Name of potential file (split on cases of initial, SRG, or Magnus)
        if method == 'initial':
        
            file_name += 'lam12.0_reg_0_3_0.out'
        
        # File name is different for SRG- or Magnus-evolved potentials
        else:
    
            if generator == 'Block-diag':
                
                file_name += \
                    f'_{method}_{generator}{lambda_bd:.2f}_lambda{lamb:.1f}'

            else:
            
                # Added this bit to load potentials that were evolved to \lambda
                # values that specify two decimal places (e.g. lamb = 1.35)
                if lamb == round(lamb, 1): # Standard is one decimal place
                
                    file_name += f'_{method}_{generator}_lambda{lamb:.1f}'
                
                else:
            
                    file_name += f'_{method}_{generator}_lambda{lamb:.2f}'
            
            if method == 'magnus':
                
                file_name += f'_k{k_magnus:d}_ds{ds:.1e}'

        file_name += '.out'
        
        # Get momenta and weights
        k_array, k_weights = self.load_mesh()
        
        f = open(self.potential_directory + file_name, 'w')
    
        # Coupled-channel potentials are written differently - write each sub-
        # block as a column
        if self.coupled_channel_bool:
        
            header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format(
                'k', 'kp', 'V11', 'V12', 'V21', 'V22'
            )
            f.write('#' + header + '\n')
        
            for i, (ik, iw) in enumerate(zip(k_array, k_weights)):
                for j, (jk, jw) in enumerate(zip(k_array, k_weights)):

                    # Divide out the integration measure 2/\pi dk k^2
                    factor = np.pi / (2.0*ik*jk*np.sqrt(iw*jw))
                    v11 = V_matrix[i, j] * factor
                    v12 = V_matrix[i, j+self.ntot] * factor
                    v21 = V_matrix[i+self.ntot, j] * factor
                    v22 = V_matrix[i+self.ntot, j+self.ntot] * factor
                
                    line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
                        ik, jk, v11, v12, v21, v22
                    )
                
                    f.write(line + '\n')
                
        else:
        
            header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'V')
            f.write('#' + header + '\n')
        
            for i, (ik, iw) in enumerate( zip(k_array, k_weights) ):
                for j, (jk, jw) in enumerate( zip(k_array, k_weights) ):

                    # Divide out the integration measure 2/\pi dk k^2
                    factor = np.pi / (2.0*ik*jk*np.sqrt(iw*jw))
                
                    v = V_matrix[i, j] * factor
                
                    line = '{:^15.6f}{:^15.6f}{:^23e}'.format(ik, jk, v)
                
                    f.write(line+'\n')
                
        f.close()
    
    def convert_V_to_MeV(self, V_matrix_fm):
        """
        Converts the potential from units fm to MeV with respect to the
        momentum array and weights of the integration mesh. That is, the
        returned potential includes integration factors along with dividing
        out M / h-bar^2.
        
        Parameters
        ----------
        V_matrix_fm : 2-D ndarray
            Potential matrix [fm].

        Returns
        -------
        V_matrix_MeV : 2-D ndarray
            Potential matrix with integration measure [MeV].

        """
        
        k_array, k_weights = self.load_mesh()
    
        # If coupled channel, double the length of the mesh arrays
        if self.coupled_channel_bool:
            
            k_array = np.concatenate((k_array, k_array))
            k_weights = np.concatenate((k_weights, k_weights))
            
        # Build grids of momentum and weights factor
        col, row = np.meshgrid(k_array * np.sqrt(2/np.pi*k_weights),
                               k_array * np.sqrt(2/np.pi*k_weights))
    
        # Multiply the potential by 2/pi * row * col -> fm^-2 conversion and 
        # multiplying by hbar_sq_over_m gives the MeV conversion
        V_matrix_MeV = V_matrix_fm * row * col * self.hbar_sq_over_m
    
        return V_matrix_MeV
    
    def convert_V_to_fm(self, V_matrix_MeV):
        """
        Converts the potential from units MeV to fm with respect to the
        momentum array and weights of the integration mesh. This is V(k,k') 
        in the following equation:
            k^2 \psi(k) + \int dk k'^2 V(k,k') \psi(k') = E \psi(k).
        
        Parameters
        ----------
        V_matrix_MeV : 2-D ndarray
            Potential matrix with integration measure [MeV].

        Returns
        -------
        V_matrix_fm : 2-D ndarray
            Potential matrix without integration measure [fm].

        """
    
        k_array, k_weights = self.load_mesh()
    
        # If coupled channel, double the length of the mesh arrays
        if self.coupled_channel_bool:
            
            k_array = np.concatenate((k_array, k_array))
            k_weights = np.concatenate((k_weights, k_weights))
            
        # Build grids of momentum and weights factor
        col, row = np.meshgrid(k_array * np.sqrt(2/np.pi*k_weights),
                               k_array * np.sqrt(2/np.pi*k_weights))
            
        # Build grids of momentum and weights factor
        col, row = np.meshgrid(k_array * np.sqrt(2/np.pi*k_weights),
                               k_array * np.sqrt(2/np.pi*k_weights))
    
        # Dividing the potential by 2/pi * row * col gives MeV fm^3 conversion
        # and dividing by hbar_sq_over_m gives the fm conversion
        V_matrix_fm = V_matrix_MeV / (2/np.pi*row*col*self.hbar_sq_over_m)
    
        return V_matrix_fm
    
    def run_srg(
            self, generator, lambda_array, lambda_bd_array=np.empty(0),
            save=True):
        """
        SRG evolves the specified Hamiltonian to several values of \lambda,
        and possibly \Lambda_BD, and has the option to save the evolved
        potentials.

        Parameters
        ----------
        generator : str
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lambda_array : 1-D ndarray
            SRG evolution parameters \lambda [fm^-1].
        lambda_bd_array : 1-D ndarray, optional
            \Lambda_BD values for block-diagonal generator [fm^-1].
        save : bool, optional
            If true, saves the evolved potentials.

        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each \lambda value (and possibly an additional
            key for \Lambda_BD). E.g., d[1.5] returns the evolved Hamiltonian
            at \lambda = 1.5 fm^-1.

        """
    
        # Initial value of \lambda [fm^-1]
        # Technically, this value should be infinity, but we can take 20 fm^-1
        # which is reasonably large
        lambda_initial = 20.0

        # Initialize SRG class (we're feeding in the Potential object as the
        # first argument to the SRG class)
        srg = SRG(self, generator)

        # Time the evolution and return dictionary d
        t0 = time.time() # Start time
        d = srg.evolve(lambda_initial, lambda_array, lambda_bd_array)
        t1 = time.time() # End time
    
        # Print details
        mins = round((t1-t0)/60.0, 4)  # Minutes elapsed evolving H(s)
        print('_'*85)
        print(
            f'Done evolving to final \lambda = {lambda_array[-1]:.2f} fm^-1 after {mins:.4f} minutes'
        )
        print('_'*85)
        print('\nSpecifications:\n')
        print(f'kvnn = {self.kvnn:d}, channel = {self.channel}')
        print(f'kmax = {self.kmax:.1f}, kmid = {self.kmid:.1f}, ntot = {self.ntot:d}')
        print(f'method = srg, generator = {generator}')
        if generator == 'Block-diag':
            print(f'Final \Lambda_BD = {lambda_bd_array[-1].2f} fm^-1')
    
        # Save evolved potentials
        if save:
            
            # Get relative kinetic energy and convert to [fm^-2]
            T_matrix = self.load_kinetic_energy() / self.hbar_sq_over_m

            if generator == 'Block-diag':
                
                # Additionally loop over \Lambda_BD
                for lambda_bd in lambda_bd_array:
                    for lamb in lambda_array:

                        # Scattering units here [fm^-2]
                        H_matrix = d[lambda_bd][lamb]
                    
                        # Subtract off kinetic energy [fm^-2]
                        V_matrix = H_matrix - T_matrix
                    
                        # Save evolved potential
                        self.save_potential(V_matrix, 'srg', generator, lamb,
                                            lambda_bd)
                
            # Only need to loop over \lambda for band-diagonal generators
            else:
            
                for lamb in lambda_array:

                    # Scattering units here [fm^-2]
                    H_matrix = d[lamb]
                    
                    # Subtract off kinetic energy [fm^-2]
                    V_matrix = H_matrix - T_matrix
                    
                    # Save evolved potential
                    self.save_potential(V_matrix, 'srg', generator, lamb)
                
        return d
    
    def load_omega(self):
        return None
    
    def save_omega(self):
        return None
    
    def run_magnus(self):
        return None


# def load_omega(kvnn, channel, kmax=0.0, kmid=0.0, ntot=0, generator='Wegner', 
#                lamb=1.5, lambda_bd=0.00, k_magnus=6, ds=1e-5):
#     """
#     Loads a Magnus evolved omega matrix.
    
#     Parameters
#     ----------
#     kvnn : int
#         This number specifies the potential.
#     channel : str
#         The partial wave channel (e.g. '1S0').
#     kmax : float, optional
#         Maximum value in the momentum mesh [fm^-1].
#     kmid : float, optional
#         Mid-point value in the momentum mesh [fm^-1].
#     ntot : int, optional
#         Number of momentum points in mesh.
#     generator : str, optional
#         SRG generator 'Wegner', 'T', or 'Block-diag'.
#     lamb : float, optional
#         Evolution parameter lambda [fm^-1].
#     lambda_bd : float, optional
#         Lambda value for block-diagonal decoupling [fm^-1].
#     k_magnus : int, optional
#         Number of terms to include in Magnus sum (that is,
#         dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
#     ds : float, optional
#         Step-size in the flow parameter s (for Magnus only) [fm^4].
        
#     Returns
#     -------
#     O : 2-D ndarray
#         Magnus evolved omega matrix.
        
#     """

#     # Convert kvnn to string
#     if kvnn < 10:
#         kvnn_string = '0'+str(kvnn)
#     else:
#         kvnn_string = str(kvnn)
    
#     # Get potential directory
#     potential_directory = 'potentials/vsrg_macos/vsrg_kvnn_' + \
#                           '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
#                           (kvnn_string, kmax, kmid, ntot)
    
#     # Name of omega file
#     omega_file = 'omega_%s_kvnn_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
#                  (channel, kvnn_string, generator, lamb, k_magnus, ds)

#     # Load output file
#     data = np.loadtxt(potential_directory + '/' + omega_file)

#     # Coupled channel potential
#     if coupled_channel(channel):
        
#         o11 = np.reshape( data[: ,2], (ntot, ntot) )
#         o12 = np.reshape( data[:, 3], (ntot, ntot) )
#         o21 = np.reshape( data[:, 4], (ntot, ntot) )
#         o22 = np.reshape( data[:, 5], (ntot, ntot) )
#         O = np.vstack( ( np.hstack( (o11, o12) ), np.hstack( (o21, o22) ) ) )
        
#     else:
        
#         O = np.reshape( data[:, 2], (ntot, ntot) )

#     return O


# def save_potential(k_array, k_weights, V, kvnn, channel, kmax=0.0, kmid=0.0, 
#                    ntot=0, method='srg', generator='Wegner', lamb=1.5, 
#                    lambda_bd=0.00, k_magnus=6, ds=1e-5):
#     """
#     Saves an SRG or Magnus evolved potential.
    
#     Parameters
#     ----------
#     k_array : 1-D ndarray
#         Momentum array [fm^-1].
#     k_weights: 1-D ndarray
#         Momentum weights [fm^-1].
#     V : 2-D ndarray
#         Potential matrix [fm^-2].
#     kvnn : int
#         This number specifies the potential.
#     channel : str
#         The partial wave channel (e.g. '1S0').
#     kmax : float, optional
#         Maximum value in the momentum mesh [fm^-1].
#     kmid : float, optional
#         Mid-point value in the momentum mesh [fm^-1].
#     ntot : int, optional
#         Number of momentum points in mesh.
#     method : str, optional
#         The evolution method 'srg' or 'magnus'. Choose 'initial' if unevolved.
#     generator : str, optional
#         SRG generator 'Wegner', 'T', or 'Block-diag'.
#     lamb : float, optional
#         Evolution parameter lambda [fm^-1].
#     lambda_bd : float, optional
#         Lambda value for block-diagonal decoupling [fm^-1].
#     k_magnus : int, optional
#         Number of terms to include in Magnus sum (that is,
#         dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
#     ds : float, optional
#         Step-size in the flow parameter s (for Magnus only) [fm^4].
        
#     """

#     # Get current working directory
#     cwd = getcwd()
    
#     # Go to directory of specified potential
#     # Convert kvnn to string
#     if kvnn < 10:
#         kvnn_string = '0'+str(kvnn)
#     else:
#         kvnn_string = str(kvnn)
        
#     potential_directory = 'potentials/vsrg_macos/vsrg_kvnn_' + \
#                           '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
#                           (kvnn_string, kmax, kmid, ntot)
    
#     chdir(potential_directory)
    
#     # Name potential file and save
#     if method == 'srg': 
    
#         if generator == 'Block-diag':
            
#             vnn_file = 'vnn_%s_kvnn_%s_%s_%s%.2f_lambda%.1f.out' % (channel, 
#                         kvnn_string, method, generator, lambda_bd, lamb)
            
#         else:
            
#             # Added this bit to load potentials that were evolved to \lambda
#             # values that specify two decimal places (e.g. lamb = 1.35)
#             if lamb == round(lamb, 1): # Standard is one decimal place
                
#                 vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f.out' % (channel,
#                             kvnn_string, method, generator, lamb)
                
#             else:
            
#                 vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.2f.out' % (channel,
#                             kvnn_string, method, generator, lamb)
            
#     elif method == 'magnus': 

#         vnn_file = 'vnn_%s_kvnn_%s_%s_%s_lambda%.1f_k%d_ds%.1e.out' % (channel, 
#                     kvnn_string, method, generator, lamb, k_magnus, ds)
    
#     else:
        
#         vnn_file = 'vnn_%s_kvnn_%s_lam12.0_reg_0_3_0.out' % (channel, 
#                                                              kvnn_string)
    
#     f = open(vnn_file,'w')
    
#     # Coupled-channel potentials are written differently - write each sub-
#     # block as a column
#     if coupled_channel(channel):
        
#         header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format('k',
#                  'kp', 'V11', 'V12', 'V21', 'V22')
#         f.write('#' + header + '\n')
        
#         for i, (ik, iw) in enumerate( zip(k_array, k_weights) ):
#             for j, (jk, jw) in enumerate( zip(k_array, k_weights) ):

#                 # Divide out the integration measure 2/\pi dk k^2
#                 factor = np.pi / ( 2.0 * ik * jk * np.sqrt(iw * jw) )
#                 v11 = V[i, j] * factor
#                 v12 = V[i, j+ntot] * factor
#                 v21 = V[i+ntot, j] * factor
#                 v22 = V[i+ntot, j+ntot] * factor
                
#                 line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
#                         ik, jk, v11, v12, v21, v22)
                
#                 f.write(line + '\n')
                
#     else:
        
#         header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'V')
#         f.write('#' + header + '\n')
        
#         for i, (ik, iw) in enumerate( zip(k_array, k_weights) ):
#             for j, (jk, jw) in enumerate( zip(k_array, k_weights) ):

#                 # Divide out the integration measure 2/\pi dk k^2
#                 factor = np.pi / ( 2.0 * ik * jk * np.sqrt(iw * jw) )
                
#                 v = V[i,j] * factor
                
#                 line = '{:^15.6f}{:^15.6f}{:^23e}'.format(ik, jk, v)
                
#                 f.write(line+'\n')
                
#     f.close()

#     chdir(cwd)
    
    
# def save_omega(k_array, O, kvnn, channel, kmax=0.0, kmid=0.0, ntot=0, 
#                generator='Wegner', lamb=1.5, lambda_bd=0.00, k_magnus=6, 
#                ds=1e-5):
#     """
#     Saves a Magnus evolved omega matrix.
    
#     Parameters
#     ----------
#     k_array : 1-D ndarray
#         Momentum array [fm^-1].
#     O : 2-D ndarray
#         Magnus evolved omega matrix.
#     kvnn : int
#         This number specifies the potential.
#     channel : str
#         The partial wave channel (e.g. '1S0').
#     kmax : float, optional
#         Maximum value in the momentum mesh [fm^-1].
#     kmid : float, optional
#         Mid-point value in the momentum mesh [fm^-1].
#     ntot : int, optional
#         Number of momentum points in mesh.
#     generator : str, optional
#         SRG generator 'Wegner', 'T', or 'Block-diag'.
#     lamb : float, optional
#         Evolution parameter lambda [fm^-1].
#     lambda_bd : float, optional
#         Lambda value for block-diagonal decoupling [fm^-1].
#     k_magnus : int, optional
#         Number of terms to include in Magnus sum (that is,
#         dOmega / ds ~ \sum_0^k_magnus ... for Magnus only).
#     ds : float, optional
#         Step-size in the flow parameter s (for Magnus only) [fm^4].
        
#     """

#     # Get current working directory
#     cwd = getcwd()
    
#     # Go to directory of specified potential
#     # Convert kvnn to string
#     if kvnn < 10:
#         kvnn_string = '0'+str(kvnn)
#     else:
#         kvnn_string = str(kvnn)
        
#     potential_directory = 'potentials/vsrg_macos/vsrg_kvnn_' + \
#                           '%s_lam12.0_kmax%d_kmid%d_ntot%d' % \
#                           (kvnn_string, kmax, kmid, ntot)
    
#     chdir(potential_directory)
    
#     # Name omega file and save
#     omega_file = 'omega_%s_kvnn_%s_%s_lambda%.1f_k%d_ds%.1e.out' % \
#                  (channel, kvnn_string, generator, lamb, k_magnus, ds)
    
#     f = open(omega_file,'w')
    
#     if coupled_channel(channel):
        
#         header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format('k', 'kp', 
#                   'O11', 'O12', 'O21', 'O22')
#         f.write('#' + header + '\n')
        
#         for i, ik in enumerate(k_array):
#             for j, jk in enumerate(k_array):

#                 o11 = O[i, j]
#                 o12 = O[i, j+ntot]
#                 o21 = O[i+ntot, j]
#                 o22 = O[i+ntot, j+ntot]
                
#                 line = '{:^15.6f}{:^15.6f}{:^23e}{:^23e}{:^23e}{:^23e}'.format(
#                         ik, jk, o11, o12, o21, o22)
                
#                 f.write(line + '\n')
                
#     else:
        
#         header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'O')
#         f.write('#' + header + '\n')
        
#         for i, ik in enumerate(k_array):
#             for j, jk in enumerate(k_array):
                
#                 o = O[i,j]
                
#                 line = '{:^15.6f}{:^15.6f}{:^23e}'.format(ik, jk, o)
                
#                 f.write(line + '\n')
                
#     f.close()

#     chdir(cwd)