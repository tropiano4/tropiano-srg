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

Last update: March 25, 2022

"""

# To-do: Change file naming convention to be more consistent across scripts.
#        E.g., use convert_number_to_string() for labeling \lambda, etc.
# To-do: Use np.savetxt() to save potentials? Definitely a simpler way.
# To-do: Should probably incorporate dU/ds = \eta U method?
# To-do: Should probably incorporate Magnus_split method?
# To-do: Combine run_srg and run_magnus into a one method.

# Python imports
import numpy as np
import time

# Imports from A.T. codes
from modules.magnus import Magnus
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
        kmax_int = int(kmax)
        kmid_int = int(kmid)
        potential_directory = ('../potentials/vsrg_macos/vsrg_kvnn'
                               f'_{kvnn_string}_lam12.0_kmax{kmax_int:d}'
                               f'_kmid{kmid_int:d}_ntot{ntot:d}/')
                            
        self.potential_directory = potential_directory
        
    def get_file_name(
            self, file_type, generator='', lamb=0.0, lambda_bd=0.0, k_max=0,
            ds=0.0):
        """
        Handles file names for momentum meshes, initial potentials, and SRG-
        or Magnus-evolved potentials.

        Parameters
        ----------
        file_type : str
            Mesh, initial potential, SRG-, or Magnus-evolved potential file:
                'mesh', 'initial', 'srg', or 'magnus'.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float, optional
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max (for Magnus evolution).
        ds : float, optional
            Step-size in the flow parameter s [fm^4] (for Magnus evolution).

        Returns
        -------
        file_name : str
            Name of the file.

        """
        
        if file_type == 'mesh':
            
            file_name = (f'vsrg_{self.channel}_kvnn_{self.kvnn_string}'
                         '_lam12.0_reg_0_3_0_mesh.out')
            
        elif file_type == 'initial':
            
            file_name = (f'vnn_{self.channel}_kvnn_{self.kvnn_string}'
                         '_lam12.0_reg_0_3_0.out')
            
        elif file_type == 'srg' or file_type == 'magnus':
            
            file_name = (f'vnn_{self.channel}_kvnn_{self.kvnn_string}'
                         f'_{file_type}_{generator}')
            
            # Get \lambda with correct number of decimals
            if lamb == round(lamb, 1):
                lamb_str = str(round(lamb, 1))
            else:
                lamb_str = str(round(lamb, 2))
            
            # Add \Lambda_BD to name for block-diagonal generator
            if generator == 'Block-diag':  
                file_name += f'{lambda_bd:.2f}_lambda{lamb_str}'
            else:
                file_name += f'_lambda{lamb_str}'
            
            # Add k_max and ds for Magnus evolution
            if file_type == 'magnus':
                file_name += f'_k{k_max:d}_ds{ds:.1e}'

            file_name += '.out'
            
        else:
            
            raise RuntimeError(
                "Invalid file type. Please specify one of the following:\n"
                "\t'mesh', 'initial', 'srg', or 'magnus'.")
        
        return file_name
        
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
        mesh_file = self.get_file_name('mesh')
        data = np.loadtxt(self.potential_directory + mesh_file)
    
        # Momentum is the first column and weights are the second
        k_array = data[:, 0]  # fm^-1
        k_weights = data[:, 1]  # fm^-1

        return k_array, k_weights
    
    def load_potential(
            self, method='initial', generator='', lamb=0.0,
            lambda_bd=0.0, k_max=0, ds=0.0):
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
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max (for Magnus evolution).
        ds : float, optional
            Step-size in the flow parameter s [fm^4] (for Magnus evolution).
        
        Returns
        -------
        V_matrix : 2-D ndarray
            Potential matrix [fm].
        
        """
        
        # Get file name for the potential matrix elements
        file_name = self.get_file_name(method, generator, lamb, lambda_bd,
                                       k_max, ds)

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
            self, method='initial', generator='', lamb=0.0, lambda_bd=0.0,
            k_max=0, ds=0.0):
        """
        Loads the NN Hamiltonian. This function is capable of loading SRG-
        or Magnus-evolved Hamiltonians as well.
    
        Parameters
        ----------
        method : str, optional
            The evolution method 'srg' or 'magnus'. Choose 'initial' if you
            want the unevolved Hamiltonian.
        generator : str, optional
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float, optional
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max (for Magnus evolution).
        ds : float, optional
            Step-size in the flow parameter s [fm^4] (for Magnus evolution).
        
        Returns
        -------
        H_matrix : 2-D ndarray
            Hamiltonian matrix [MeV].
        
        """
        
        # Load relative kinetic energy [MeV]
        T_matrix = self.load_kinetic_energy()
    
        # Load potential [fm]
        V_matrix_fm = self.load_potential(method, generator, lamb, lambda_bd,
                                          k_max, ds)

        # Convert potential to units MeV
        V_matrix = self.convert_V_to_MeV(V_matrix_fm)
    
        H_matrix = T_matrix + V_matrix
    
        return H_matrix
    
    def save_potential(
            self, V_matrix, method, generator, lamb, lambda_bd=0.0, k_max=0,
            ds=0.0):
        """
        Saves the potential. This function is used for saving SRG- or Magnus-
        evolved potentials.
    
        Parameters
        ----------
        V_matrix : 2-D ndarray
            Potential matrix [fm^-2].
        method : str
            The evolution method 'srg' or 'magnus'.
        generator : str
            SRG generator 'Wegner', 'T', or 'Block-diag'.
        lamb : float
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max (for Magnus evolution).
        ds : float, optional
            Step-size in the flow parameter s [fm^4] (for Magnus evolution).
        
        """
        
        # Get file name for the potential matrix elements
        file_name = self.get_file_name(method, generator, lamb, lambda_bd,
                                       k_max, ds)

        # Get momenta and weights
        k_array, k_weights = self.load_mesh()
        
        f = open(self.potential_directory + file_name, 'w')
    
        # Write each sub-block as a column for coupled-channel potentials
        if self.coupled_channel_bool:
        
            header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format(
                'k', 'kp', 'V11', 'V12', 'V21', 'V22'
            )
            f.write('#' + header + '\n')
        
            for i, (ik, iw) in enumerate(zip(k_array, k_weights)):
                for j, (jk, jw) in enumerate(zip(k_array, k_weights)):

                    # Divide out the integration measure 2/\pi dk k^2
                    # Units go from fm^-2 to fm
                    factor = np.pi / (2.0*ik*jk*np.sqrt(iw*jw))
                    
                    v11 = V_matrix[i, j] * factor
                    v12 = V_matrix[i, j+self.ntot] * factor
                    v21 = V_matrix[i+self.ntot, j] * factor
                    v22 = V_matrix[i+self.ntot, j+self.ntot] * factor
                
                    line = ('{:^15.6f}{:^15.6f}'.format(ik, jk)
                            + '{:^23e}{:^23e}'.format(v11, v12)
                            + '{:^23e}{:^23e}'.format(v21, v22))
                
                    f.write(line + '\n')
                
        else:
        
            header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'V')
            f.write('#' + header + '\n')
        
            for i, (ik, iw) in enumerate( zip(k_array, k_weights) ):
                for j, (jk, jw) in enumerate( zip(k_array, k_weights) ):

                    # Divide out the integration measure 2/\pi dk k^2
                    # Units go from fm^-2 to fm
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
        factor_array = k_array * np.sqrt(2/np.pi*k_weights)
    
        # If coupled channel, double the length of the mesh arrays
        if self.coupled_channel_bool:
            factor_array = np.concatenate((factor_array, factor_array))
            
        # Build meshgrids of momentum and weights factor
        col, row = np.meshgrid(factor_array, factor_array)
    
        # Multiply the potential by the integration measure giving the fm^-2 
        # conversion and multiplying by hbar_sq_over_m gives the MeV conversion
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
        factor_array = k_array * np.sqrt(2/np.pi*k_weights)
    
        # If coupled channel, double the length of the mesh arrays
        if self.coupled_channel_bool:
            factor_array = np.concatenate((factor_array, factor_array))
            
        # Build meshgrids of momentum and weights factor
        col, row = np.meshgrid(factor_array, factor_array)
    
        # Dividing the potential by the integration measure gives the MeV fm^3 
        # conversion and dividing by hbar_sq_over_m gives the fm conversion
        V_matrix_fm = V_matrix_MeV / (row*col*self.hbar_sq_over_m)
    
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
        # Technically, this value should be infinity, but we can take 
        # 20 fm^-1 which is reasonably large
        lambda_initial = 20.0

        # Initialize SRG class (we're feeding in the Potential object as the
        # first argument to the SRG class)
        srg = SRG(self, generator)

        # Time the evolution and return dictionary d
        t0 = time.time() # Start time
        d = srg.srg_evolve(lambda_initial, lambda_array, lambda_bd_array)
        t1 = time.time() # End time
    
        # Print details
        mins = round((t1-t0)/60.0, 4)  # Minutes elapsed evolving H(s)
        print('_'*85)
        lamb_str = tl.convert_number_to_string(lambda_array[-1])
        print(f'Done evolving to final \lambda = {lamb_str} fm^-1 after'
              f' {mins:.4f} minutes.')
        print('_'*85)
        print('\nSpecifications:\n')
        print(f'kvnn = {self.kvnn:d}, channel = {self.channel}')
        print(f'kmax = {self.kmax:.1f}, kmid = {self.kmid:.1f}, '
              f'ntot = {self.ntot:d}')
        print(f'method = SRG, generator = {generator}')
        if generator == 'Block-diag':
            lambda_bd_str = tl.convert_number_to_string(lambda_bd_array[-1])
            print(f'Final \Lambda_BD = {lambda_bd_str} fm^-1')
    
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
                    
                        # Save evolved potential in units [fm]
                        self.save_potential(V_matrix, 'srg', generator, lamb,
                                            lambda_bd)
                
            # Only need to loop over \lambda for band-diagonal generators
            else:
            
                for lamb in lambda_array:

                    # Scattering units here [fm^-2]
                    H_matrix = d[lamb]
                    
                    # Subtract off kinetic energy [fm^-2]
                    V_matrix = H_matrix - T_matrix
                    
                    # Save evolved potential in units [fm]
                    self.save_potential(V_matrix, 'srg', generator, lamb)
                
        return d

    def run_magnus(
            self, generator, lambda_array, lambda_bd_array=np.empty(0),
            k_max=6, ds=1e-5, save=True):
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
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max.
        ds : float, optional
            Step-size in the flow parameter s [fm^4].
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
    
        # Initialize Magnus class (we're feeding in the Potential object as
        # the first argument to the Magnus class)
        magnus = Magnus(self, generator)

        # Time the evolution and return dictionary d
        t0 = time.time() # Start time
        d = magnus.magnus_evolve(lambda_array, lambda_bd_array, k_max, ds)
        t1 = time.time() # End time
    
        # Print details
        mins = round((t1-t0)/60.0, 4)  # Minutes elapsed evolving H(s)
        print('_'*85)
        lamb_str = tl.convert_number_to_string(lambda_array[-1])
        print(f'Done evolving to final \lambda = {lamb_str} fm^-1 after'
              f' {mins:.4f} minutes.')
        print('_'*85)
        print('\nSpecifications:\n')
        print(f'kvnn = {self.kvnn:d}, channel = {self.channel}')
        print(f'kmax = {self.kmax:.1f}, kmid = {self.kmid:.1f}, '
              f'ntot = {self.ntot:d}')
        print(f'method = Magnus, generator = {generator}')
        if generator == 'Block-diag':
            lambda_bd_str = tl.convert_number_to_string(lambda_bd_array[-1])
            print(f'Final \Lambda_BD = {lambda_bd_str} fm^-1')
        print(f'k_max = {k_max:d}, ds = {ds:.1e}')
    
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
                    
                        # Save evolved potential in units [fm]
                        # For large \lambda, Magnus-evolved potentials 
                        # automatically have ds <= 1e-6
                        if lamb >= 10.0 and ds > 1e-6:
                            self.save_potential(V_matrix, 'magnus', generator,
                                                lamb, lambda_bd, k_max, 1e-6)
                        else:
                            self.save_potential(V_matrix, 'magnus', generator,
                                                lamb, lambda_bd, k_max, ds)
                
            # Only need to loop over \lambda for band-diagonal generators
            else:
            
                for lamb in lambda_array:

                    # Scattering units here [fm^-2]
                    H_matrix = d[lamb]
                    
                    # Subtract off kinetic energy [fm^-2]
                    V_matrix = H_matrix - T_matrix
                    
                    # Save evolved potential in units [fm]
                    # For large \lambda, Magnus-evolved potentials 
                    # automatically have ds <= 1e-6
                    if lamb >= 10.0 and ds > 1e-6:
                        self.save_potential(V_matrix, 'magnus', generator,
                                            lamb, lambda_bd, k_max, 1e-6)
                    else:
                        self.save_potential(V_matrix, 'magnus', generator,
                                            lamb, lambda_bd, k_max, ds)
                
        return d