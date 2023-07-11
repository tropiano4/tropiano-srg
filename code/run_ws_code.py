#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run the Woods-Saxon Fortran code for some nucleus.

WORK IN PROGRESS!
"""

# Python imports
import numpy as np
import shutil

# Imports from scripts
from scripts.woodsaxon import ws


def main(nucleus_name, Z, N, n_max, l_max, rmax, ntab):
    """Run Woods-Saxon code to generate data."""
        
    # Total number of nucleons
    A = Z + N
        
    # Type of orbitals: 1 - nucleons with no Coulomb
    #                   2 - distinguish protons and neutrons
    ntau = 2
        
    # Orbitals to consider (note, we track 2*j not j)
    norb, lorb, jorb = [], [], []
    for n in range(1, n_max+1):
        for l in range(0, l_max+1):
            norb.append(n)
            lorb.append(l)
            jorb.append(int(2*(l+1/2)))
            if int(2*(l-1/2)) > 0:  # Don't append negative j
                norb.append(n)
                lorb.append(l)
                jorb.append(int(2*(l-1/2)))
    nrad = len(jorb)
    orbws = np.zeros(shape=(2,nrad,ntab), order='F')
    
    # Divide orbital by r? -> get R(r); false: get u(r)=r R(r)
    rdiv = False
    dens = True
    
    # Set parameters of the Woods-Saxon potential
    prm = np.zeros(shape=(2,9), order='F')
    
    # Starting with vws (p & n)
    if nucleus_name == 'He4':
        prm[:,0] = 76.8412
    elif nucleus_name == 'Be9':
        prm[:,0] = 66.6397
    elif nucleus_name == 'C12':
        prm[:,0] = 60.1478
    elif nucleus_name == 'O16':
        prm[:,0] = 58.0611
    elif nucleus_name == 'Ca40':
        prm[:,0] = 54.3051
    elif nucleus_name == 'Ca48':
        prm[0,0] = 59.4522
        prm[1,0] = 46.9322
    elif nucleus_name == 'Fe56':
        prm[0,0] = 55.9744
        prm[1,0] = 50.0125
    elif nucleus_name == 'Pb208':
        prm[0,0] = 59.3452
        prm[1,0] = 44.899
    
    # Other parameters of the Woods-Saxon
    prm[:,1] = 1.275
    prm[:,2] = 0.7
    prm[:,3] = 0.
    prm[:,4] = 1.
    prm[:,5] = 36
    prm[:,6] = 1.32
    prm[:,7] = 0.7
    prm[:,8] = 1.275
        
    # Print summary, potentials, and densities
    prnt = True
    prntorb = True

    # Run Fortran subroutine
    ws(ntau, A, Z, rmax, orbws, norb, lorb, jorb, prm, rdiv, prnt, prntorb,
       dens)
    
    
if __name__ == '__main__':
    
    # He4
    nucleus_name, Z, N = 'He4', 2, 2
    
    # Woods-Saxon directory
    woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}/"
                
    # Order single-particle states with lowest energy first
    self.order_sp_states(Z, N)
        
    # Organize wave functions in dictionary with the file name as the key
    self.sp_wfs = {}
        
    for sp_state in self.sp_states:
            
        # Wave functions are independent of m_j, so fix m_j=j
        if sp_state.m_j == sp_state.j:
                
            file_name = self.get_orbital_file_name(sp_state)
                
            if run_woods_saxon:
                    
                shutil.move(file_name, self.woods_saxon_directory + file_name)
                    
                data = np.loadtxt(self.woods_saxon_directory + file_name)

    
    main()
    
    # Move output files to relevant directory
    shutil.move("ws_log", self.woods_saxon_directory + "ws_log")
    shutil.move("ws_pot", self.woods_saxon_directory + "ws_pot")
    shutil.move("ws_rho", self.woods_saxon_directory + "ws_rho")