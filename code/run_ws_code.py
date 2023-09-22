#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: single_particle_states.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: July 26, 2023

Run the Woods-Saxon Fortran code for some nucleus. See
https://nucracker.volya.net or arXiv:0709.3525 [nucl-th]
to set parameters of nuclei.

Last update: September 22, 2023

"""

# Python imports
import numpy as np
import shutil

# Imports from scripts
from scripts.woodsaxon import ws
from single_particle_states import SingleParticleState


def get_orbitals():
    """Return lists of orbitals."""

    # Track 2*j not j
    norb, lorb, jorb = [], [], []
    # These must match the nmax and lmax in woodsaxon.f90
    nmax, lmax = 3, 6
    for n in range(1, n_max+1):
        for l in range(0, l_max+1):
            norb.append(n)
            lorb.append(l)
            # j = l + 1/2
            jorb.append(int(2*(l+1/2)))
            # j = l - 1/2 (don't append negative j)
            if int(2*(l-1/2)) > 0:
                norb.append(n)
                lorb.append(l)
                jorb.append(int(2*(l-1/2)))

    return norb, lorb, jorb


def set_ws_parameters(nucleus_name):
    """Set the Woods-Saxon parameters given the nucleus."""

    prm = np.zeros(shape=(2, 9), order='F')
    
    # Central potential strength
    if nucleus_name == 'He4':
        prm[:, 0] = 76.8412
    elif nucleus_name == 'Be9':
        prm[:, 0] = 66.6397
    elif nucleus_name == 'C12':
        prm[:, 0] = 60.1478
    elif nucleus_name == 'O16':
        prm[:, 0] = 58.0611
    elif nucleus_name == 'Ca40':
        prm[:, 0] = 54.3051
    elif nucleus_name == 'Ca48':
        prm[0, 0] = 59.4522  # Proton
        prm[1, 0] = 46.9322  # Neutron
    elif nucleus_name == 'Fe56':
        prm[0, 0] = 55.9744  # Proton
        prm[1, 0] = 50.0125  # Neutron
    elif nucleus_name == 'Sn118':
        prm[0, 0] = 57.7428  # Proton
        prm[1, 0] = 46.9911  # Neutron
    elif nucleus_name == 'Pb208':
        prm[0, 0] = 59.3452  # Proton
        prm[1, 0] = 44.899  # Neutron
    else:
        raise RuntimeError("Don't have this nucleus yet.")
    
    # Other parameters of the Woods-Saxon
    # Central R_0
    # prm[:,1] = 1.275
    prm[0,1] = 1.275  # Proton
    prm[1,1] = 1.347  # Neutron
    # Central surface diffuseness a
    prm[:,2] = 0.7
    # Wine-Bottle: kwb
    prm[:,3] = 0.
    # Wine-Bottle: awb
    prm[:,4] = 1.
    # Spin-orbit strength \lambda
    # prm[:,5] = 36
    prm[0,5] = 36  # Proton
    prm[1,5] = 35  # Neutron
    # Spin-orbit R_0
    # prm[:,6] = 1.32
    prm[0,6] = 1.32  # Proton
    prm[1,6] = 1.31  # Neutron
    # Spin-orbit surface diffuseness a
    prm[:,7] = 0.7
    # Coulomb R_0
    prm[:,8] = 1.275

    return prm


def order_sp_states(directory, Z, N):
    """Get organized lists of all s.p. states and occupied states."""

    sp_states = []
    occ_states = []
    proton_count = 0
    neutron_count = 0
        
    # File with details of the orbitals
    ws_file = directory + "ws_log"
    
    # Order single-particle states using the ws_log file
    with open(ws_file, 'r') as f:
        for line in f:
            unit = line.strip().split()
                
            # Protons
            if len(unit) > 0 and unit[0] == '1':

                j = int(unit[3])/2
                for m_j in np.arange(-j, j+1, 1):
                    sp_state = SingleParticleState(
                        int(unit[1])+1, int(unit[2]), j, m_j, 1/2
                    )  # n, l, j, m_j, m_t
                    
                    sp_states.append(sp_state)
                
                    if proton_count < Z:
                        occ_states.append(sp_state)
                        # Add up filled proton states
                        proton_count += 1
                
            # Neutrons
            elif len(unit) > 0 and unit[0] == '2':

                j = int(unit[3])/2
                for m_j in np.arange(-j, j+1, 1):
                    sp_state = SingleParticleState(
                        int(unit[1])+1, int(unit[2]), j, m_j, -1/2
                    )  # n, l, j, m_j, m_t
                    
                    sp_states.append(sp_state)
                
                    if neutron_count < N:
                        occ_states.append(sp_state)
                        # Add up filled neutron states
                        neutron_count += 1
                        
    return sp_states, occ_states


def get_orbital_file_name(sp_state):
    """Returns the file name of the orbital."""
            
    # Proton
    if sp_state.m_t == 1/2:
        file_name = (f"p.n{int(sp_state.n-1)}.l{int(sp_state.l)}"
                     f".j{int(2*sp_state.j)}.orb")
    # Neutron
    elif sp_state.m_t == -1/2:
        file_name = (f"n.n{int(sp_state.n-1)}.l{int(sp_state.l)}"
                     f".j{int(2*sp_state.j)}.orb")
            
    return file_name


def main(nucleus_name, Z, N, rmax=40, ntab=2000):
    """Run Woods-Saxon code to generate orbital files."""
        
    # Total number of nucleons
    A = Z + N
        
    # Type of orbitals: 1 - nucleons with no Coulomb
    #                   2 - distinguish protons and neutrons
    ntau = 2
        
    norb, lorb, jorb = get_orbitals()
    nrad = len(jorb)
    orbws = np.zeros(shape=(2, nrad, ntab), order='F')
    
    # Divide orbital by r? -> get R(r); false: get u(r)=r*R(r)
    rdiv = False
    dens = True
    
    # Get the Woods-Saxon parameters for the input nucleus
    prm = set_ws_parameters(nucleus_name)
  
    # Print summary, potentials, and densities
    prnt = True
    prntorb = True

    # Run Fortran subroutine
    ws(ntau, A, Z, rmax, orbws, norb, lorb, jorb, prm, rdiv, prnt, prntorb,
       dens)

    # Woods-Saxon directory
    woods_saxon_directory = f"../data/woods_saxon/{nucleus_name}/"
    
    # Move output files to relevant directory
    shutil.move("ws_log", woods_saxon_directory + "ws_log")
    shutil.move("ws_pot", woods_saxon_directory + "ws_pot")
    shutil.move("ws_rho", woods_saxon_directory + "ws_rho")
                
    # Get all s.p. states in order of lowest energy first
    sp_states, _ = order_sp_states(woods_saxon_directory, Z, N)

    # Move orbital files to Woods-Saxon directory
    for sp_state in sp_states:
            
        # Wave functions are independent of m_j, so fix m_j=j
        if sp_state.m_j == sp_state.j:
                
            file_name = get_orbital_file_name(sp_state)
            shutil.move(file_name, woods_saxon_directory + file_name)
    
    
    
if __name__ == '__main__':
    
    # nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'Be9', 4, 5
    # nucleus_name, Z, N = 'C12', 6, 6
    # nucleus_name, Z, N = 'O12', 8, 8
    # nucleus_name, Z, N = 'Ca40', 20, 20
    # nucleus_name, Z, N = 'Ca48', 20, 28
    # nucleus_name, Z, N = 'Fe56', 26, 30
    # nucleus_name, Z, N = 'Sn118', 50, 68
    nucleus_name, Z, N = 'Pb208', 82, 126
    
    # Generate orbital files
    main(nucleus_name, Z, N)