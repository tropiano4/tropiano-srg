#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: generate_woods_saxon.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: July 26, 2023

Run the Woods-Saxon Fortran code for some nucleus. See
https://nucracker.volya.net or arXiv:0709.3525 [nucl-th]
to set parameters of nuclei.

Last update: February 16, 2024

"""

# Python imports
import numpy as np
import shutil

# Imports from scripts
from .momentum_distributions import SingleParticleState
from .woodsaxon import ws


def get_orbitals():
    """Return lists of orbitals."""

    # Track 2*j not j
    norb, lorb, jorb = [], [], []
    # These must match the nmax and lmax in woodsaxon.f90
    nmax, lmax = 3, 6
    for n in range(1, nmax+1):
        for l in range(0, lmax+1):
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


def set_seminole_parameters(nucleus_name, A):
    """Sets the Woods-Saxon parameters according to the 'Seminole'
    parametrization. The first index of the parameters array 'prm' corresponds
    to proton (0) or neutron (1). The second index corresponds to the following:
        
        0 : Central strength V_0 [MeV]
        1 : Central radius R_0 [fm]
        2 : Central surface diffuseness a [fm]
        3 : Wine-Bottle kwb (set to 0 for no Wine-Bottle potential)
        4 : Wine-Bottle awb
        5 : Spin-orbit strength \lambda [MeV]
        6 : Spin-orbit radius R_SO [fm]
        7 : Spin-orbit surface diffuseness (same as a) [fm]
        8 : Coulomb radius (same as R_0) [fm]
    
    """

    prm = np.zeros(shape=(2, 9), order='F')
    
    # These parameters are the same for every nucleus
    a = 0.644174
    kwb = 0.0
    awb = 1.0
    prm[:, 2] = a
    prm[:, 3] = kwb
    prm[:, 4] = awb
    prm[:, 7] = a
    
    # Nucleus-dependent parameters
    if nucleus_name == 'He4':
        prm[:, 0] = 76.8412
        prm[:, 1] = 2.00013 / (A ** (1/3))
        prm[0, 5] = 51.8575
        prm[1, 5] = 51.7507
        prm[:, 6] = 1.90409 / (A ** (1/3))
        prm[:, 8] = 2.00013 / (A ** (1/3))
    elif nucleus_name == 'C12':
        prm[:, 0] = 60.1478
        prm[:, 1] = 2.88468 / (A ** (1/3))
        prm[0, 5] = 34.6305
        prm[1, 5] = 34.5432
        prm[:, 6] = 2.74617 / (A ** (1/3))
        prm[:, 8] = 2.88468 / (A ** (1/3))
    elif nucleus_name == 'O16':
        prm[:, 0] = 58.0611
        prm[:, 1] = 3.175 / (A ** (1/3))
        prm[0, 5] = 33.0985
        prm[1, 5] = 33.0131
        prm[:, 6] = 3.02255 / (A ** (1/3))
        prm[:, 8] = 3.175 / (A ** (1/3))
    elif nucleus_name == 'Ca40':
        prm[:, 0] = 54.3051
        prm[:, 1] = 4.30914 / (A ** (1/3))
        prm[0, 5] = 30.5847
        prm[1, 5] = 30.5027
        prm[:, 6] = 4.10223 / (A ** (1/3))
        prm[:, 8] = 4.30914 / (A ** (1/3))
    else:
        raise RuntimeError("Don't have this nucleus yet.")
        
    return prm
        
        
def set_universal_parameters():
    """Sets the Woods-Saxon parameters according to the 'Universal'
    parametrization."""
    
    prm = np.zeros(shape=(2, 9), order='F')
    
    # Central potential strength V_0 [MeV]
    prm[:, 0] = 49.6
    # Central radius R_0 [fm]
    prm[0, 1] = 1.275  # Proton
    prm[1, 1] = 1.347  # Neutron
    # Central surface diffuseness a [fm]
    prm[:, 2] = 0.7
    # Wine-Bottle kwb (set to 0 for no Wine-Bottle potential)
    prm[:, 3] = 0.0
    # Wine-Bottle awb
    prm[:, 4] = 1.0
    # Spin-orbit potential strength \lambda [MeV]
    prm[0, 5] = 36  # Proton
    prm[1, 5] = 35  # Neutron
    # Spin-orbit radius R_SO [fm]
    prm[0, 6] = 1.32  # Proton
    prm[1, 6] = 1.31  # Neutron
    # Spin-orbit surface diffuseness (same as a) [fm]
    prm[:, 7] = 0.7
    # Coulomb radius (same as R_0) [fm]
    prm[0, 8] = 1.275  # Proton
    prm[1, 8] = 1.347  # Neutron
    
    return prm


def set_match_parameters(nucleus_name, A):
    """Use this to function to match to VMC and CVMC momentum distributions by
    tuning the Woods-Saxon central potential strength and central radius.
    """
    
    prm = np.zeros(shape=(2, 9), order='F')

    # Matching to VMC with AV18 + UX for A < 16
    if A < 16:
        
        # Central potential strength V_0 [MeV] (from Seminole)
        if nucleus_name == 'He4':
            prm[:, 0] = 76.8412
        elif nucleus_name == 'C12':
            prm[:, 0] = 60.1478
            
        # Central radius R_0 [fm] (from Universal)
        prm[0, 1] = 1.275  # Proton
        prm[1, 1] = 1.347  # Neutron
        
    # Matching to CVMC with AV18 + UIX for O16 and up
    else:
        
        # Central potential strength V_0 [MeV]
        prm[:, 0] = 52.0
        
        # Central radius R_0 [fm]
        prm[0, 1] = 1.55  # Proton
        # Neutron radius is scaled up by same amount as Universal
        prm[1, 1] = 1.55 * 1.0564705882352943  # Neutron
    
    # These parameters are the same as the Universal parametrization
    # Central surface diffuseness a [fm]
    prm[:, 2] = 0.7
    # Wine-Bottle kwb (set to 0 for no Wine-Bottle potential)
    prm[:, 3] = 0.0
    # Wine-Bottle awb
    prm[:, 4] = 1.0
    # Spin-orbit potential strength \lambda [MeV]
    prm[0, 5] = 36  # Proton
    prm[1, 5] = 35  # Neutron
    # Spin-orbit radius R_SO [fm]
    prm[0, 6] = 1.32  # Proton
    prm[1, 6] = 1.31  # Neutron
    # Spin-orbit surface diffuseness (same as a) [fm]
    prm[:, 7] = 0.7
    # Coulomb radius (same as R_0) [fm]
    prm[0, 8] = prm[0, 1]  # Proton
    prm[1, 8] = prm[1, 1]  # Neutron

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


def main(nucleus_name, Z, N, rmax=40, ntab=2000, parametrization='seminole'):
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
    if parametrization == 'seminole':
        prm = set_seminole_parameters(nucleus_name, A)
    elif parametrization == 'universal':
        prm = set_universal_parameters()
    elif parametrization == 'match':
        prm = set_match_parameters(nucleus_name, A)
  
    # Print summary, potentials, and densities
    prnt = True
    prntorb = True

    # Run Fortran subroutine
    ws(ntau, A, Z, rmax, orbws, norb, lorb, jorb, prm, rdiv, prnt, prntorb,
       dens)

    # Woods-Saxon directory
    directory = f"../data/woods_saxon/{parametrization}/{nucleus_name}/"
    
    # Move output files to relevant directory
    shutil.move("ws_log", directory + "ws_log")
    shutil.move("ws_pot", directory + "ws_pot")
    shutil.move("ws_rho", directory + "ws_rho")
                
    # Get all s.p. states in order of lowest energy first
    sp_states, _ = order_sp_states(directory, Z, N)

    # Move orbital files to Woods-Saxon directory
    for sp_state in sp_states:
            
        # Wave functions are independent of m_j, so fix m_j=j
        if sp_state.m_j == sp_state.j:
                
            file_name = get_orbital_file_name(sp_state)
            shutil.move(file_name, directory + file_name)
    
    
    
if __name__ == '__main__':
    
    # nucleus_name, Z, N = 'He4', 2, 2
    # nucleus_name, Z, N = 'C12', 6, 6
    # nucleus_name, Z, N = 'O16', 8, 8
    nucleus_name, Z, N = 'Ca40', 20, 20
    
    # Woods-Saxon parametrization
    prm = 'seminole'
    # prm = 'universal'
    # prm = 'match'
    
    # Generate orbital files
    main(nucleus_name, Z, N, parametrization=prm)