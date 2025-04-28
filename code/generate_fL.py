#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: generate_longitudinal_sf.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: April 25, 2025

Script for running the deuteron electrodisintegration longitudinal structure
function code.

Last update: April 25, 2025

"""

# Python imports
import numpy as np
import time

# JAX imports
from jax import config
import jax.numpy as jnp

# Imports from scripts
from scripts.deuteron_electrodisintegration import DeuteronElectrodisintegration


# Enable double precision
config.update("jax_enable_x64", True)


def write_data_wrt_thetap(kvnn, kmax, kmid, ntot, lamb, Ep, thetap_array, q,
                          L_max=2):
    """Writes data for f_L at the given kinematics w.r.t. \theta'."""

    # Initialize classes
    de_ia_unevolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=jnp.inf, L_max=L_max, option=1
    )
    de_ia_evolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=lamb, L_max=L_max, option=2
    )
    de_fsi_unevolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=jnp.inf, L_max=L_max, option=3
    )
    de_fsi_evolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=lamb, L_max=L_max, option=4
    )
    
    # Compute f_L for each case
    f_L_ia_unevolved = de_ia_unevolved.fL_vmap_theta(Ep, thetap_array, q)
    print("Done with IA + unevolved.")
    f_L_ia_evolved = de_ia_evolved.fL_vmap_theta(Ep, thetap_array, q)
    print("Done with IA + evolved.")
    f_L_fsi_unevolved = de_fsi_unevolved.fL_vmap_theta(Ep, thetap_array, q)
    print("Done with FSI + unevolved.")
    f_L_fsi_evolved = de_fsi_evolved.fL_vmap_theta(Ep, thetap_array, q)
    print("Done with FSI + evolved.")
    
    # Save data
    directory = '../data/longitudinal_sf/'
    filename = f"fL_lamb{lamb:.1f}_Ep{Ep:d}_q{q:.2f}_Lmax{L_max:d}.txt"
    
    data = np.vstack((
        thetap_array, f_L_ia_unevolved, f_L_ia_evolved, f_L_fsi_unevolved,
        f_L_fsi_evolved
    )).T
    hdr = ("theta' [deg], IA + Unevolved [fm], IA + Evolved [fm],"
           " FSI + Unevolved [fm], FSI + Evolved [fm]\n")
    np.savetxt(directory + filename, data, header=hdr)
    
    
def load_data_wrt_thetap(lamb, Ep, q, L_max=2):
    """Load data for f_L at the given kinematics w.r.t. \theta'."""
    
    # Load data
    directory = '../data/longitudinal_sf/'
    filename = f"fL_lamb{lamb:.1f}_Ep{Ep:d}_q{q:.2f}_Lmax{L_max:d}.txt"
    
    data = np.loadtxt(directory + filename)
    
    thetap_array = data[:, 0]
    f_L_ia_unevolved = data[:, 1]
    f_L_ia_evolved = data[:, 2]
    f_L_fsi_unevolved = data[:, 3]
    f_L_fsi_evolved = data[:, 4]
    
    return (thetap_array, f_L_ia_unevolved, f_L_ia_evolved, f_L_fsi_unevolved,
            f_L_fsi_evolved)
    
    
def write_data_wrt_Ep(kvnn, kmax, kmid, ntot, lamb, Ep_array, thetap, L_max=2):
    """Writes data for f_L at the given kinematics w.r.t. E'."""

    # Initialize classes
    de_ia_unevolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=jnp.inf, L_max=L_max, option=1
    )
    de_ia_evolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=lamb, L_max=L_max, option=2
    )
    de_fsi_unevolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=jnp.inf, L_max=L_max, option=3
    )
    de_fsi_evolved = DeuteronElectrodisintegration(
        kvnn, kmax, kmid, ntot, lamb=lamb, L_max=L_max, option=4
    )
    
    # Compute f_L for each case
    t0 = time.time()
    f_L_ia_unevolved = de_ia_unevolved.fL_quasifree_ridge(Ep_array, thetap)
    t1 = time.time()
    mins = (t1 - t0) / 60
    print(f"Done with IA + unevolved after {mins:.2f}.")
    
    t0 = time.time()
    f_L_ia_evolved = de_ia_evolved.fL_quasifree_ridge(Ep_array, thetap)
    t1 = time.time()
    mins = (t1 - t0) / 60
    print(f"Done with IA + evolved after {mins:.2f}.")
    
    t0 = time.time()
    f_L_fsi_unevolved = de_fsi_unevolved.fL_quasifree_ridge(Ep_array, thetap)
    t1 = time.time()
    mins = (t1 - t0) / 60
    print(f"Done with FSI + unevolved after {mins:.2f}.")
    
    t0 = time.time()
    f_L_fsi_evolved = de_fsi_evolved.fL_quasifree_ridge(Ep_array, thetap)
    t1 = time.time()
    mins = (t1 - t0) / 60
    print(f"Done with FSI + evolved after {mins:.2f}.")
    
    # Save data
    directory = '../data/longitudinal_sf/'
    filename = f"fL_lamb{lamb:.1f}_qfr_thetap{thetap:.1f}_Lmax{L_max:d}.txt"
    
    data = np.vstack((
        Ep_array, f_L_ia_unevolved, f_L_ia_evolved, f_L_fsi_unevolved,
        f_L_fsi_evolved
    )).T
    hdr = ("E' [MeV], IA + Unevolved [fm], IA + Evolved [fm],"
           " FSI + Unevolved [fm], FSI + Evolved [fm]\n")
    np.savetxt(directory + filename, data, header=hdr)
    
def load_data_wrt_Ep(lamb, thetap, L_max=2):
    """Load data for f_L at the given kinematics w.r.t. E'."""
    
    # Load data
    directory = '../data/longitudinal_sf/'
    filename = f"fL_lamb{lamb:.1f}_qfr_thetap{thetap:.1f}_Lmax{L_max:d}.txt"
    
    data = np.loadtxt(directory + filename)
    
    Ep_array = data[:, 0]
    f_L_ia_unevolved = data[:, 1]
    f_L_ia_evolved = data[:, 2]
    f_L_fsi_unevolved = data[:, 3]
    f_L_fsi_evolved = data[:, 4]
    
    return (Ep_array, f_L_ia_unevolved, f_L_ia_evolved, f_L_fsi_unevolved,
            f_L_fsi_evolved)


if __name__ == '__main__':
    
    # Potential
    kvnn, kmax, kmid, ntot = 6, 25.0, 4.0, 120
    
    # # SRG \lambda [fm^-1]
    # lamb = 1.5
    # # lamb = 15.0

    # # Truncation on sums involving total orbital angular momentum L
    # L_max = 2
    # # L_max = 3

    # # \theta values for fixed q and E'
    # thetap_array = jnp.linspace(0.01, 179.9, 100)
    # Ep_values = [100, 10, 30]
    # q_values = [jnp.sqrt(10), 2, 4]
    # for iEp, jq in zip(Ep_values, q_values):
    #     write_data_wrt_thetap(kvnn, kmax, kmid, ntot, lamb, iEp, thetap_array,
    #                           jq, L_max)
        
    # # Quasifree ridge where \omega = 0
    # Ep_array = np.linspace(10, 115, 100)
    # thetap = 15.0
    # write_data_wrt_Ep(kvnn, kmax, kmid, ntot, lamb, Ep_array, thetap, L_max)
    
    ### TESTING JUST B + F terms
    L_max = 2
    ded = DeuteronElectrodisintegration(kvnn, kmax, kmid, ntot, lamb=jnp.inf,
                                        L_max=L_max, option=4)
    f_L = ded.fL(10.0, 15.0, 2.0)
    print(f_L)
