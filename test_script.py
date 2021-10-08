#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a separate file with an
# extension _test.py. For example, momentum_projection_operator_testv1.py (v1
# means 'version 1'). Use the revision history below to document when and why
# these files are created.
#
# Revision history:
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py based off last tests in
#                this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py.
#   01/26/21 --- Renamed to test_script.py.
#   03/16/21 --- Testing pair momentum distributions for different nuclei
#                using LDA with simple filled Fermi sea single-particle
#                momentum distributions. Created pmd_simple_test.py.
#   03/16/21 --- Testing pair momentum distribution for deuteron starting from
#                second quantization derivation and using expansion of
#                U(k, k'). Created pmd_deuteron_test.py.
#   03/16/21 --- Testing single-nucleon momentum distributions for different
#                nuclei using LDA. Created single_particle_momentum_dist.py.
#   04/14/21 --- Creating AV18 SRG evolution figure for APS April Meeting
#                presentation.
#                potential_contours_kvnn_6_channel_1P1_Wegner_lamb_1p5.png in
#                figures/operator_evolution/old/potential_contours.
#   04/28/21 --- Testing normalization and contributions of \delta U, etc. or
#                pp/pn to single-nucleon momentun distributions. Created
#                lda_normalizations_test.py.
#   05/04/21 --- Testing higher partial waves of SRG transformations: 3P2-3F2
#                and 3D3-3G3 have numerical artifacts. Created
#                high_partial_waves_srg_test.py.
#   06/10/21 --- Verifying \theta functions averaging in snmd.py and dmd.py by
#                comparing numerical functions to analytic evaluation of 
#                \int d3K \int d3k \theta(kF-|K/2+k|) \theta(kF-|K/2-k|).
#                Created theta_functions_test.py.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Testing Fourier transformation of s.p. wave functions from SLy4 and HFBRAD
#   code.


import matplotlib.pyplot as plt
import numpy as np
from scipy.special import spherical_jn
# Scripts made by A.T.
from misc.integration import gaussian_quadrature_mesh
from potentials.vsrg_macos import vnn


def hankel_transformation(L, k_array, r_array, dr):
    """
    <k|r> transformation matrix for given orbital angular momentum L. If
    len(r_array) = m and len(k_array) = n, then this function returns an 
    n x m matrix.
    
    Parameters
    ----------
    L : int
        Orbital angular momentum.
    k_array : 1-D ndarray
        Momentum array [fm^-1].
    r_array : 1-D ndarray
        Coordinates array [fm].
    dr : float
        Coordinates step-size (weight) [fm].
        
    Returns
    -------
    M : 2-D ndarray
        Hankel transformation matrix [fm^3\2].

    """

    # r_array column vectors and k_array row vectors where both grids are
    # n x m matrices
    r_cols, k_rows = np.meshgrid(r_array, k_array)
        
    # M = dr * r_cols**2 * spherical_jn(L, k_rows * r_cols)
    M = dr * r_cols * spherical_jn(L, k_rows * r_cols)

    return M

# Set nucleus
nucleus = 'O16'

# Load phi_\alpha(r)
data_directory = 'densities/HFBRAD_SLY4/%s/wfs/' % nucleus

phi_data = np.loadtxt(data_directory+'qp0006_8_8.gfx')
r_array = phi_data[:, 0]
dr = r_array[1] - r_array[0]
phi_r_array = phi_data[:, 2]
# for ir, iphi in zip(r_array, phi_r_array):
#     print(ir, iphi)
    
# Compute normalization of coordinate-space WF
# norm_r = np.sum(dr*r_array**2*phi_r_array)
norm_r = np.sum(dr*phi_r_array**2)
print(norm_r)
# Normalization is \int dr u(r)^2 = 2j+1?

# Do I use momentum values possibly negative for p_missing? or same as
# \delta U?
# k_array, k_weights = vnn.load_momentum(6, '1S0', 15.0, 3.0, 120)
k_array, k_weights = gaussian_quadrature_mesh(10.0, 200)

# Compute FT of WF
ft_matrix = hankel_transformation(0, k_array, r_array, dr)
phi_k_array = ft_matrix @ phi_r_array

# Compute normalization of momentum-space WF
norm_k = 2/np.pi * np.sum(k_weights*k_array**2*phi_k_array**2)
print(norm_k)

# Plot wave functions
plt.clf()
plt.plot(r_array, phi_r_array, label=r'$1s_{1/2}$')
plt.ylabel(r'$\phi(r)$' + ' [fm' + r'$^{-1/2}$' + ']')
plt.xlabel('r [fm]')
plt.legend(loc='upper right')
plt.show()

plt.clf()
# plt.plot(k_array, phi_k_array, label=r'$1s_{1/2}$')
plt.semilogy(k_array, phi_k_array**2, label=r'$1s_{1/2}$')
plt.ylabel(r'$|\phi(k)|^2$' + ' [fm' + r'$^{3}$' + ']')
plt.xlabel('k [fm' + r'$^{-1}$' + ']')
plt.legend(loc='upper right')
plt.show()