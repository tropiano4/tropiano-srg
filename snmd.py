#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: snmd.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     March 25, 2021
# 
# Calculates SRG-evolved single-nucleon momentum distributions for nuclei
# assuming the evolved wave function is given by a free Fermi gas. (Note,
# 'snmd' stands for single-nucleon momentum distribution.) Combine these
# functions with LDA.py for nuclear-averaged momentum distributions. This
# script is based off single_particle_momentum_dist.py in Old_codes.
#
# Revision history:
#   XX/XX/XX --- ...
#
#------------------------------------------------------------------------------


import numpy as np
from scipy.interpolate import RectBivariateSpline
# Scripts made by A.T.
from Potentials.vsrg_macos import vnn
from SRG.srg_unitary_transformation import SRG_unitary_transformation


class single_nucleon_momentum_distributions(object):
    
    
    def __init__(self):
        """
        """
        
        return None
    
    
    # Interpolation of \delta U's
    # Evaluation of matrix elements
    # Grids for \theta functions and angle-averaging
    # Return n(q, kF_1, kF_2) total or 1, \delta U, \delta U\delta U^{\dagger}
    # or pp, nn, pn, and np contributions
    
    def n_lambda(self, q, kF_1, kF_2):
        """
        """
        
        return None
    
    
    
    
if __name__ == '__main__':
    
    x = None