#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 09:06:21 2025

@author: anthonytropiano
"""

# Python imports
import numpy as np
from scripts.potentials import Potential
from scripts.srg import SRG


# kvnn, kmax, kmid, ntot = 6, 25.0, 4.0, 120
# kvnn, kmax, kmid, ntot = 7, 15.0, 3.0, 120
kvnn, kmax, kmid, ntot = 999, 25.0, 4.0, 120
# channels = ['1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3']
# channel = '3S1'
channels = ['1S0', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3']
lamb = 1.5
lambda_array = np.array([lamb])
# lambda_array = np.array([4.0, 3.0, 2.0, 1.5, 1.2])

for channel in channels:
    
    potential = Potential(kvnn, channel, kmax, kmid, ntot)
    d = SRG(potential, 'T')(lambda_array, save=True)

    # # Zero out off-diagonl blocks
    # V_matrix = potential.load_potential()
    # if potential.coupled_channel_bool:
    #     V_matrix[:ntot, ntot:] = np.zeros((ntot, ntot))
    #     V_matrix[ntot:, :ntot] = np.zeros((ntot, ntot))

    # potential_new = Potential(999, channel, kmax, kmid, ntot)

    # potential_new.save_potential(V_matrix, 'initial', 'T', lamb)

# potential = Potential(kvnn, channel, kmax, kmid, ntot)
# d = SRG(potential, 'T')(lambda_array, save=True)

# kvnn, kmax, kmid, ntot = 6, 25.0, 4.0, 120
# channel = '3S1'
# lamb = 1.5

# potential = Potential(kvnn, channel, kmax, kmid, ntot)



# # Zero out off-diagonl blocks
# V_matrix[:ntot, ntot:] = np.zeros((ntot, ntot))
# V_matrix[ntot:, :ntot] = np.zeros((ntot, ntot))

# potential_new = Potential(999, channel, kmax, kmid, ntot)

# potential_new.save_potential(V_matrix, 'initial', 'T', lamb)