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


kvnn, kmax, kmid, ntot = 6, 25.0, 4.0, 120
# kvnn, kmax, kmid, ntot = 7, 15.0, 3.0, 120
channels = ['1S0', '3S1', '3P0', '1P1', '3P1', '3P2', '1D2', '3D2', '3D3']
lambda_array = np.array([6.0])
# channels = ['3F3', '1F3', '3F4']
# lambda_array = np.array([4.0, 3.0, 2.0, 1.5, 1.2])

for channel in channels:
    
    potential = Potential(kvnn, channel, kmax, kmid, ntot)
    d = SRG(potential, 'T')(lambda_array, save=True)