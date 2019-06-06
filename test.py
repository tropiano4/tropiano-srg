#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes.
#
# Last thing tested:
#   Making sure Magnus files are correct by over-writing the current files in
#   Potentials with old ones from the local folder (kvnn = 900, 901 for G = 
#   Wegner and T). DID NOT WORK. Will try over-writing the new momenta/weight
#   files with the old ones which go to several more decimals places. DID NOT
#   WORK. Try re-evolving these magnus files with less precise momentum/weights
#   file.
#
#------------------------------------------------------------------------------


import os
import numpy as np


cwd = os.getcwd()

os.chdir('../../')

gp = np.loadtxt('gp.dat')
gw = np.loadtxt('gw.dat')

os.chdir(cwd)
os.chdir('Potentials/vsrg_macos')

#kvnn = 900
kvnn = 901
#kvnn = 902

folder_name = 'vsrg_kvnn_%d_lam12.0_kmax30_kmid4_ntot120'%kvnn
os.chdir(folder_name)

file_name = 'vsrg_3S1_kvnn_%d_lam12.0_reg_0_3_0_mesh.out'%kvnn
f = open(file_name, 'w')
for i in range(len(gp)):
    
    k = gp[i]
    w = gw[i]
    line = '{:^25.18e}\t{:^25.18e}'.format(k, w)
    f.write(line+'\n')
    
f.close()

os.chdir(cwd)