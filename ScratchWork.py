# Scratch Work

# Test codes in progress in this script.


import numpy as np
import numpy.linalg as la
from Potentials.vsrg_macos import load_save_potentials as lp
from Magnus_codes.magnus_wegner_may10 import magnus_wegner as mw


hbar_sq_over_M = 41.47

# Hamiltonian specifications
kvnn = 119
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120
# Evolve to lambda 
lamb = 3.0

H0_matrix = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
eig_init = np.sort( la.eig(H0_matrix)[0] )

Hs_matrix = mw.evolve_hamiltonian(H0_matrix, lamb)
eig_final = np.sort( la.eig(Hs_matrix*hbar_sq_over_M)[0] )

print(eig_init[:5])
print(eig_final[:5])