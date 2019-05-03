# Scratch Work


import matplotlib.pyplot as plt
from Potentials.vsrg_macos import load_save_potentials as lp
from deuteron import Deuteron
from SRG_codes.srg_unitary_transformation import SRG_unitary_transformation


# Test new deuteron script with EM N3LO potential
kvnn = 10 
channel = '3S1'
kmax = 30.0
kmid = 4.0
ntot = 120

method = 'srg'
generator = 'Wegner'
lamb = 1.2

k_array, k_weights = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
H_matrix = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)

D = Deuteron(H_matrix, k_array, k_weights)

# Unevolved
u, w = D.wave_func()

# Test momentum distribution
psi_squared_units = D.momentum_distribution(u, w)

plt.semilogy(k_array, psi_squared_units)
plt.xlim([0.0,3.0])
plt.ylim([1e-5,1e3])
plt.xlabel(r'$k [fm^{-1}]$')
plt.ylabel(r'$|\phi(k)|^2$')
plt.show()