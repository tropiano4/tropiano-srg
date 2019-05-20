# Scratch Work

# Test codes in progress in this script.
# Last thing tested: Scattering state momentum distribution for 3S1 potential


from os import chdir, getcwd
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
# Python scripts made by A.T.
from Figures import figures_functions as ff
from Potentials.vsrg_macos import load_save_potentials as lp


def wave_func(H_matrix, index, U=np.empty(0)):
    '''Diagonalizes the Hamiltonian and returns the ith wave function as u 
    and w corresponding to the 3S1 and 3D1 channels where index is the ith 
    sorted eigenvalue. The wave function is unitless, that is, the momenta 
    and weights are factored in such that \sum_i { u(k_i)^2 + w(k_i)^2 } = 1. 
    For an evolved wave function, enter in a unitary transformation U. Also
    returns the eigenvalue in MeV.'''
        
    # Arguments
        
    # U (2-D NumPy array): Option to apply an SRG or Magnus unitary
    # transformation to the wave function

    # Diagonalize Hamiltonian
    eigenvalues,eigenvectors = la.eig(H_matrix)
    
    eigenvalue = np.sort(eigenvalues)[index]
    eigenvector_index = list(eigenvalues).index(eigenvalue)
    psi = eigenvectors[:,eigenvector_index]

    # Evolve wave function by applying unitary transformation U
    if U.any():
        psi = U @ psi
        
    u = psi[:120] # 3S1 part 
    w = psi[120:] # 3D1 part
        
    # Check normalization in momentum-space
    #normalization = np.sum((u**2+w**2))
    #print('Normalization = %.4f (k-space)'%normalization)
            
    return u, w, eigenvalue
    
    
def main(index, kvnn, kmax, kmid, method, generator, lamb, lambda_bd=0.00):
    
    cwd = getcwd()
    
    # Specify the potential
    channel = '3S1'
    ntot = 120
    
    # Load initial and evolved Hamiltonian
    H0_matrix = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot)
    Hs_matrix = lp.load_hamiltonian(kvnn, channel, kmax, kmid, ntot, method, 
                                    generator, lamb, lambda_bd)
    
    # Load momentum and weights
    gp, gw = lp.load_momentum(kvnn, channel, kmax, kmid, ntot)
    
    # Compute initial and evolved wave functions in S and D wave components: u and w
    u_0, w_0, epsilon_0 = wave_func(H0_matrix, index)
    u_s, w_s, epsilon_s = wave_func(Hs_matrix, index) 
    # Initial and evolved momentum distribution (divide by momenta and weights for mesh-independent result)
    phi_squared_0 = ( u_0**2 + w_0**2 ) / ( gp**2 * gw )
    phi_squared_s = ( u_s**2 + w_s**2 ) / ( gp**2 * gw )

    # Limits of x and y axes
    xlim = [0.0,4.0]
    ylim = [1e-6,1e2]
    
    # Generator label
    if generator == 'Block-diag':
        generator_label = r'$G=H_{BD}$'+'\n'+r'$\Lambda=%.2f \/ fm^{-1}$'%lambda_bd
    elif generator == 'T':
        generator_label = r'$G=T_{rel}$'
    elif generator == 'Wegner':
        generator_label = r'$G=H_{D}$'
    
    # Fontsize for labels
    legend_label_size = 16
    x_label_size = 18
    y_label_size = 20
    generator_label_size = 18
    eigenvalue_label_size = 20
    
    # Location of legend and generator labels
    legend_label_location = 2
    generator_label_location = 6
    
    # Label for eigenvalue
    eigenvalue_label = r'$\epsilon = %.3f$'%epsilon_0+' MeV'
    
    # Plot psi^2 [fm^3] as a function of k [fm^-1]
    plt.close('all')
    
    f, ax = plt.subplots()
    
    ax.semilogy(gp, phi_squared_s, 'r-', linewidth=1.5, label=r'$\lambda=%.1f \/ \/ fm^{-1}$'%lamb)
    ax.semilogy(gp, phi_squared_0, 'k:', label=r'$\lambda=\infty \/ \/ fm^{-1}$')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc=legend_label_location, frameon=False, fontsize=legend_label_size)
    ax.set_title(eigenvalue_label, fontsize=eigenvalue_label_size)
    ax.set_xlabel(r'$k \/ [fm^{-1}]$', fontsize=x_label_size)
    ax.set_ylabel(r'$\phi^2 \/ \/ [fm^3]$', fontsize=y_label_size)
    anchored_text = AnchoredText(generator_label, prop=dict(size=generator_label_size), 
                                 loc=generator_label_location, frameon=False)
    ax.add_artist(anchored_text)
    
    # Name of jpeg file
    if generator == 'Block-diag':
        name = 'scattering_state_momentum_distribution_kvnn%s_%s_%s%.2f_lamb%.1f'%(
                kvnn, method, generator, lambda_bd, lamb)
    else:
        name = 'scattering_state_distribution_kvnn%s_%s_%s_lamb%.1f'%(kvnn, 
                method, generator, lamb)
        
    # Save figure
    folder = ff.current_date() # Gets current month and year
    chdir('Figures/%s'%folder)
    f.savefig(name+'.jpg', bbox_inches='tight')
    chdir(cwd)    
    

if __name__ == '__main__':
    
    #index = 30 # epsilon = 5.358 MeV
    #index = 60 # epsilon = 58.562 MeV
    index = 120 # epsilon = 191.947 MeV
    
    main(index, 112, 8.0, 2.0, 'magnus', 'Wegner', 1.5)
    main(index, 112, 8.0, 2.0, 'srg', 'Block-diag', 1.5, lambda_bd=2.00)
    main(index, 112, 8.0, 2.0, 'srg', 'Block-diag', 1.5, lambda_bd=3.00)