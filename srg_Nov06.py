# Created 11/06/18 by A.T. (tropiano.4@osu.edu)

# SRG code: Evolves Hamiltonian to diagonal or block-diagonal, decoupled form 
# with parameter s.
# Main function: SRG evolves specified potential and saves evolved, vectorized 
# Hamiltonian.

import numpy as np
from os import chdir,getcwd
from numba import jit
from scipy.integrate import odeint
import time


# --------------------------------------------------------------------------- #
# Initialize global variables here by setting placeholders
# These variables are used by multiple functions and must be declared first

potential_type = 'potential_type' # Specifies the potential
T0_matrix = np.zeros((240,240)) # Initialize relative kinetic energy
srg_generator = 0 # Specifies the SRG generator
decouple = 0.0 # Block-diagonal lambda for decoupling
# Placeholders for block-diagonal projection matrices
P = np.zeros((240,240))
Q = np.zeros((240,240))
# --------------------------------------------------------------------------- #


def load_hamiltonian(lambda_eft):
    '''Load initial Hamiltonian for given potential type ('EGM','RKE','Wendt',
    etc.) and EFT cutoff lambda_eft [fm^-1] ('4', '9', '20', etc.). Returns
    H0_matrix in units fm^-2.'''
    
    cwd = getcwd()
    # potential_type is a global variable
    chdir('Potentials_in_MeV/%sFiles'%(potential_type))
    # Matrices are in units MeV: convert to fm^-2
    hb2m = 41.47 # Units are MeV fm^2
        
    # Relative kinetic energy
    global T0_matrix 
    T0_matrix = np.loadtxt('T.dat')/hb2m 
    # Argonne v18 potential
    if lambda_eft == '0': 
        V0_matrix = np.loadtxt('V.dat')/hb2m
    # Potential corresponding to EFT cutoff lambda_eft
    else:
        V0_matrix = np.loadtxt('V_%s.000000.dat'%(lambda_eft))/hb2m
    
    chdir(cwd)
    
    # Initial Hamiltonian
    H0_matrix = T0_matrix+V0_matrix
    
    return H0_matrix   

@jit(nopython=True)     
def commutator(A,B):
    '''Returns commutator of A and B, [A,B].'''
     
    return np.dot(A,B)-np.dot(B,A)   

@jit(nopython=True) 
def matrix2vector(A):
    '''Takes the top right of the matrix A (including the diagonal) and 
    reshapes it into a vector B of dimension N*(N+1)/2.'''
    
    # Dimension of matrix
    N = len(A)
    # Dimension of vectorized matrix
    n = int(N*(N+1)/2)
    # Initialize vectorized matrix
    B = np.zeros(n)
    
    first = 0
    last = N
    
    for i in range(N):
        
        B[first:last] = A[i][i:]
        first = last
        last += N-i-1

    return B

@jit(nopython=True)    
def vector2matrix(B):
    '''Takes the vector of a top right matrix and returns the full matrix. (USE
    ONLY FOR HERMITIAN MATRICES.)'''
        
    # Dimension of the vectorized matrix
    n = len(B)
    # Dimension of matrix (given by solving N*(N+1)/2 = n)
    N = int( (-1 + np.sqrt(1+8*n) ) / 2 )
    # Initialize matrix
    A = np.zeros((N,N))
    
    # Build upper half of A with diagonal

    first = 0
    last = N

    for i in range(N):

        A[i,i:] = B[first:last]
        first = last
        last += N-i-1

    # Reflect upper half to lower half to build full matrix
    # [np.transpose(A)-np.diag(np.diag(A))] is the lower half of A excluding 
    # the diagonal
    return A+(np.transpose(A)-np.diag(np.diag(A)))

def derivs(Hs_vector,s):
    '''Returns RHS of SRG flow equation where Hs_vector is the solution vector.'''
        
    # Matrix of the solution vector
    Hs_matrix = vector2matrix(Hs_vector)

    # Note: srg_generator is a global variable
    # SRG generator, eta = [G,H] 
    if srg_generator == 1:
        # G = Diagonal matrix of H(s)
        G = np.diag(np.diag(Hs_matrix))
            
    elif srg_generator == 2:
        # G = T_rel
        G = T0_matrix # Global variable
            
    elif srg_generator == 3:
        # G = P*H(s)*P+Q*H(s)*Q
        G = np.dot(P,np.dot(Hs_matrix,P))+np.dot(Q,np.dot(Hs_matrix,Q))
            
    eta = commutator(G,Hs_matrix)
            
    dH_matrix = commutator(eta,Hs_matrix)
        
    # Returns vector form of RHS of flow equation
    dH_vector = matrix2vector(dH_matrix)
        
    return dH_vector
    
def evolve_hamiltonian(H0_matrix,lambda_array):
    '''Returns evolved vector Hs_vector at several values of lambda for given
    lambda array.'''

    # Set-up ODE
    # Reshape initial hamiltonian to a vector
    H0_vector = matrix2vector(H0_matrix)
        
    # Evaluate H(s) at the following values of lambda (or s)
    s_array = np.zeros(len(lambda_array)+1)
    s_array[1:] = 1.0/lambda_array**4.0 # This array includes s = 0

    # Solve the flow equations
    # Sol returns a vectorized H(s) at several values of s
    sol = odeint(derivs,H0_vector,s_array,atol=1e-06,rtol=1e-06,mxstep=5000000)

    # Return a dictionary of H(s) at the values of s (i.e., d[1.2] returns 
    # H(lambda=1.2)) which is a vector
    d = {}
    i = 1
    for l in lambda_array:
            
        d[l] = sol[i] # H( s = s_array[i] )
        i += 1
        
    return d
    
def write_data(lambda_eft,gen,lambda_array,d):
    '''Write H(s) data to files where H(s) is a vector of the top right part of
    H(s) matrix.'''
    
    cwd = getcwd()
    chdir('Evolved_Hamiltonians')
    
    p = potential_type # Name of potential (global)
    l = lambda_eft # EFT cutoff lambda
    g = gen # SRG generator
    
    # Loop over values of lambda
    for lamb in lambda_array:
        
        Hs_vector = d[lamb] # Top right vector
        
        # Name of file
        if g == 'Block-diag':
            # File name (e.g., 'srg_Wendt_Lamb4_Block-diag2.00_lamb1.2.dat')
            if p == 'AV18' or p == 'EM':
                name = 'srg_%s_%s%.2f_lamb%.1f'%(p,g,decouple,lamb) 
            else:
                name = 'srg_%s_Lamb%s_%s%.2f_lamb%.1f'%(p,l,g,decouple,lamb) 
        else:
            # File name (e.g., 'srg_RKE_Lamb9_Wegner_lamb1.2.dat')
            if p == 'AV18' or p == 'EM':
                name = 'srg_%s_%s_lamb%.1f'%(p,g,lamb) 
            else:
                name = 'srg_%s_Lamb%s_%s_lamb%.1f'%(p,l,g,lamb)
           
        # Write to file
        f = open(name+'.dat','w')
        
        for val in Hs_vector:
            f.write(str(val)+'\n')
            
        f.close()

    chdir(cwd)  
    
def main(p_type,lambda_eft,gen,dcpl=2.00,save=True):
    '''
    Main program: SRG evolve potential of given type (p_type which is a string) 
    at EFT cutoff lambda_eft (string). Use SRG generator specified by gen 
    (string). If srg_generator = 'Block-diag', need to specify where to split 
    sub-blocks with float dcpl. The main function saves the evolved vectorized 
    Hamiltonian in units fm^-2 for lambda = 10, 2.8, 2.0, and 1.2 fm^-1 in the 
    folder Evolved_Hamiltonians if save = True.
    '''
    
    # Save arguments globally
    global potential_type
    potential_type = p_type
    global srg_generator
    if gen == 'Wegner':
        srg_generator = 1
    elif gen == 'T':
        srg_generator = 2
    elif gen == 'Block-diag':
        srg_generator = 3
        global decouple # For block-diagonal flow only
        decouple = dcpl

    l = lambda_eft # EFT cutoff lambda
    
    # Initial Hamiltonian H(0)
    H0_matrix = load_hamiltonian(lambda_eft)
    # Dimension of matrix
    N = len(H0_matrix)
    
    # For block-diagonal generator, initialize projection matrices
    if gen == 'Block-diag':
        
        cwd = getcwd()
        chdir('Potentials_in_MeV/%sFiles'%(p_type))
    
        # Length of momentum array
        n = int(N/2) 
        # Load momentum array, gp
        gp = np.loadtxt('gp.dat')
        
        chdir(cwd)

        # Array of True and False
        bool_array = gp < decouple
        # Build projection operators
        # Identity matrix for k < decouple
        p = np.diag(np.ones(n)*bool_array)
        # Opposite of p: Identity matrix for k > decouple
        q = np.identity(n)-p
        o = np.zeros((n,n))
        # Construct full matrices for coupled channel vnn
        global P
        P = np.vstack((np.hstack((p,o)),np.hstack((o,p))))
        global Q
        Q = np.vstack((np.hstack((q,o)),np.hstack((o,q))))
    
    # Lambda values in units fm^-1
    lambda_array = np.array([10.0,2.8,2.0,1.2])
    
    # Evolve Hamiltonian at values in lambda_array
    t0 = time.time() # Start time
    d = evolve_hamiltonian(H0_matrix,lambda_array)
    t1 = time.time() # End time
    
    mins = round((t1-t0)/60.0,2) # Minutes elapsed evolving H(s)
    print('_'*80)
    if gen == 'Block-diag':
        print('SRG: %s potential, Lambda = %s fm^-1, generator: %s%.2f'% \
        (p_type,l,gen,decouple))
    else:
        print('SRG: %s potential, Lambda = %s fm^-1, generator: %s'% \
        (p_type,l,gen))
    print('H(s) done evolving to lambda = 1.2 fm^-1 after %f minutes'%mins)
    
    if save:
        # Write and save data 
        write_data(lambda_eft,gen,lambda_array,d)
    else:
        # Return dictionary
        return d
        
        
if __name__ == '__main__':
    
    #p = 'AV18'
    p = 'Wendt'
    #p = 'EM'
    #l_list = ['4','9','20']
    #l = '0'
    l = '4'
    #l = '9'
    #g_list  = ['Wegner','T']
    #g = 'Wegner'
    g = 'T'

    #for l in l_list:
        #for g in g_list:
    
            #main(p,l,g)
    main(p,l,g)
