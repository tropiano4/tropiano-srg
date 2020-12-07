#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a seperate file with an
# extension _testv#.py where v# corresponds to the version number. For example,
# momentum_projection_operator_testv1.py. Use the revision history below to
# document when and why these files are created.
#
# Revision history:
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#   10/01/19 --- Testing SRG-evolution of an operator which is a constant at
#                all values of k and k' (this is a delta function in
#                coordinate-space). Created toy_operator_srg_evolution_v1.py in
#                Old_codes.
#   10/14/19 --- Making a couple of plots for DNP 2019 meeting.
#   03/16/20 --- Tested generalized NN operator conventions and mesh-
#                dependence. Created operators_test.py in Old_codes.
#   04/08/20 --- Tested SRG changes in r^2 operator.
#   05/22/20 --- Testing mesh-dependence in several operators. Created
#                mesh_dependence_test.py in Old_codes.
#
#------------------------------------------------------------------------------


# Description of this test:
#   Test solve_ivp vs ode with SRG code.


import numpy as np
import numpy.linalg as la
from scipy.integrate import ode, solve_ivp
import time
# Scripts made by A.T.
from Potentials.vsrg_macos import vnn


# --- Load Hamiltonian --- #

# h-bar^2 / M [MeV fm^2]
hbar_sq_over_m = 41.47

# Set potential
kvnn = 111
channel = '3S1'
kmax = 10.0
kmid = 2.0
ntot = 120
N = 2 * ntot
H_matrix_MeV = vnn.load_hamiltonian(kvnn, channel, kmax, kmid, ntot) # MeV
H_matrix = H_matrix_MeV / hbar_sq_over_m # fm^-2


# --- Functions --- #

def commutator(A, B):
    """
    Commutator of A and B, [A,B] where A and B are square matrices.
    
    Parameters
    ----------
    A : 2-D ndarray
        First input square matrix.
    B : 2-D ndarray
        Second input square matrix.
            
    Returns
    -------
    out : 2-D ndarray
        Commutator of the two input matrices.
            
    """
        
    return A @ B - B @ A

def matrix2vector(A):
    """
    Takes the upper triangle of the matrix A (including the diagonal) 
    and reshapes it into a vector B of length N*(N+1)/2.
        
    Parameters
    ----------
    A : 2-D ndarray
        Input matrix.
            
    Returns
    -------
    B : 1-D ndarray
        Output vector.
            
    """
    
    # Length of vectorized matrix
    n = int( N * (N+1) / 2 )
        
    # Initialize vectorized matrix
    B = np.zeros(n)
    
    a = 0
    b = N
    
    for i in range(N):
        
        B[a:b] = A[i][i:]
        a = b
        b += N-i-1

    return B

def vector2matrix(B):
    """
    Takes the vector of a upper triangle matrix and returns the full 
    matrix. Use only for symmetric matrices.
        
    Parameters
    ----------
    B : 1-D ndarray
        Input vector.
        
    Returns
    -------
    A : 2-D ndarray
        Output matrix.
            
    """
        
    # Initialize matrix
    A = np.zeros( (N, N) )
    
    # Build upper half of A with diagonal
    a = 0
    b = N

    for i in range(N):

        A[i, i:] = B[a:b]
        a = b
        b += N-i-1

    # Reflect upper half to lower half to build full matrix
    # [np.transpose(A)-np.diag(np.diag(A))] is the lower half of A 
    # excluding the diagonal
    return A + ( np.transpose(A) - np.diag( np.diag(A) ) )

def derivative_method_1(lamb, H_evolved):
    """
    Right-hand side of the SRG flow equation using the Wegner generator.
        
    Parameters
    ----------
    lamb : float
        Evolution parameter lambda [fm^-1].
    H_evolved : 1-D ndarray
        Evolving Hamiltonian which is a vector and function of lambda [fm^-2].
        
    Returns
    -------
    dH_vector : 1-D ndarray
        Derivative with respect to lambda of the evolving Hamiltonian which 
        is a vector [fm^-2].

    """
        
    # Matrix form of the evolving Hamiltonian
    H_evolved_matrix = vector2matrix(H_evolved)

    # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
    # evolving Hamiltonian)
    eta = commutator( np.diag( np.diag(H_evolved_matrix) ), H_evolved_matrix)
            
    # RHS of flow equation in matrix form
    dH_matrix = -4.0 / lamb**5 * commutator(eta, H_evolved_matrix)
        
    # Returns vector form of RHS of flow equation
    dH_vector = matrix2vector(dH_matrix)
        
    return dH_vector

def derivative_method_2(s, H_evolved):
    """
    Right-hand side of the SRG flow equation using the Wegner generator.
        
    Parameters
    ----------
    s : float
        Evolution parameter s [fm^4].
    H_evolved : 1-D ndarray
        Evolving Hamiltonian which is a vector and function of s [fm^-2].
        
    Returns
    -------
    dH_vector : 1-D ndarray
        Derivative with respect to s of the evolving Hamiltonian which 
        is a vector [fm^-2].

    """
        
    # Matrix form of the evolving Hamiltonian
    H_evolved_matrix = vector2matrix(H_evolved)

    # Wegner SRG generator, eta = [G,H] where G = H_D (diagonal of the 
    # evolving Hamiltonian)
    eta = commutator( np.diag( np.diag(H_evolved_matrix) ), H_evolved_matrix)
            
    # RHS of flow equation in matrix form
    dH_matrix = commutator(eta, H_evolved_matrix)
        
    # Returns vector form of RHS of flow equation
    dH_vector = matrix2vector(dH_matrix)
        
    return dH_vector


# --- Evolve and compare methods --- #

# Set final \lambda value
lambda_final = 6.0

# Method 1: Use ode and flow parameter \lambda

# Set-up ODE

# Initial Hamiltonian as a vector and dictionary
H_initial = matrix2vector(H_matrix)

# # Use SciPy's ode function to solve flow equation
# solver = ode(derivative_method_1)
# # Following the example in Hergert:2016iju with modifications to nsteps
# # and error tolerances
# solver.set_integrator('vode', method='bdf', order=5)
# # Set initial value of Hamiltonian at lambda = lambda_initial
# solver.set_initial_value(H_initial, 20.0)
    
# t0 = time.time() # Start time
# # Solve ODE up to lamb and store in dictionary
# while solver.successful() and solver.t > lambda_final:
            
#     # Select step-size depending on extent of evolution
#     if solver.t >= 6.0:
#         dlamb = 1.0
#     elif solver.t < 6.0 and solver.t >= 2.5:
#         dlamb = 0.5
#     elif solver.t < 2.5 and solver.t >= lambda_final:
#         dlamb = 0.1
                
#     # This if statement prevents the solver from over-shooting 
#     # lambda and takes a step in lambda equal to the exact amount
#     # necessary to reach the specified lambda value
#     if solver.t - dlamb < lambda_final:
                
#         dlamb = solver.t - lambda_final
                
#     # Integrate to next step in lambda
#     H_evolved_vec = solver.integrate(solver.t - dlamb)

# t1 = time.time() # End time            

# # Print time method 1 took
# mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
# print('_'*85)
# print( 'Method 1 complete after %f minutes' % mins )
# print('_'*85)
 
# # Calculate eigenvalues
# H_final = vector2matrix(H_evolved_vec)
# eigenvalues_1 = la.eig(H_final * hbar_sq_over_m)[0]


# # Method 2: Use solve_ivp and flow parameter s

# # Set final s value
# s_final = 1.0 / lambda_final**4.0
# s_span = (0.0, s_final)

# t0 = time.time() # Start time
# solution = solve_ivp(derivative_method_2, s_span, H_initial, method='BDF',
#                       t_eval=np.array( [s_final] ))
# t1 = time.time() # End time 

# # Print time method 2 took
# mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
# print('_'*85)
# print( 'Method 2 complete after %f minutes' % mins )
# print('_'*85)

# # Calculate eigenvalues
# H_evolved_vec = solution.y[:, -1]
# H_final = vector2matrix(H_evolved_vec)
# eigenvalues_2 = la.eig(H_final * hbar_sq_over_m)[0]


# Method 3: Use solve_ivp and flow parameter \lambda

# Set final \lambda value
lambda_span = (20.0, lambda_final)

t0 = time.time() # Start time
solution = solve_ivp(derivative_method_1, lambda_span, H_initial,
                     method='BDF', t_eval=np.array( [lambda_final] ))
t1 = time.time() # End time 

# Print time method 3 took
mins = round( (t1 - t0) / 60.0, 2) # Minutes elapsed evolving H(s)
print('_'*85)
print( 'Method 3 complete after %f minutes' % mins )
print('_'*85)

# Calculate eigenvalues
H_evolved_vec = solution.y[:, -1]
H_final = vector2matrix(H_evolved_vec)
eigenvalues_3 = la.eig(H_final * hbar_sq_over_m)[0]