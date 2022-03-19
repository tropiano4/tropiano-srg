#!/usr/bin/env python3

"""
File: srg.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: May 1, 2019

Evolves Hamiltonian to band-diagonal or block-diagonal decoupled form with
flow parameter \lambda [fm^-1]. This class relies on the Potential class from
vnn.py. The function can solve either the usual flow equation
    dH(s)/ds = [\eta(s), H(s)],
or the differential equation for U(s) directly,
    dU(s)/ds = \eta(s) U(s).

Last update: March 18, 2022

"""

# To-do: Implement solve_ivp in evolve() method.

# Python imports
import numpy as np
from scipy.integrate import ode

# Imports from A.T. codes
from .tools import build_coupled_channel_matrix


class SRG:
    
    def __init__(self, potential, generator):
        """
        Loads the initial Hamiltonian and other relevant operators depending
        on the SRG generator.
        
        Parameters
        ----------
        potential : Potential
            Potential object from vnn.py. Contains information and useful
            potential-related methods.
        generator : str
            SRG generator 'Wegner', 'T', or 'Block-diag'.
            
        """
        
        # Get initial Hamiltonian associated with the input potential
        H_initial_MeV = potential.load_hamiltonian()  # [MeV]
        
        # Convert Hamiltonian to scattering units [fm^-2] where hbar^2/M = 1
        self.H_initial = H_initial_MeV / potential.hbar_sq_over_m
        
        # Save length of Hamiltonian
        self.Ntot = len(self.H_initial)
        
        # Get relative kinetic energy
        if generator == 'T':
            
            T_rel_MeV = potential.load_kinetic_energy()  # [MeV]
            
            # Convert to scattering units [fm^-2]
            self.T_rel = T_rel_MeV / potential.hbar_sq_over_m
        
        # Need relative momenta for construction of projection operators for
        # the block-diagonal generator
        elif generator == 'Block-diag':
            
            # Get momentum array (don't worry about weights) in [fm^-1]
            self.k_array, _ = potential.load_mesh()
            self.ntot = len(self.k_array)
            
        # Save generator for evaluation of \eta
        self.generator = generator
    
    def set_projection_operators(self, lambda_bd):
        """
        Sets sub-block projection operators for block-diagonal generator where
        P_matrix (Q_matrix) corresponds to the low (high) sub-block of the
        Hamiltonian. Both operators have the same shape as the Hamiltonian.

        Parameters
        ----------
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].

        """
        
        # Use booleans to partition the Hamiltonian
        bool_array = self.k_array < lambda_bd
        
        ntot = self.ntot

        # Matrix of ones along the diagonal up to k > lambda_bd
        P_matrix = np.diag(np.ones(ntot)*bool_array)
        
        # Opposite of p
        Q_matrix = np.identity(ntot) - P_matrix
        
        # Projection operators for coupled-channel potentials
        if ntot != self.Ntot:
            
            zeros = np.zeros( (ntot, ntot) )
            P_matrix = build_coupled_channel_matrix(P_matrix, zeros, zeros,
                                                    P_matrix)
            Q_matrix = build_coupled_channel_matrix(Q_matrix, zeros, zeros,
                                                    Q_matrix)
        
        # Save both operators
        self.P_matrix = P_matrix
        self.Q_matrix = Q_matrix
        
    def matrix_to_vector(self, M):
        """
        Takes the upper triangle of the matrix M (including the diagonal) 
        and reshapes it into a vector v.
        
        Parameters
        ----------
        M : 2-D ndarray
            Input matrix of shape (N, N).
            
        Returns
        -------
        v : 1-D ndarray
            Output vector of shape (N*(N+1)/2,).
            
        """
    
        # Length of matrix
        N = len(M)
        # Length of vectorized matrix
        n = int(N*(N+1)/2)
        
        # Initialize vectorized matrix
        v = np.zeros(n)
    
        # Algorithm for reshaping M to the vector v
        i = 0
        j = N
        for k in range(N):
            v[i:j] = M[k][k:]
            i = j
            j += N-k-1

        return v
    
    def vector_to_matrix(self, v):
        """
        Takes the vector of a upper triangle matrix v and returns the full 
        matrix M. Use only for symmetric matrices.
        
        Parameters
        ----------
        v : 1-D ndarray
            Input vector of shape (N*(N+1)/2,).
        
        Returns
        -------
        A : 2-D ndarray
            Output matrix of shape (N, N).
            
        """
        
        # Length of matrix (we know the length from the Hamiltonian)
        N = self.Ntot
    
        # Initialize matrix
        M = np.zeros((N, N))
    
        # Build the upper half of M with the diagonal included

        # Algorithm for reshaping v to the matrix M
        i = 0
        j = N
        for k in range(N):
            M[k, k:] = v[i:j]
            i = j
            j += N-k-1

        # Now reflect the upper half to lower half to build full matrix
        # M.T - np.diag(np.diag(M)) is the lower half of M excluding diagonal
        return M + (M.T-np.diag(np.diag(M)))
    
    def commutator(self, A, B):
        """Commutator of square matrices A and B."""
        
        return A @ B - B @ A
    
    def H_deriv(self, lamb, H_vector):
        """
        Right-hand side of the SRG flow equation.
        
        Parameters
        ----------
        lamb : float
            SRG evolution parameter \lambda [fm^-1].
        H_vector : 1-D ndarray
            Evolving Hamiltonian [fm^-2].
        
        Returns
        -------
        dH_vector : 1-D ndarray
            Derivative with respect to \lambda of the evolving Hamiltonian
            which is also a vector in scattering units [fm^-2].

        """
        
        # Matrix form of the evolving Hamiltonian
        H_matrix = self.vector_to_matrix(H_vector)

        # Compute \eta = [G,H] which is dependent on SRG generator
        if self.generator == 'Wegner':
        
            # G = H_D (diagonal of the evolving Hamiltonian)
            G_matrix = np.diag(np.diag(H_matrix))
            
        elif self.generator == 'T':
            
            # G = T_rel
            G_matrix = self.T_rel
            
        elif self.generator == 'Block-diag':
            
            # G = H_BD (block-diagonal of the evolving Hamiltonian)
            G_matrix = (self.P_matrix @ H_matrix @ self.P_matrix
                        + self.Q_matrix @ H_matrix @ self.Q_matrix)

        eta_matrix = self.commutator(G_matrix, H_matrix)
            
        # RHS of the flow equation in matrix form
        dH_matrix = -4.0 / lamb**5 * self.commutator(eta_matrix, H_matrix)
        
        # Returns vector form of RHS of flow equation
        dH_vector = self.matrix_to_vector(dH_matrix)
        
        return dH_vector
    
    def U_deriv(self, lamb, U_vector):
        """
        Right-hand side of differential equation for U(\lambda)
            dU(\lambda)/d\lambda = -4/\lambda^5 \eta(\lambda) U(\lambda).
        
        Parameters
        ----------
        lamb : float
            SRG evolution parameter \lambda [fm^-1].
        U_vector : 1-D ndarray
            SRG transformation [unitless] as a function of \lambda.
        
        Returns
        -------
        dU_vector : 1-D ndarray
            Derivative with respect to \lambda of the SRG transformation which
            is also a unitless vector.

        """

        # Matrix form of the SRG transformation
        U_matrix = np.reshape(U_vector, (self.Ntot, self.Ntot))
        
        # Evolve the Hamiltonian to compute \eta
        H_matrix = U_matrix @ self.H_initial @ U_matrix.T

        # Compute \eta = [G,H] which is dependent on SRG generator
        if self.generator == 'Wegner':
        
            # G = H_D (diagonal of the evolving Hamiltonian)
            G_matrix = np.diag(np.diag(H_matrix))
            
        elif self.generator == 'T':
            
            # G = T_rel
            G_matrix = self.T_rel
            
        elif self.generator == 'Block-diag':
            
            # G = H_BD (block-diagonal of the evolving Hamiltonian)
            G_matrix = (self.P_matrix @ H_matrix @ self.P_matrix
                        + self.Q_matrix @ H_matrix @ self.Q_matrix)

        eta_matrix = self.commutator(G_matrix, H_matrix)
             
        # RHS in matrix form
        dU_matrix = -4.0 / lamb**5 * eta_matrix @ U_matrix
        
        # Returns vector form of RHS
        dU_vector = np.reshape(dU_matrix, -1)
        
        return dU_vector
    
    def evolve(
            self, lambda_initial, lambda_array, lambda_bd_array=np.empty(0),
            method='hamiltonian'):
        """
        Evolve the Hamiltonian at each value of \lambda, and possibly
        \Lambda_BD for block-diagonal decoupling.
        
        Parameters
        ----------
        lambda_initial : float
            Initial value of lambda [fm^-1]. Technically this should be
            infinity but a large value (~20 fm^-1) is sufficient.
        lambda_array : 1-D ndarray
            SRG evolution parameters \lambda [fm^-1].
        lambda_bd_array : 1-D ndarray, optional
            \Lambda_BD values for block-diagonal generator [fm^-1].
        method : str, optional
            The default method is to solve the flow equation
                dH(s)/ds = [\eta(s), H(s)],
            specified by method='hamiltonian'. Alternatively, one can set
            method='srg_transformation' to solve for U(s) directly:
                dU(s)/ds = \eta(s) U(s).
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian (or SRG transformation)
            with keys (floats) corresponding to each \lambda value (and
            possibly an additional key for \Lambda_BD). E.g., d[1.5] returns
            the evolved Hamiltonian at \lambda = 1.5 fm^-1.
            
        Notes
        -----
        Could implement the scipy.integrate.solve_ivp(...) ODE solver.
            
        """

        # Store solutions of the ODE in the following dictionary
        d = {}
        
        # Initial Hamiltonian as a vector
        H_initial = self.matrix_to_vector(self.H_initial)
        
        # Set-up ODE with SciPy's integrate.ode function
        if method == 'hamiltonian':
            
            solver = ode(self.H_deriv)
            
            # Set initial value of Hamiltonian at \lambda = lambda_initial
            solver.set_initial_value(H_initial, lambda_initial)
            
        elif method == 'srg_transformation':
            
            solver = ode(self.U_deriv)
            
            # Set initial value: U = identity matrix and make vector
            U_initial = np.eye(self.Ntot).reshape(-1)
            solver.set_initial_value(U_initial, lambda_initial)
        
        # Print an error message if method is invalid
        else:
            
            print('Need to specify a valid method.')
            return None
            
        # Following the example in Hergert:2016iju with modifications to
        # nsteps and error tolerances
        solver.set_integrator('vode', method='bdf', order=5, atol=1e-6,
                              rtol=1e-6, nsteps=5000000)

        # Evolve the Hamiltonian (or U) to each value of \lambda and store in
        # the dictionary (loop over \Lambda_BD as well for block-diagonal)
        if self.generator == 'Block-diag':
            
            for lambda_bd in lambda_bd_array:
                
                # Set the projection operators P and Q
                self.set_projection_operators(lambda_bd)
                
                # Set first key as \Lambda_BD
                d[lambda_bd] = {}
                
                for lamb in lambda_array:
                    
                    # Solve ODE up to lamb
                    while solver.successful() and solver.t > lamb:
                        
                        # Select step-size depending on extent of evolution
                        if solver.t >= 6.0:
                            dlamb = 1.0
                        elif solver.t < 6.0 and solver.t >= 2.5:
                            dlamb = 0.5
                        elif solver.t < 2.5 and solver.t >= lamb:
                            dlamb = 0.1
                            
                        # This if statement prevents the solver from over-
                        # shooting \lambda and takes a step in \lambda equal
                        # to the exact amount necessary to reach the final
                        # \lambda value
                        if solver.t - dlamb < lamb:
                
                            dlamb = solver.t - lamb
                            
                        # Integrate to next step in \lambda
                        solution_vector = solver.integrate(solver.t - dlamb)
                
                # Store evolved Hamiltonian (or U) matrix in dictionary
                d[lambda_bd][lamb] = self.vector_to_matrix(solution_vector)
        
        # Band-diagonal generators
        else:
            
            for lamb in lambda_array:
            
                # Solve ODE up to lamb and store in dictionary
                while solver.successful() and solver.t > lamb:
            
                    # Select step-size depending on extent of evolution
                    if solver.t >= 6.0:
                        dlamb = 1.0
                    elif solver.t < 6.0 and solver.t >= 2.5:
                        dlamb = 0.5
                    elif solver.t < 2.5 and solver.t >= lamb:
                        dlamb = 0.1
                
                    # See previous comment on this
                    if solver.t - dlamb < lamb:
                
                        dlamb = solver.t - lamb
                
                    # Integrate to next step in lambda
                    solution_vector = solver.integrate(solver.t - dlamb)
                
                # Store evolved Hamiltonian (or U) matrix in dictionary
                d[lamb] = self.vector_to_matrix(solution_vector)
            
        return d