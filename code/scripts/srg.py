#!/usr/bin/env python3

"""
File: srg.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 1, 2019

The SRG class evolves potentials to band-diagonal or block-diagonal decoupled
form with respect to the flow parameter \lambda [fm^-1], and possibly
\Lambda_BD. This class is a subclass of the Potential class from potentials.py.
The SRG class can solve either the usual flow equation

    dH(s)/ds = [\eta(s), H(s)],
    
or the differential equation for U(s) directly,

    dU(s)/ds = \eta(s) U(s).
    
Additionally, this script includes a function for computing the SRG unitary
transformation itself, given the initial and SRG-evolved Hamiltonians.

Last update: April 3, 2023

"""

# Python imports
import numpy as np
import numpy.linalg as la
from scipy.integrate import ode, solve_ivp
import time

# Imports from A.T. codes
from .potentials import Potential
from .tools import build_coupled_channel_matrix, convert_number_to_string


class SRG(Potential):
    """
    Evolves potentials to band-diagonal or block-diagonal decoupled form with
    respect to the flow parameter \lambda [fm^-1], and possibly \Lambda_BD.

    Parameters
    ----------
    kvnn : int
        This number specifies the potential.
    channel : str
        The partial wave channel (e.g. '1S0').
    kmax : float
        Maximum value in the momentum mesh [fm^-1].
    kmid : float
        Mid-point value in the momentum mesh [fm^-1].
    ntot : int
        Number of momentum points in mesh.
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
            
    """

    def __init__(self, kvnn, channel, kmax, kmid, ntot, generator):
        """Loads the initial Hamiltonian and other relevant operators depending
        on the specifications of the potential and the SRG generator.
        """

        # Call Potential class given the potential specifications
        super().__init__(kvnn, channel, kmax, kmid, ntot)

        # Get initial Hamiltonian associated with the potential
        H_initial_MeV = self.load_hamiltonian()  # [MeV]

        # Convert Hamiltonian to scattering units [fm^-2]
        self.H_initial = H_initial_MeV / Potential.hbar_sq_over_m

        # Set length of Hamiltonian
        self.Ntot = len(self.H_initial)

        # Get relative kinetic energy
        if generator == 'T':

            T_rel_MeV = self.load_kinetic_energy()  # [MeV]

            # Convert to scattering units [fm^-2]
            self.T_rel = T_rel_MeV / Potential.hbar_sq_over_m

        # Need relative momenta for construction of projection operators for
        # the block-diagonal generator
        elif generator == 'Block-diag':

            # Get momentum array (don't worry about weights) in [fm^-1]
            self.k_array, _ = self.load_mesh()
            self.ntot = len(self.k_array)

        # Set generator for evaluation of \eta
        self.generator = generator

    def set_projection_operators(self, lambda_bd):
        """
        Sets sub-block projection operators for block-diagonal generator where
        P_matrix (Q_matrix) corresponds to the low (high) sub-block of the
        Hamiltonian.

        Parameters
        ----------
        lambda_bd : float
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].

        """

        # Use booleans to partition the Hamiltonian
        bool_array = self.k_array < lambda_bd

        # Number of points in the momentum mesh
        ntot = self.ntot

        # Matrix of ones along the diagonal up to k > lambda_bd
        P_matrix = np.diag(np.ones(ntot) * bool_array)

        # Opposite of P
        Q_matrix = np.identity(ntot) - P_matrix

        # Projection operators for coupled-channel potentials
        if ntot != self.Ntot:
            zeros = np.zeros((ntot, ntot))
            P_matrix = build_coupled_channel_matrix(P_matrix, zeros, zeros,
                                                    P_matrix)
            Q_matrix = build_coupled_channel_matrix(Q_matrix, zeros, zeros,
                                                    Q_matrix)

        # Save both operators for evaluations of \eta
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
        n = int(N * (N + 1) / 2)

        # Initialize vectorized matrix
        v = np.zeros(n)

        # Algorithm for reshaping M to the vector v
        i = 0
        j = N
        for k in range(N):
            v[i:j] = M[k][k:]
            i = j
            j += N - k - 1

        return v

    def vector_to_matrix(self, v):
        """
        Takes the vector of an upper triangle matrix v and returns the full
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
            j += N - k - 1

        # Now reflect the upper half to lower half to build full matrix
        # M.T - np.diag(np.diag(M)) is the lower half of M excluding diagonal
        return M + (M.T - np.diag(np.diag(M)))

    def commutator(self, A, B):
        """Commutator of square matrices A and B."""

        return A @ B - B @ A

    def eta(self, H_matrix):
        """
        SRG generator \eta = [G, H] where G is specified by the decoupling
        scheme.
        
        Parameters
        ----------
        H_matrix : 2-D ndarray
            Evolving Hamiltonian in matrix form [fm^-2].
        
        Returns
        -------
        output : 2-D ndarray
            SRG generator \eta in matrix form [fm^-4].
        
        """

        # G = H_D (diagonal of the evolving Hamiltonian)
        if self.generator == 'Wegner':

            G_matrix = np.diag(np.diag(H_matrix))

        # G = T_rel (relative kinetic energy)
        elif self.generator == 'T':

            G_matrix = self.T_rel

        # G = H_BD (block-diagonal of the evolving Hamiltonian)
        elif self.generator == 'Block-diag':

            G_matrix = (self.P_matrix @ H_matrix @ self.P_matrix
                        + self.Q_matrix @ H_matrix @ self.Q_matrix)

        # \eta = [G, H]
        return self.commutator(G_matrix, H_matrix)

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
            which is also a vector [fm^-1].

        """

        # Matrix form of the evolving Hamiltonian
        H_matrix = self.vector_to_matrix(H_vector)

        # Get SRG generator \eta = [G, H]
        eta_matrix = self.eta(H_matrix)

        # RHS of the flow equation in matrix form
        dH_matrix = -4.0 / lamb ** 5 * self.commutator(eta_matrix, H_matrix)

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
            is also a vector [fm].

        """

        # Matrix form of the SRG transformation
        U_matrix = np.reshape(U_vector, (self.Ntot, self.Ntot))

        # Evolve the Hamiltonian to compute \eta
        H_matrix = U_matrix @ self.H_initial @ U_matrix.T

        # Get SRG generator \eta = [G, H]
        eta_matrix = self.eta(H_matrix)

        # RHS in matrix form
        dU_matrix = -4.0 / lamb ** 5 * eta_matrix @ U_matrix

        # Returns vector form of RHS
        dU_vector = np.reshape(dU_matrix, -1)

        return dU_vector

    def select_step_size(self, solver_lambda, lambda_final):
        """
        Select ODE solver step-size depending on the extent of evolution.
        We can take bigger steps at large values of \lambda.

        Parameters
        ----------
        solver_lambda : float
            The ODE solver's interal value of \lambda [fm^-1].
        lambda_final : float
            The external value of \lambda [fm^-1].

        Returns
        -------
        dlamb : float
            Step-size in terms of \lambda [fm^-1].

        """
            
        if solver_lambda >= 6.0:
            dlamb = 1.0
        elif 2.5 <= solver_lambda < 6.0:
            dlamb = 0.5
        elif 1.5 <= solver_lambda < 2.5:
            dlamb = 0.1
        else:
            dlamb = 0.05

        return dlamb

    def get_ode_solver(self, lambda_initial, method='hamiltonian'):
        """Sets up the ODE solver with SciPy's integrate.ode function."""

        # Solving for H(s)
        if method == 'hamiltonian':

            solver = ode(self.H_deriv)

            # Initial Hamiltonian as a vector
            H_initial = self.matrix_to_vector(self.H_initial)

            # Set initial conditions
            solver.set_initial_value(H_initial, lambda_initial)

        # Solving for U(s)
        elif method == 'srg_transformation':

            solver = ode(self.U_deriv)

            # Set initial value: U = identity matrix and make vector
            U_initial = np.eye(self.Ntot).reshape(-1)
            solver.set_initial_value(U_initial, lambda_initial)

        # Print an error message if method is invalid
        else:

            raise RuntimeError("Need to specify a valid method.")

        # Following the example in Hergert:2016iju with modifications to
        # nsteps and error tolerances
        # solver.set_integrator('vode', method='bdf', order=5, atol=1e-10,
        #                       rtol=1e-10, nsteps=5000000)
        solver.set_integrator('lsoda', atol=1e-10, rtol=1e-10, nsteps=5000000)

        return solver

    def srg_evolve(
            self, lambda_array, lambda_bd_array=None, lambda_initial=20.0,
            method='hamiltonian', save=False):
        """
        Evolve the Hamiltonian at each value of \lambda, and possibly
        \Lambda_BD for block-diagonal decoupling.
        
        Parameters
        ----------
        lambda_array : 1-D ndarray
            SRG evolution parameters \lambda [fm^-1].
        lambda_bd_array : 1-D ndarray, optional
            \Lambda_BD values for block-diagonal generator [fm^-1].
        lambda_initial : float, optional
            Initial value of lambda [fm^-1]. Technically this should be
            infinity but a large value (~20 fm^-1) is sufficient.
        save : bool, optional
            If true, saves data files.
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
            
        """
        
        # TESTING
        # Limits of \lambda
        lamb_limits = [lambda_initial, lambda_array[-1]]

        # Initial Hamiltonian as a vector
        if method == 'hamiltonian':
            
            H_initial_vector = self.matrix_to_vector(self.H_initial)
            
        # Initial SRG transformation as a matrix
        elif method == 'srg_transformation':
            
            U_initial = np.eye(self.Ntot)

            # Initial SRG transformation as a vector
            U_initial_vector  = np.reshape(U_initial, -1)

        # Start time
        t0 = time.time()

        # Evolve the Hamiltonian (or U) to each value of \lambda and store in
        # a dictionary (loop over \Lambda_BD as well for block-diagonal)
        d = {}
        if self.generator == 'Block-diag':

            for lambda_bd in lambda_bd_array:

                # Set the projection operators P and Q
                self.set_projection_operators(lambda_bd)

                # Set first key as \Lambda_BD
                d[lambda_bd] = {}

                    # # Set-up ODE solver
                    # solver = self.get_ode_solver(lambda_initial, method)
    
                    # for lamb in lambda_array:
    
                    #     # Solve ODE up to lamb
                    #     while solver.successful() and round(solver.t, 2) > lamb:
                            
                    #         # Get ODE solver step-size in \lambda
                    #         dlamb = self.select_step_size(solver.t, lamb)
    
                    #         # Integrate to next step in \lambda
                    #         solution_vector = solver.integrate(solver.t - dlamb)
                
                
                # TESTING
                if method == 'hamiltonian':
                    
                    result = solve_ivp(
                        self.H_deriv, lamb_limits, H_initial_vector,
                        method='BDF', t_eval=lambda_array, atol=1e-10,
                        rtol=1e-10
                    )
                    
                elif method == 'srg_transformation':
                    
                    result = solve_ivp(
                        self.U_deriv, lamb_limits, U_initial_vector,
                        method='BDF', t_eval=lambda_array, atol=1e-10,
                        rtol=1e-10
                    )
                
                # for lamb in lambda_array:
                for i, lamb in enumerate(lambda_array):

                    # Store evolved Hamiltonian (or U) matrix in dictionary
                    if method == 'hamiltonian':
                        
                        # d[lambda_bd][lamb] = self.vector_to_matrix(
                        #     solution_vector
                        # )
                        d[lambda_bd][lamb] = self.vector_to_matrix(
                            result.y[:, i]
                        )
                        
                        if save:  # Save evolved potential?
                            
                            self.save_srg_potential(d[lambda_bd][lamb], lamb,
                                                    lambda_bd)
                        
                    elif method == 'srg_transformation':
                        
                        # d[lambda_bd][lamb] =  np.reshape(
                        #     solution_vector, (self.Ntot, self.Ntot)
                        # )
                        d[lambda_bd][lamb] =  np.reshape(
                            result.y[:, i], (self.Ntot, self.Ntot)
                        )
                        
                        if save:  # Save SRG transformation?
                        
                            self.save_srg_transformation(d[lambda_bd][lamb],
                                                         lamb, lambda_bd)
                        
        # Band-diagonal generators
        else:
            
            # if method == 'hamiltonian':
                    
            #     result = solve_ivp(
            #         self.H_deriv, lamb_limits, H_initial_vector, method='BDF',
            #         t_eval=lambda_array, atol=1e-10, rtol=1e-10
            #     )
                    
            # elif method == 'srg_transformation':
                    
            #     result = solve_ivp(
            #         self.U_deriv, lamb_limits, U_initial_vector, method='BDF',
            #         t_eval=lambda_array, atol=1e-10, rtol=1e-10
            #     )
                
            # for i, lamb in enumerate(lambda_array):

            #     # Store evolved Hamiltonian (or U) matrix in dictionary
            #     if method == 'hamiltonian':
                        
            #         d[lamb] = self.vector_to_matrix(result.y[:, i])
                        
            #         if save:  # Save evolved potential?
                            
            #             self.save_srg_potential(d[lamb], lamb)
                        
            #     elif method == 'srg_transformation':

            #         U_vector = result.y[:, i]
            #         d[lamb] = np.reshape(U_vector, (self.Ntot, self.Ntot))
            #         # d[lamb] = np.reshape(result.y[:, i],
            #         #                      (self.Ntot, self.Ntot))
                        
            #         if save:  # Save SRG transformation?
                        
            #             self.save_srg_transformation(d[lamb], lamb)

            # Set-up ODE solver
            solver = self.get_ode_solver(lambda_initial, method)

            for lamb in lambda_array:

                # Solve ODE up to lamb and store in dictionary
                while solver.successful() and round(solver.t, 2) > lamb:
                    
                    # Get ODE solver step-size in \lambda
                    dlamb = self.select_step_size(solver.t, lamb)

                    # Integrate to next step in lambda
                    solution_vector = solver.integrate(solver.t - dlamb)

                # Store evolved Hamiltonian (or U) matrix in dictionary
                if method == 'hamiltonian':
                    
                    d[lamb] = self.vector_to_matrix(solution_vector)
                    
                    if save:  # Save evolved potential?
                    
                        self.save_srg_potential(d[lamb], lamb)
                    
                elif method == 'srg_transformation':
                    
                    d[lamb] =  np.reshape(solution_vector,
                                          (self.Ntot, self.Ntot))
                    
                    if save:  # Save SRG transformation?
                    
                        self.save_srg_transformation(d[lamb], lamb)

        # End time
        t1 = time.time()

        # Print details
        mins = round((t1 - t0) / 60.0, 4)  # Minutes elapsed evolving H(s)
        print("_" * 85)
        lamb_str = convert_number_to_string(lambda_array[-1])
        print(f"Done evolving to final \lambda = {lamb_str} fm^-1 after"
              f" {mins:.4f} minutes.")
        print("_" * 85)
        print("\nSpecifications:\n")
        print(f"kvnn = {self.kvnn:d}, channel = {self.channel}")
        print(f"kmax = {self.kmax:.1f}, kmid = {self.kmid:.1f}, "
              f"ntot = {self.ntot:d}")
        print(f"method = SRG, generator = {self.generator}")
        if self.generator == 'Block-diag':
            lambda_bd_str = convert_number_to_string(lambda_bd_array[-1])
            print(f"Final \Lambda_BD = {lambda_bd_str} fm^-1")

        return d
    
    def save_srg_potential(self, H_matrix, lamb, lambda_bd=None):
        """Saves the SRG-evolved potential."""
        
        # Get relative kinetic energy and convert to [fm^-2]
        T_matrix = self.load_kinetic_energy() / Potential.hbar_sq_over_m
        
        # Subtract off kinetic energy [fm^-2]
        V_matrix = H_matrix - T_matrix
        
        self.save_potential(V_matrix, 'srg', self.generator, lamb, lambda_bd)
        
    def srg_transformation_file_name(self, lamb, lambda_bd=None):
        """Get file name for SRG transformation."""
        
        file_name = (f'u_{self.channel}_kvnn_{self.kvnn_string}'
                     f'_{self.generator}')

        # Get \lambda with correct number of decimals
        if lamb == round(lamb, 1):
            lamb_str = str(round(lamb, 1))
        else:
            lamb_str = str(round(lamb, 2))

        # Add \Lambda_BD to name for block-diagonal generator
        if self.generator == 'Block-diag':
            file_name += f'{lambda_bd:.2f}_lambda{lamb_str}'
        else:
            file_name += f'_lambda{lamb_str}'

        file_name += '.out'
        
        return file_name

    def save_srg_transformation(self, U_matrix, lamb, lambda_bd=None):
        """Saves the SRG transformation."""
        
        # Get file name for the SRG transformation matrix elements
        file_name = self.srg_transformation_file_name(lamb, lambda_bd)

        # Get momenta and weights
        k_array, k_weights = self.load_mesh()

        f = open(self.potential_directory + file_name, 'w')

        # Write each sub-block as a column for coupled-channel transformations
        if self.coupled_channel_bool:

            header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format(
                'k', 'kp', 'U11', 'U12', 'U21', 'U22'
            )
            f.write('#' + header + '\n')

            for i, ik in enumerate(k_array):
                for j, jk in enumerate(k_array):

                    u11 = U_matrix[i, j]
                    u12 = U_matrix[i, j+self.ntot]
                    u21 = U_matrix[i+self.ntot, j]
                    u22 = U_matrix[i+self.ntot, j+self.ntot]

                    line = ('{:^15.6f}{:^15.6f}'.format(ik, jk)
                            + '{:^23e}{:^23e}'.format(u11, u12)
                            + '{:^23e}{:^23e}'.format(u21, u22))

                    f.write(line + '\n')

        else:

            header = '{:^15s}{:^15s}{:^23s}'.format('k', 'kp', 'U')
            f.write('#' + header + '\n')

            for i, ik in enumerate(k_array):
                for j, jk in enumerate(k_array):

                    u = U_matrix[i, j]

                    line = '{:^15.6f}{:^15.6f}{:^23e}'.format(ik, jk, u)

                    f.write(line + '\n')

        f.close()
        
    def load_srg_transformation(self, lamb, lambda_bd=None):
        """
        Loads SRG unitary transformation from data file generated from
        solving dU/d\lambda =  -4/\lambda^5 \eta(\lambda) U(\lambda).
        
        Parameters
        ----------
        lamb : float
            SRG evolution parameter \lambda [fm^-1].
        lambda_bd : float, optional
            SRG \Lambda_BD value for block-diagonal generator [fm^-1].
            
        Returns
        -------
        U_matrix : 2-D ndarray
            SRG unitary transformation matrix with integration weights
            [unitless].
            
        """

        # Get file name for the SRG transformation matrix elements
        file_name = self.srg_transformation_file_name(lamb, lambda_bd)

        # Load output file
        data = np.loadtxt(self.potential_directory + file_name)

        # Coupled-channel potential?
        if self.coupled_channel_bool:

            u11 = np.reshape(data[:, 2], (self.ntot, self.ntot))
            u12 = np.reshape(data[:, 3], (self.ntot, self.ntot))
            u21 = np.reshape(data[:, 4], (self.ntot, self.ntot))
            u22 = np.reshape(data[:, 5], (self.ntot, self.ntot))
            U_matrix = build_coupled_channel_matrix(u11, u12, u21, u22)

        else:

            U_matrix = np.reshape(data[:, 2], (self.ntot, self.ntot))

        return U_matrix


def compute_srg_transformation(H_initial, H_evolved):
    """
    SRG unitary transformation built out of eigenvectors of the initial and 
    evolved Hamiltonians.
    
    Parameters
    ----------
    H_initial : 2-D ndarray
        Initial Hamiltonian matrix [MeV].
    H_evolved : 2-D ndarray
        Evolved Hamiltonian matrix [MeV].
        
    Returns
    -------
    U_matrix : 2-D ndarray
        SRG unitary transformation matrix with integration weights [unitless].
        
    """

    Ntot = len(H_initial)

    # Get the eigenvectors of the initial and SRG-evolved Hamiltonians
    _, vecs_initial = la.eigh(H_initial)
    _, vecs_evolved = la.eigh(H_evolved)

    # Initialize unitary transformation U with same size as Hamiltonians
    U_matrix = np.zeros((Ntot, Ntot))

    # Transformation is given by summing over the outer product of evolved and
    # initial eigenvectors
    for alpha in range(Ntot):

        # Individual eigenvectors (these are already sorted correctly from 
        # numpy.linalg.eigh)
        psi_alpha_initial = vecs_initial[:, alpha]
        psi_alpha_evolved = vecs_evolved[:, alpha]

        # Make sure the phases match
        if psi_alpha_initial.T @ psi_alpha_evolved < 0:
            psi_alpha_evolved = -psi_alpha_evolved

        # Outer product of eigenvectors
        U_matrix += np.outer(psi_alpha_evolved, psi_alpha_initial)

    return U_matrix