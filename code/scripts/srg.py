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

Last update: May 15, 2023

"""

# Python imports
import numpy as np
import numpy.linalg as la
from scipy.integrate import ode
import time

# Imports from scripts
from .potentials import Potential
from .tools import build_coupled_channel_matrix, convert_number_to_string


class SRG:
    """
    Evolves potentials to band-diagonal or block-diagonal decoupled form with
    respect to the flow parameter \lambda [fm^-1], and possibly \Lambda_BD.

    Parameters
    ----------
    potential : Potential
        NN potential projected onto a partial wave channel in relative momentum
        space.
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
            
    """

    
    def __init__(self, potential, generator):
        """Loads the initial Hamiltonian and other relevant operators depending
        on the potential and the SRG generator.
        """
        
        # Set potential as instance attribute
        self.potential = potential
        
        # Get initial Hamiltonian associated with the potential
        H_initial_MeV = potential.load_hamiltonian()  # [MeV]

        # Convert Hamiltonian to scattering units [fm^-2]
        self.H_initial = H_initial_MeV / Potential.hbar_sq_over_m

        # Set length of Hamiltonian
        self.Ntot = len(self.H_initial)

        # Set relative kinetic energy as instance attribute
        T_rel_MeV = potential.load_kinetic_energy()  # [MeV]
        self.T_rel = T_rel_MeV / Potential.hbar_sq_over_m  # [fm^-2]

        # Need relative momenta for construction of projection operators for
        # the block-diagonal generator
        if generator == 'Block-diag':

            # Get momentum array (don't worry about weights) in [fm^-1]
            self.k_array, _ = potential.load_mesh()
            self.ntot = len(self.k_array)

        # Set generator for evaluation of \eta
        self.generator = generator
        
        
    def __call__(
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
        
        # Start time
        t0 = time.time()

        # Evolve the Hamiltonian (or U) to each value of \lambda and store in
        # a dictionary (loop over \Lambda_BD as well for block-diagonal)
        d = {}
        
        # Set-up ODE solver
        solver = self.get_ode_solver(lambda_initial, method)
        
        # Block-diagonal generator
        if self.generator == 'Block-diag':
            
            for lambda_bd in lambda_bd_array:

                # Set the projection operators P and Q
                self.set_projection_operators(lambda_bd)

                # Set first key as \Lambda_BD
                d[lambda_bd] = {}

                for lamb in lambda_array:

                    # Solve ODE up to lamb and store in dictionary
                    while solver.successful() and round(solver.t, 2) > lamb:
                    
                        # Get ODE solver step-size in \lambda
                        dlamb = self.select_step_size(solver.t, lamb)

                        # Integrate to next step in lambda
                        solution_vector = solver.integrate(solver.t - dlamb)

                    # Store evolved Hamiltonian (or U) matrix in dictionary
                    if method == 'hamiltonian':
                    
                        d[lambda_bd][lamb] = self.vector_to_matrix(
                            solution_vector)
                    
                        if save:  # Save evolved potential?
                    
                            V_evolved = d[lambda_bd][lamb] - self.T_rel  # fm^-2
                            self.potential.save_potential(
                                V_evolved, 'srg', self.generator, lamb,
                                lambda_bd
                            )
                    
                    elif method == 'srg_transformation':
                    
                        d[lambda_bd][lamb] =  np.reshape(solution_vector,
                                                         (self.Ntot, self.Ntot))
                    
                        if save:  # Save SRG transformation?

                            save_srg_transformation(
                                d[lambda_bd][lamb], self.potential,
                                self.generator, lamb, lambda_bd
                            )
     
        # Band-diagonal generators
        else:
            
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
                    
                        V_evolved = d[lamb] - self.T_rel  # fm^-2
                        self.potential.save_potential(V_evolved, 'srg',
                                                      self.generator, lamb)
                    
                elif method == 'srg_transformation':
                    
                    d[lamb] =  np.reshape(solution_vector,
                                          (self.Ntot, self.Ntot))
                    
                    if save:  # Save SRG transformation?
                    
                        save_srg_transformation(d[lamb], self.potential,
                                                self.generator, lamb)

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
        print(f"kvnn = {self.potential.kvnn:d}")
        print(f"channel = {self.potential.channel}")
        print(
            f"kmax = {self.potential.kmax:.1f}, "
            f"kmid = {self.potential.kmid:.1f}, "
            f"ntot = {self.potential.ntot:d}"
        )
        print(f"method = SRG, generator = {self.generator}")
        if self.generator == 'Block-diag':
            lambda_bd_str = convert_number_to_string(lambda_bd_array[-1])
            print(f"Final \Lambda_BD = {lambda_bd_str} fm^-1")

        return d

    
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

        # Matrix of ones along the diagonal up to k > lambda_bd
        P_matrix = np.diag(np.ones(self.ntot) * bool_array)

        # Opposite of P
        Q_matrix = np.identity(self.ntot) - P_matrix

        # Projection operators for coupled-channel potentials
        if self.ntot != self.Ntot:
            zeros = np.zeros((self.ntot, self.ntot))
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
        solver.set_integrator('vode', method='bdf', order=5, atol=1e-10,
                              rtol=1e-10, nsteps=5000000)

        return solver


def srg_transformation_file_name(potential, generator, lamb, lambda_bd=None):
    """File name for SRG transformation."""
        
    file_name = (f"u_{potential.channel}_kvnn_{potential.kvnn_string}"
                 f"_{generator}")

    # Get \lambda with correct number of decimals
    if lamb == round(lamb, 1):
        lamb_str = str(round(lamb, 1))
    else:
        lamb_str = str(round(lamb, 2))

    # Add \Lambda_BD to name for block-diagonal generator
    if generator == 'Block-diag':
        file_name += f'{lambda_bd:.2f}_lambda{lamb_str}'
    else:
        file_name += f'_lambda{lamb_str}'

    file_name += '.out'
        
    return file_name


def save_srg_transformation(
        U_matrix, potential, generator, lamb, lambda_bd=None
):
    """Saves the SRG transformation."""
        
    # Get file name for the SRG transformation matrix elements
    file_name = srg_transformation_file_name(potential, generator, lamb,
                                             lambda_bd)

    # Get momenta and weights
    k_array, k_weights = potential.load_mesh()

    f = open(potential.potential_directory + file_name, 'w')

    # Write each sub-block as a column for coupled-channel transformations
    if potential.coupled_channel_bool:

        header = '{:^15s}{:^15s}{:^23s}{:^23s}{:^23s}{:^23s}'.format(
            'k', 'kp', 'U11', 'U12', 'U21', 'U22'
        )
        f.write('#' + header + '\n')

        for i, ik in enumerate(k_array):
            for j, jk in enumerate(k_array):

                u11 = U_matrix[i, j]
                u12 = U_matrix[i, j+potential.ntot]
                u21 = U_matrix[i+potential.ntot, j]
                u22 = U_matrix[i+potential.ntot, j+potential.ntot]

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
    
    
def load_srg_transformation(potential, generator, lamb, lambda_bd=None):
    """
    Loads SRG unitary transformation from data file generated from solving
        dU/d\lambda =  -4/\lambda^5 \eta(\lambda) U(\lambda).
        
    Parameters
    ----------
    potential : Potential
        NN potential projected onto a partial wave channel in relative momentum
        space.
    generator : str
        SRG generator 'Wegner', 'T', or 'Block-diag'.
    lamb : float
        SRG evolution parameter \lambda [fm^-1].
    lambda_bd : float, optional
        SRG \Lambda_BD value for block-diagonal generator [fm^-1].
            
    Returns
    -------
    U_matrix : 2-D ndarray
        SRG unitary transformation matrix with integration weights [unitless].
            
    """

    # Get file name for the SRG transformation matrix elements
    file_name = srg_transformation_file_name(potential, generator, lamb,
                                             lambda_bd)

    # Load output file
    data = np.loadtxt(potential.potential_directory + file_name)

    # Coupled-channel potential?
    if potential.coupled_channel_bool:

        u11 = np.reshape(data[:, 2], (potential.ntot, potential.ntot))
        u12 = np.reshape(data[:, 3], (potential.ntot, potential.ntot))
        u21 = np.reshape(data[:, 4], (potential.ntot, potential.ntot))
        u22 = np.reshape(data[:, 5], (potential.ntot, potential.ntot))
        U_matrix = build_coupled_channel_matrix(u11, u12, u21, u22)

    else:

        U_matrix = np.reshape(data[:, 2], (potential.ntot, potential.ntot))

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