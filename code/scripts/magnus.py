#!/usr/bin/env python3

"""
File: magnus.py

Author: A. J. Tropiano (tropiano.4@osu.edu)
Date: June 3, 2019

The Magnus class evolves potentials to band-diagonal or block-diagonal
decoupled form with respect to the flow parameter \lambda [fm^-1], and 
possibly \Lambda_BD using the Magnus expansion implementation of the SRG.
This class is a sub-class of the SRG class from srg.py. In the Magnus 
framework, we solve the ODE for \Omega(s) using the first-order Euler step-
method, and directly apply transformations to evolve the Hamiltonian:
    
    H(s) = e^\Omega(s) H(0) e^-\Omega(s).
    
Note, tried to solve with respect to \lambda similar to SRG codes but kept
getting infinity errors in computing \Omega matrix. Thus, we evaluate with
respect to the flow parameter s, which has worked before.

Last update: June 1, 2022

"""

# To-do: Add print_info optional argument and save during first \lambda loop
# in magnus_evolve (make Potential.save_potential() take H not V).
# To-do: Might be a better way of evaluating B_k and k!
# To-do: Could maybe clean-up looping over s (make neater)?

# Python imports
from math import factorial
import numpy as np
import numpy.linalg as la
# from scipy.linalg import expm  # Option to use expm instead of BCH formula
from sympy import bernoulli
import time

# Imports from A.T. codes
from .srg import SRG
from .tools import convert_number_to_string


class Magnus(SRG):
    """
    Evolves potentials to band-diagonal or block-diagonal decoupled form with
    respect to the flow parameter \lambda [fm^-1], and possibly \Lambda_BD
    using the Magnus expansion implementation of the SRG.
    
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
        
        # Call SRG class given the potential specifications and SRG generator
        super().__init__(kvnn, channel, kmax, kmid, ntot, generator)
        
        # Initialize factorial and Bernoulli number arrays for summations up
        # to 30 terms
        self.factorial_array = np.zeros(31)
        bernoulli_array = np.zeros(31)
        for i in range(31):
            self.factorial_array[i] = factorial(i)
            bernoulli_array[i] = bernoulli(i)
        self.magnus_factors = bernoulli_array / self.factorial_array
    
    def bch(self, H_initial, O_matrix, k_max):
        """
        Evolved Hamiltonian evaluated using the BCH formula.

        Parameters
        ----------
        H_initial : 2-D ndarray
            Initial Hamiltonian matrix [fm^-2].
        O_matrix : 2-D ndarray
            Evolved \Omega matrix [unitless].
        k_max : int
            BCH sum from 0 to k_max.
        
        Returns
        -------
        H_matrix : 2-D ndarray
            Evolved Hamiltonian matrix [fm^-2].
            
        """
        
        # Initial nested commutator ad_\Omega^0(H)
        ad_matrix = H_initial
        
        # k = 0 term
        H_matrix = ad_matrix / self.factorial_array[0]
        
        # Sum from k = 1 to k = k_truncate
        for k in range(1, k_max+1):
            
            ad_matrix = self.commutator(O_matrix, ad_matrix)
            H_matrix += ad_matrix / self.factorial_array[k]
            
        return H_matrix
    
    def O_deriv(self, O_matrix, return_eta_norm=False):
        """
        Right-hand side of the Magnus \Omega(s) equation.
        
        Parameters
        ----------
        O_matrix : 2-D ndarray
            Evolved \Omega matrix [unitless].
        return_eta_norm : bool, optional
            Option to return the Frobenius norm of \eta(s).
        
        Returns
        -------
        dO_matrix : 2-D ndarray
            Derivative with respect to s of the evolving \Omega matrix [fm^-4].
        eta_norm : float, optional
            Can additionally return the Frobenius norm of \eta(s) [fm^-4].
            
        Notes
        -----
        We could use the BCH formula to evaluate H(s), though it doesn't
        change anything.

        """
        
        # Compute the evolving Hamiltonian with the BCH formula
        H_matrix = self.bch(self.H_initial, O_matrix, 25)
        # Use scipy.linalg.expm to exponentiate the matrices
        # H_matrix = expm(O_matrix) @ self.H_initial @ expm(-O_matrix)
        
        # Get SRG generator \eta = [G, H]
        eta_matrix = self.eta(H_matrix)
        
        # Initial nested commutator ad_\Omega^0(eta)
        ad_matrix = eta_matrix
        
        # k = 0 term
        dO_matrix = self.magnus_factors[0] * ad_matrix
        
        # Sum from k = 1 to k = k_max
        for k in range(1, self.k_max+1):
            
            ad_matrix = self.commutator(O_matrix, ad_matrix)
            dO_matrix += self.magnus_factors[k] * ad_matrix
    
        if return_eta_norm:
            eta_norm = la.norm(eta_matrix)
            return dO_matrix, eta_norm
        else:
            return dO_matrix
    
    def euler_method(self, s_initial, s_final, ds, O_initial):
        """
        Solves the Magnus \Omega(s) equation using the first-order Euler
        method.
        
        Parameters
        ----------
        s_initial : float
            Initial s value [fm^4].
        s_final : float
            Final s value [fm^4].
        ds : float
            Step-size in the flow parameter s [fm^4].
        O_initial : 2-D ndarray
            Initial \Omega matrix [unitless].
            
        Returns
        -------
        O_matrix : 2-D ndarray
            Evolved \Omega matrix [unitless].
            
        """

        # Set initial values
        s = s_initial
        O_matrix = O_initial
        
        # Step in s until s_final is reached
        while s <= s_final:

            # Next step in s
            O_matrix += self.O_deriv(O_matrix) * ds

            # Step to next s value
            s += ds
                
        # To ensure omega stops at s = s_final step-size is fixed to go from 
        # last value of s < s_final to s_final where ds_exact = s_final-(s-ds)
        ds_exact = s_final - s + ds
        O_matrix += self.O_deriv(O_matrix) * ds_exact
        
        return O_matrix
    
    def magnus_evolve(
            self, lambda_array, lambda_bd_array=None, k_max=6, ds=1e-5,
            save=False):
        """
        Evolve the Hamiltonian using the Magnus expansion to each value of
        \lambda, and possibly \Lambda_BD for block-diagonal decoupling. Here
        we solve for \Omega(s) using the first-order Euler method.

        Parameters
        ----------
        lambda_array : 1-D ndarray
            SRG evolution parameters \lambda [fm^-1].
        lambda_bd_array : 1-D ndarray, optional
            \Lambda_BD values for block-diagonal generator [fm^-1].
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max.
        ds : float, optional
            Step-size in the flow parameter s [fm^4].
        save : bool, optional
            If true, saves the evolved potentials.
            
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each \lambda value (and possibly an additional 
            key for \Lambda_BD). E.g., d[1.5] returns the evolved Hamiltonian
            at \lambda = 1.5 fm^-1.
            
        Notes
        -----
        Could implement an adaptive step-size in s. I'll leave ds as an
        argument to the method euler_method for this reason.
            
        """
        
        # Start flow from s = 0
        s_initial = 0.0
        
        # Initial \Omega matrix is entirely zero
        O_initial = np.zeros((self.Ntot, self.Ntot))
        
        # Evaluate 0 to k_max terms in d\Omega(s)/ds sum
        self.k_max = k_max
        
        # Start time
        t0 = time.time()
    
        # Evolve the Hamiltonian to each value of \lambda and store in a
        # dictionary (loop over \Lambda_BD as well for block-diagonal)
        d = {}
        if self.generator == 'Block-diag':
            
            for lambda_bd in lambda_bd_array:
                
                # Set the projection operators P and Q
                self.set_projection_operators(lambda_bd)
                
                # Set first key as \Lambda_BD
                d[lambda_bd] = {}
                
                for lamb in lambda_array:
            
                    # Convert lamb [fm^-1] to s value [fm^4]
                    s_final = 1.0 / lamb**4.0
            
                    # Solve for \Omega at s_final using the Euler method
                    # Automatically take smaller steps in s for large \lambda
                    if lamb >= 10.0 and ds > 1e-6:
                        O_matrix = self.euler_method(s_initial, s_final, 1e-6,
                                                     O_initial)
                    else:
                        O_matrix = self.euler_method(s_initial, s_final, ds,
                                                     O_initial)

                    # Compute the evolving Hamiltonian with the BCH formula
                    H_matrix = self.bch(self.H_initial, O_matrix, 25)
                    # Use scipy.linalg.expm to exponentiate the matrices
                    # H_matrix = (expm(O_matrix) @ self.H_initial
                    #             @ expm(-O_matrix))
            
                    # Store evolved Hamiltonian matrix in dictionary
                    d[lambda_bd][lamb] = H_matrix
                    
                    # Set starting point for next \lambda value
                    s_initial = s_final
                    O_initial = O_matrix
                    
        # Band-diagonal generators
        else:
    
            for lamb in lambda_array:
                
                # Convert lamb [fm^-1] to s value [fm^4]
                s_final = 1.0 / lamb**4.0

                # Solve for \Omega at s_final using the Euler method
                # Automatically take smaller steps in s for large \lambda
                if lamb >= 10.0 and ds > 1e-6:
                    O_matrix = self.euler_method(s_initial, s_final, 1e-6,
                                                 O_initial)
                else:
                    O_matrix = self.euler_method(s_initial, s_final, ds,
                                                 O_initial)

                # Compute the evolving Hamiltonian with the BCH formula
                H_matrix = self.bch(self.H_initial, O_matrix, 25)
                # Use scipy.linalg.expm to exponentiate the matrices
                # H_matrix = (expm(O_matrix) @ self.H_initial
                #             @ expm(-O_matrix))

                # Store evolved Hamiltonian matrix in dictionary
                d[lamb] = H_matrix

                # Set starting point for next \lambda value
                s_initial = s_final
                O_initial = O_matrix
        
        # End time
        t1 = time.time()
        
        # Print details
        mins = round((t1-t0)/60.0, 4)  # Minutes elapsed evolving H(s)
        print("_"*85)
        lamb_str = convert_number_to_string(lambda_array[-1])
        print(f"Done evolving to final \lambda = {lamb_str} fm^-1 after"
              f" {mins:.4f} minutes.")
        print("_"*85)
        print("\nSpecifications:\n")
        print(f"kvnn = {self.kvnn:d}, channel = {self.channel}")
        print(f"kmax = {self.kmax:.1f}, kmid = {self.kmid:.1f}, "
              f"ntot = {self.ntot:d}")
        print(f"method = Magnus, generator = {self.generator}")
        if self.generator == 'Block-diag':
            lambda_bd_str = convert_number_to_string(lambda_bd_array[-1])
            print(f"Final \Lambda_BD = {lambda_bd_str} fm^-1")
        print(f"k_max = {k_max:d}, ds = {ds:.1e}")
        
        # Save evolved potentials
        if save:
            
            # Get relative kinetic energy and convert to [fm^-2]
            T_matrix = self.load_kinetic_energy() / SRG.hbar_sq_over_m

            if self.generator == 'Block-diag':
                
                # Additionally loop over \Lambda_BD
                for lambda_bd in lambda_bd_array:
                    for lamb in lambda_array:
                        
                        # Scattering units here [fm^-2]
                        H_matrix = d[lambda_bd][lamb]
                    
                        # Subtract off kinetic energy [fm^-2]
                        V_matrix = H_matrix - T_matrix
                    
                        # Save evolved potential in units [fm]
                        # For large \lambda, Magnus-evolved potentials 
                        # automatically have ds <= 1e-6
                        if lamb >= 10.0 and ds > 1e-6:
                            self.save_potential(
                                V_matrix, 'magnus', self.generator, lamb,
                                lambda_bd, k_max, 1e-6)
                        else:
                            self.save_potential(
                                V_matrix, 'magnus', self.generator, lamb,
                                lambda_bd, k_max, ds)
                
            # Only need to loop over \lambda for band-diagonal generators
            else:
            
                for lamb in lambda_array:

                    # Scattering units here [fm^-2]
                    H_matrix = d[lamb]
                    
                    # Subtract off kinetic energy [fm^-2]
                    V_matrix = H_matrix - T_matrix
                    
                    # Save evolved potential in units [fm]
                    # For large \lambda, Magnus-evolved potentials 
                    # automatically have ds <= 1e-6
                    if lamb >= 10.0 and ds > 1e-6:
                        self.save_potential(V_matrix, 'magnus', self.generator,
                                            lamb, lambda_bd, k_max, 1e-6)
                    else:
                        self.save_potential(V_matrix, 'magnus', self.generator,
                                            lamb, lambda_bd, k_max, ds)

        return d
    
    def get_norms(self, s_array, k_max=6):
        """
        This function evolves over a linearly-spaced array keeping track of
        the norms ||\eta(s)|| and ||\Omega(s)||.
        
        Parameters
        ----------
        s_array : 1-D ndarray
            Linearly-spaced array of the SRG flow parameter s [fm^4]. This
            array must start from s = 0.
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max.

        Returns
        -------
        eta_norms_array : 1-D ndarray
            Frobenius norms of \eta(s) [fm^-4].
        O_norms_array : 1-D ndarray
            Frobenius norms of \Omega(s) [unitless].

        """
        
        # Evaluate 0 to k_max terms in d\Omega(s)/ds sum
        self.k_max = k_max
        
        ntot = len(s_array)
        eta_norms_array = np.zeros(ntot)
        O_norms_array = np.zeros(ntot)
        
        # Initial \Omega matrix is entirely zero
        O_matrix = np.zeros((self.Ntot, self.Ntot))
        
        # Linearly-spaced s_array means ds is constant
        ds = s_array[1]-s_array[0]

        # Step in s until s_final is reached
        for i, s in enumerate(s_array):

            # Next step in s
            try:
                
                dO_matrix, eta_norm = self.O_deriv(O_matrix,
                                                   return_eta_norm=True)
                O_matrix += dO_matrix * ds
                eta_norms_array[i] = eta_norm
                O_norms_array[i] = la.norm(O_matrix)
                
            # Check for OverflowError (\Omega(s) -> \infty)
            except OverflowError:
                
                print("_"*85)
                print("Infinities or NaNs encountered in \Omega(s).\n"
                      "Try using a smaller step-size.")
                return eta_norms_array[:i], O_norms_array[:i]
            
            except ValueError:
                
                print("_"*85)
                print("Infinities or NaNs encountered in \Omega(s).\n"
                      "Try using a smaller step-size.")
                return eta_norms_array[:i], O_norms_array[:i]

        return eta_norms_array, O_norms_array
    
    
class MagnusSplit(Magnus):
    """
    Evolves Hamiltonian to band-diagonal, decoupled form with flow parameter
    s using the Magnus "split thing" implementation. The split thing method
    involves solving the Euler method equation for the difference in \Omega
    matrices,

        \delta \Omega = Omega(s_{i+1})-Omega(s_i) = d\Omega/ds * ds.
        
    We take the matrix delta_omega and use the BCH formula to update the
    evolving Hamiltonian at each step in s,

        H(s_{i+1}) = exp^(delta_omega) * H(s_i) * exp^(-delta_omega).

    Since Omega(s_{i+1})-Omega(s_i) is directly proportional to ds, we can
    tune ds to a small enough value to prevent delta_omega from diverging to
    infinity. This method prevents the Magnus evolution from running into NaN
    or infinity errors.
    
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

    Notes
    -----
    * Evolving kvnn = 902 or 6 takes an enormous amount of time but it is
      possible using this method.
    * This class only works with band-diagonal generators! Will raise an
      error with the block-diagonal generator.

    """
    
    def __init__(self, kvnn, channel, kmax, kmid, ntot, generator):
        """Loads the initial Hamiltonian and other relevant operators depending
        on the specifications of the potential and the SRG generator.
        """
        
        # Call Magnus class given the potential specifications and SRG 
        # generator
        super().__init__(kvnn, channel, kmax, kmid, ntot, generator)
        
        # This class will not work with this generator
        if self.generator == 'Block-diag':
            raise RuntimeError('Invalid generator. '
                               'Please specify a band-diagonal generator.')
    
    def delta_O_deriv(self, delta_O_matrix, H_matrix):
        """
        Right-hand side of the Magnus \Omega(s) equation using the "split
        thing" approach.
        
        Parameters
        ----------
        delta_O_matrix : 2-D ndarray
            Difference in evolved \Omega matrices from flow parameter s_{i+1}
            to s_i [unitless].
        H_matrix : 2-D ndarray
            Evolving Hamiltonian in matrix form [fm^-2].
        
        Returns
        -------
        dO_matrix : 2-D ndarray
            Derivative with respect to s of the evolving \Omega matrix [fm^-4].

        """

        # Get SRG generator \eta = [G, H]
        eta_matrix = self.eta(H_matrix)
        
        # Initial nested commutator ad_\Omega^0(eta)
        ad_matrix = eta_matrix
        
        # k = 0 term
        dO_matrix = self.magnus_factors[0] * ad_matrix
        
        # Sum from k = 1 to k = k_max
        for k in range(1, self.k_max+1):
                
            ad_matrix = self.commutator(delta_O_matrix, ad_matrix)
            dO_matrix += self.magnus_factors[k] * ad_matrix
            
        return dO_matrix
    
    def magnus_split_evolve(self, lambda_array, k_max=6, ds=1e-5):
        """
        Evolve the Hamiltonian using the Magnus expansion and "split thing"
        approach to each value of \lambda, and possibly \Lambda_BD for block-
        diagonal decoupling. Here we solve for \Omega(s) using the first-order
        Euler method.

        Parameters
        ----------
        lambda_array : 1-D ndarray
            SRG evolution parameters \lambda [fm^-1].
        lambda_bd_array : 1-D ndarray, optional
            \Lambda_BD values for block-diagonal generator [fm^-1].
        k_max : int, optional
            d\Omega(s)/ds sum from 0 to k_max.
        ds : float, optional
            Step-size in the flow parameter s [fm^4].
    
        Returns
        -------
        d : dict
            Dictionary storing each evolved Hamiltonian with keys (floats)
            corresponding to each \lambda value (and possibly an additional 
            key for \Lambda_BD). E.g., d[1.5] returns the evolved Hamiltonian
            at \lambda = 1.5 fm^-1.
        
        """

        # Initial values
        s = 0.0
        delta_O_matrix = np.zeros((self.Ntot, self.Ntot))
        H_matrix = self.H_initial
        
        # Evaluate 0 to k_max terms in d\Omega(s)/ds sum
        self.k_max = k_max
    
        # Evolve the Hamiltonian to each value of \lambda and store in a
        # dictionary (loop over \Lambda_BD as well for block-diagonal)
        d = {}
        for lamb in lambda_array:
                
            # Convert lamb [fm^-1] to s value [fm^4]
            s_final = 1.0 / lamb**4.0
            
            # Solve ODE up to s_final using the Euler method
            while s <= s_final:
                
                # Next step in s
                delta_O_matrix = (self.delta_O_deriv(delta_O_matrix, H_matrix)
                                  * ds)
                
                # Step to next s value
                s += ds
                
                # Step to the next value of evolving Hamiltonian with BCH
                H_matrix = self.bch(H_matrix, delta_O_matrix, 25)
                    
                # Check for NaN's and infinities
                if np.isnan(H_matrix).any() or np.isinf(H_matrix).any():
                    
                    print("_"*85)
                    error = "Infinities or NaNs encountered in Hamiltonian."
                    print(error)
                    print(f"s = {s:.5e}")
                    suggestion = "Try setting ds to a smaller value."
                    print(suggestion)
                    
                    return d

            # To ensure \delta \Omega stops at s = s_val step-size is fixed to
            # go from last value of s < s_max to s_max
            ds_exact = s_final - s + ds
            delta_O_matrix = (self.delta_O_deriv(delta_O_matrix, H_matrix)
                              * ds_exact)
                
            H_matrix = self.bch(H_matrix, delta_O_matrix, 25)
            
            # Store evolved Hamiltonian matrix in dictionary
            d[lamb] = H_matrix
            
            # Reset initial s value to last s_final and continue looping
            s = s_final

        return d