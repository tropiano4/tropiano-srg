#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: momentum_distribution_quantum_numbers.py

Author: A. J. Tropiano (atropiano@anl.gov)
Date: May 23, 2023

Functions that return the set of all possible combinations of s.p. and partial
wave channel quantum numbers in each term of the nuclear momentum distribution
using SRG transformations.

Last update: May 23, 2023

"""

# Python imports
import numpy as np

# Imports from scripts
from scripts.tools import coupled_channel

from clebsch_gordan import clebsch_gordan_coefficient_v1 as cg_func


class PartialWaveChannel:
    """
    Partial wave channel class. Packs together the quantum numbers of a partial
    wave channel into one object.
    
    Parameters
    ----------
    channel : str
        Name of the channel (e.g., '1S0').
    
    """

    
    def __init__(self, channel):
    
        # Set instance attributes
        self.channel = channel
        self.L, self.Lp = self.get_orbital_angular_momentum()
        self.J = self.get_angular_momentum()
        self.S = self.get_spin()
        self.T = self.get_isospin(self.L, self.S)

        
    def __eq__(self, channel):

        if (
            self.L == channel.L and self.Lp == channel.Lp
            and self.J == channel.J and self.S == channel.S
            and self.T == channel.T
        ):
            
            return True
        
        else:
            
            return False
    
    
    def get_orbital_angular_momentum(self):
        """Gets the total orbital angular momentum L and L'."""
        
        if self.channel[1] == 'S':
            L = 0
        elif self.channel[1] == 'P':
            L = 1
        elif self.channel[1] == 'D':
            L = 2
        elif self.channel[1] == 'F':
            L = 3
        elif self.channel[1] == 'G':
            L = 4
        elif self.channel[1] == 'H':
            L = 5
        else:
            raise RuntimeError("Channel L exceeds the range of the function.")
            
        # Coupled-channel
        if coupled_channel(self.channel[:3]):
            
            if self.channel[5] == 'S':
                Lp = 0
            elif self.channel[5] == 'P':
                Lp = 1
            elif self.channel[5] == 'D':
                Lp = 2
            elif self.channel[5] == 'F':
                Lp = 3
            elif self.channel[5] == 'G':
                Lp = 4
            elif self.channel[5] == 'H':
                Lp = 5
            else:
                raise RuntimeError(
                    "Channel L' exceeds the range of the function.")
        
        # L = L' if the channel is not coupled
        else:
            
            Lp = L
            
        return L, Lp
    
    
    def get_angular_momentum(self):
        """Total angular momentum J = 0, 1, 2, ..."""
        
        J = int(self.channel[2])
        
        return J
    
    
    def get_spin(self):
        """Total spin S = 0 or 1."""
        
        S = int((int(self.channel[0])-1)/2)
        
        return S
    
    
    def get_isospin(self, L, S):
        """Total isospin according to antisymmetry."""
    
        # Make sure [1-(-1)^(L+S+T)] is not zero
        if (1 - (-1) ** (L+S)) == 0:
            T = 1
        else:
            T = 0
        
        return T


def kronecker_delta(x, y):
    """Kronecker \delta function: \delta_{x,y}."""
    
    return int(x == y)


def get_delta_U_quantum_numbers(tau, occupied_states, channels, cg_table):
    """Returns a list of every combination of quantum numbers in the \delta U
    linear terms that gives a non-zero product of Clebsch-Gordan coefficients.
    """
    
    spins = np.array([1/2, -1/2])
    combinations = []
    
    # Partial wave channels
    for channel in channels:
        
        # Get quantum numbers of channel
        pwc = PartialWaveChannel(channel)
        L = pwc.L
        Lp = pwc.Lp
        S = pwc.S
        J = pwc.J
        T = pwc.T
        
        # 1 - (-1)^(L+S+T) factor
        lst_factor = 1-(-1)**(L+S+T)
        # 1 - (-1)^(L'+S+T) factor
        lpst_factor = 1-(-1)**(Lp+S+T)

        # Spin projections \sigma_1 and \sigma_2
        for sigma_1 in spins:
            for sigma_2 in spins:
                
                # < \sigma_1 \sigma_2 | S M_S >
                M_S = sigma_1 + sigma_2
                spin_12_cg = cg_func(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                     cg_table)
                
                # Spin projections \sigma and \sigma'
                for sigma in spins:
                    for sigmap in spins:
                        
                        # < S M_S' | \sigma \sigma' >
                        M_Sp = sigma + sigmap
                        spin_ssp_cg = cg_func(1/2, sigma, 1/2, sigmap, S, M_Sp,
                                              cg_table)
                        
                        # Isospin projections \tau_1 and \tau_2
                        for tau_1 in spins:
                            for tau_2 in spins:
                                
                                # < \tau_1 \tau_2 | T M_T >
                                M_T = tau_1 + tau_2
                                isospin_12_cg = cg_func(1/2, tau_1, 1/2, tau_2,
                                                        T, M_T, cg_table)
                                
                                # Isospin projection \tau'
                                for taup in spins:
                                    
                                    # < T M_T | \tau \tau' >
                                    M_Tp = tau + taup
                                    if M_Tp == M_T:
                                        isospin_ttp_cg = cg_func(
                                            1/2, tau, 1/2, taup, T, M_T,
                                            cg_table
                                        )
                                    else:
                                        isospin_ttp_cg = 0
                                    
                                    for M_J in np.arange(-J, J+1):
                                        
                                        # Channel orbital angular momentum
                                        # projection M_L  
                                        M_L = M_J - M_S
                                        M_Lp = M_J - M_Sp
                                        
                                        # < L M_L S M_S | J M_J >
                                        lsj_cg = cg_func(L, M_L, S, M_S, J, M_J,
                                                         cg_table)
                                        
                                        # < J M_J | L' M_L' S M_S' >
                                        lpsj_cg = cg_func(Lp, M_Lp, S, M_Sp, J,
                                                          M_J, cg_table)
                                        
                                        # Single-particle states
                                        for alpha in occupied_states:
                                            
                                            # \delta_{\tau_1, \tau_\alpha}
                                            chi_1a = kronecker_delta(tau_1,
                                                                     alpha.m_t)
                                            
                                            # \delta_{\tau, \tau_\alpha}
                                            chi_ta = kronecker_delta(tau,
                                                                     alpha.m_t)
                                            
                                            # \delta_{\tau', \tau_\alpha}
                                            chi_tpa = kronecker_delta(taup,
                                                                      alpha.m_t)
                                            
                                            # l_\alpha \sigma_1 CG
                                            alpha_1_cg = cg_func(
                                                alpha.l, alpha.m_j-sigma_1, 1/2,
                                                sigma_1, alpha.j, alpha.m_j,
                                                cg_table
                                            )
                                            
                                            # l_\alpha \sigma CG
                                            alpha_s_cg = cg_func(
                                                alpha.l, alpha.m_j-sigma, 1/2,
                                                sigma, alpha.j, alpha.m_j,
                                                cg_table
                                            )
                                            
                                            # l_\alpha \sigma' CG
                                            alpha_sp_cg = cg_func(
                                                alpha.l, alpha.m_j-sigmap, 1/2,
                                                sigmap, alpha.j, alpha.m_j,
                                                cg_table
                                            )

                                            for beta in occupied_states:
                                            
                                                # \delta_{\tau_2, \tau_\beta}
                                                chi_2b = kronecker_delta(
                                                    tau_2, beta.m_t)
                                                
                                                # \delta_{\tau', \tau_\beta}
                                                chi_tpb = kronecker_delta(
                                                    taup, beta.m_t)
                                                
                                                # \delta_{\tau, \tau_\beta}
                                                chi_tb = kronecker_delta(
                                                    tau, beta.m_t)
                                                
                                                # l_\beta \sigma_2 CG
                                                beta_2_cg = cg_func(
                                                    beta.l, beta.m_j-sigma_2,
                                                    1/2, sigma_2, beta.j,
                                                    beta.m_j, cg_table
                                                )
                                                
                                                # l_\beta \sigma' CG
                                                beta_sp_cg = cg_func(
                                                    beta.l, beta.m_j-sigmap,
                                                    1/2, sigmap, beta.j,
                                                    beta.m_j, cg_table
                                                )
                                                
                                                # l_\beta \sigma CG
                                                beta_s_cg = cg_func(
                                                    beta.l, beta.m_j-sigma, 1/2,
                                                    sigma, beta.j, beta.m_j,
                                                    cg_table
                                                )
                                                
                                                # Calculate the product of all
                                                # CG's and factors
                                                product = (
                                                    spin_12_cg
                                                    * spin_ssp_cg
                                                    * isospin_12_cg
                                                    * isospin_ttp_cg
                                                    * lsj_cg
                                                    * lpsj_cg
                                                    * lst_factor
                                                    * lpst_factor
                                                    * chi_1a * alpha_1_cg
                                                    * chi_2b * beta_2_cg
                                                    * (
                                                        chi_ta * alpha_s_cg
                                                        * chi_tpb * beta_sp_cg
                                                        - chi_tpa * alpha_sp_cg
                                                        * chi_tb * beta_s_cg
                                                    )
                                                )
                                                
                                                # Add this set of quantum
                                                # numbers to the list if the
                                                # product is nonzero
                                                if product != 0:
                                                                                                
                                                    d = {
                                                        'sigma_1': sigma_1,
                                                        'sigma_2': sigma_2,
                                                        'sigma': sigma,
                                                        'sigmap': sigmap,
                                                        'tau_1': tau_1,
                                                        'tau_2': tau_2,
                                                        'taup': taup,
                                                        'n_alpha': alpha.n,
                                                        'l_alpha': alpha.l,
                                                        'j_alpha': alpha.j,
                                                        'm_j_alpha': alpha.m_j,
                                                        'm_t_alpha': alpha.m_t,
                                                        'n_beta': beta.n,
                                                        'l_beta': beta.l,
                                                        'j_beta': beta.j,
                                                        'm_j_beta': beta.m_j,
                                                        'm_t_beta': beta.m_t,
                                                        'S': S,
                                                        'M_S': M_S,
                                                        'M_Sp': M_Sp,
                                                        'L': L,
                                                        'M_L': M_L,
                                                        'Lp': Lp,
                                                        'M_Lp': M_Lp,
                                                        'J': J,
                                                        'M_J': M_J,
                                                        'T': T,
                                                        'M_T': M_T
                                                    }
                                                                                                
                                                    combinations.append(d)
                                    
    return combinations


def get_delta_U2_quantum_numbers(tau, occupied_states, channels, cg_table):
    """Returns a list of every combination of quantum numbers in the \delta U^2
    term that gives a non-zero product of Clebsch-Gordan coefficients.
    """
    
    spins = np.array([1/2, -1/2])
    combinations = []
    
    # Partial wave channel for \delta U
    for channel_1 in channels:
        
        # Get quantum numbers of channel
        pwc_1 = PartialWaveChannel(channel_1)
        L = pwc_1.L
        Lp = pwc_1.Lp
        S = pwc_1.S
        J = pwc_1.J
        T = pwc_1.T
        
        # 1 - (-1)^(L+S+T) factor
        lst_factor = 1-(-1)**(L+S+T)
        # 1 - (-1)^(L'+S+T) factor
        lpst_factor = 1-(-1)**(Lp+S+T)
        
        # Partial wave channel for \delta U^\dagger
        for channel_2 in channels:
            
            # Get quantum numbers of channel
            pwc_2 = PartialWaveChannel(channel_2)
            Lpp = pwc_2.L
            Lppp = pwc_2.Lp
            Sp = pwc_2.S
            Jp = pwc_2.J
            Tp = pwc_2.T
            
            # Evaluate Kronecker \delta function in S and S'
            if Sp == S:

                # 1 - (-1)^(L''+S+T') factor
                lppstp_factor = 1-(-1)**(Lpp+S+Tp)
                # 1 - (-1)^(L'''+S+T') factor
                lpppstp_factor = 1-(-1)**(Lppp+S+Tp)

                # Spin projections \sigma_1 and \sigma_2
                for sigma_1 in spins:
                    for sigma_2 in spins:

                        # < \sigma_1 \sigma_2 | S M_S >
                        M_S = sigma_1 + sigma_2
                        spin_12_cg = cg_func(1/2, sigma_1, 1/2, sigma_2, S, M_S,
                                             cg_table)
                        
                        # Spin projections \sigma_3 and \sigma_4
                        for sigma_3 in spins:
                            for sigma_4 in spins:
                                
                                # < S M_S''' | \sigma_3 \sigma_4 >
                                M_Sppp = sigma_3 + sigma_4
                                spin_34_cg = cg_func(1/2, sigma_3, 1/2, sigma_4,
                                                     S, M_Sppp, cg_table)
                                
                                # Isospin projections \tau_1 and \tau_2
                                for tau_1 in spins:
                                    for tau_2 in spins:
                                        
                                        # < \tau_1 \tau_2 | T M_T >
                                        M_T = tau_1 + tau_2
                                        isospin_12_cg = cg_func(
                                            1/2, tau_1, 1/2, tau_2, T, M_T,
                                            cg_table
                                        )
                                        
                                        # Isospin projections \tau_3 and \tau_4
                                        for tau_3 in spins:
                                            for tau_4 in spins:
                                                
                                                # < T' M_T' | \tau_3 \tau_4 >
                                                M_Tp = tau_3 + tau_4
                                                isospin_34_cg = cg_func(
                                                    1/2, tau_3, 1/2, tau_4, Tp,
                                                    M_Tp, cg_table
                                                )
                                        
                                                # Isospin projection \tau'
                                                for taup in spins:
                                            
                                                    # < T M_T | \tau \tau' >
                                                    isospin_ttp_cg = cg_func(
                                                        1/2, tau, 1/2, taup, T,
                                                        M_T, cg_table
                                                    )
                                                    
                                                    # < \tau \tau' | T' M_T' >
                                                    isospin_ttpp_cg = cg_func(
                                                        1/2, tau, 1/2, taup, Tp,
                                                        M_Tp, cg_table
                                                    )
                                                    
                                                    for M_J in np.arange(-J, J+1):
                                                        
                                                        M_L = M_J - M_S
                                                        
                                                        # L S J CG
                                                        lsj_cg = cg_func(
                                                            L, M_L, S, M_S,
                                                            J, M_J, cg_table
                                                        )
                                                        
                                                        # Sum over M_S'
                                                        for M_Sp in np.arange(-S, S+1):
                                                        
                                                            M_Lp = M_J - M_Sp
                                                        
                                                            # L' S J CG
                                                            lpsj_cg = cg_func(
                                                                Lp, M_Lp, S,
                                                                M_Sp, J, M_J,
                                                                cg_table
                                                            )
                                                        
                                                            for M_Jp in np.arange(-Jp, Jp+1):
                                                            
                                                                M_Lpp = M_Jp - M_Sp
                                                                M_Lppp = M_Jp - M_Sppp
                                                                
                                                                # L'' S J' CG
                                                                lppsjp_cg = cg_func(
                                                                    Lpp, M_Lpp, S, M_Sp,
                                                                    Jp, M_Jp, cg_table
                                                                )
                                                                
                                                                # L''' S J' CG
                                                                lpppsjp_cg = cg_func(
                                                                    Lppp, M_Lppp, S, M_Sppp,
                                                                    Jp, M_Jp, cg_table
                                                                )
                                                                
                                                                # Single-particle states
                                                                for alpha in occupied_states:
                                                                    
                                                                    # \delta_{\tau_1, \tau_\alpha}
                                                                    chi_1a = kronecker_delta(tau_1, alpha.m_t)
                                                                    
                                                                    # \delta_{\tau_3, \tau_\alpha}
                                                                    chi_3a = kronecker_delta(tau_3, alpha.m_t)
                                                                    
                                                                    # \delta_{\tau_4, \tau_\alpha}
                                                                    chi_4a = kronecker_delta(tau_4, alpha.m_t)
                                                        
                                                                    # l_\alpha \sigma_1 CG
                                                                    alpha_1_cg = cg_func(
                                                                        alpha.l, alpha.m_j-sigma_1, 1/2,
                                                                        sigma_1, alpha.j, alpha.m_j, cg_table
                                                                    )
                                                        
                                                                    # l_\alpha \sigma_3 CG
                                                                    alpha_3_cg = cg_func(
                                                                        alpha.l, alpha.m_j-sigma_3, 1/2,
                                                                        sigma_3, alpha.j, alpha.m_j, cg_table
                                                                    )
                                                                    
                                                                    # l_\alpha \sigma_4 CG
                                                                    alpha_4_cg = cg_func(
                                                                        alpha.l, alpha.m_j-sigma_4, 1/2,
                                                                        sigma_4, alpha.j, alpha.m_j, cg_table
                                                                    )
                                                                    
                                                                    for beta in occupied_states:
                                                                        
                                                                        # \delta_{\tau_2, \tau_\beta}
                                                                        chi_2b = kronecker_delta(tau_2, beta.m_t)
                                                                        
                                                                        # \delta_{\tau_3, \tau_\beta}
                                                                        chi_3b = kronecker_delta(tau_3, beta.m_t)
                                                                        
                                                                        # \delta_{\tau_4, \tau_\beta}
                                                                        chi_4b = kronecker_delta(tau_4, beta.m_t)
                                                                        
                                                                        # l_\beta \sigma_2 CG
                                                                        beta_2_cg = cg_func(
                                                                            beta.l, beta.m_j-sigma_2, 1/2,
                                                                            sigma_2, beta.j, beta.m_j, cg_table
                                                                        )
                                                                        
                                                                        # l_\beta \sigma_3 CG
                                                                        beta_3_cg = cg_func(
                                                                            beta.l, beta.m_j-sigma_3, 1/2,
                                                                            sigma_3, beta.j, beta.m_j, cg_table
                                                                        )
                                                                        
                                                                        # l_\beta \sigma_4 CG
                                                                        beta_4_cg = cg_func(
                                                                            beta.l, beta.m_j-sigma_4, 1/2,
                                                                            sigma_4, beta.j, beta.m_j, cg_table
                                                                        )
                                        
                                                                        # Calculate the product of all
                                                                        # CG's and factors
                                                                        product = (
                                                                            spin_12_cg
                                                                            * spin_34_cg
                                                                            * isospin_12_cg
                                                                            * isospin_ttp_cg
                                                                            * isospin_ttpp_cg
                                                                            * isospin_34_cg
                                                                            * lsj_cg
                                                                            * lpsj_cg
                                                                            * lppsjp_cg
                                                                            * lpppsjp_cg
                                                                            * lst_factor
                                                                            * lpst_factor
                                                                            * lppstp_factor
                                                                            * lpppstp_factor
                                                                            * chi_1a * alpha_1_cg
                                                                            * chi_2b * beta_2_cg
                                                                            * (
                                                                                chi_3a * alpha_3_cg
                                                                                * chi_4b * beta_4_cg
                                                                                - chi_4a * alpha_4_cg
                                                                                * chi_3b * beta_3_cg
                                                                            )  
                                                                        )

                                                                        # Add this set of quantum
                                                                        # numbers to the list if the
                                                                        # product is nonzero
                                                                        if product != 0:
                                                                                                
                                                                            d = {
                                                                                'sigma_1': sigma_1,
                                                                                'sigma_2': sigma_2,
                                                                                'sigma_3': sigma_3,
                                                                                'sigma_4': sigma_4,
                                                                                'tau_1': tau_1,
                                                                                'tau_2': tau_2,
                                                                                'tau_3': tau_3,
                                                                                'tau_4': tau_4,
                                                                                'taup': taup,
                                                                                'n_alpha': alpha.n,
                                                                                'l_alpha': alpha.l,
                                                                                'j_alpha': alpha.j,
                                                                                'm_j_alpha': alpha.m_j,
                                                                                'm_t_alpha': alpha.m_t,
                                                                                'n_beta': beta.n,
                                                                                'l_beta': beta.l,
                                                                                'j_beta': beta.j,
                                                                                'm_j_beta': beta.m_j,
                                                                                'm_t_beta': beta.m_t,
                                                                                'S': S,
                                                                                'M_S': M_S,
                                                                                'M_Sp': M_Sp,
                                                                                'L': L,
                                                                                'M_L': M_L,
                                                                                'Lp': Lp,
                                                                                'M_Lp': M_Lp,
                                                                                'J': J,
                                                                                'M_J': M_J,
                                                                                'T': T,
                                                                                'M_T': M_T,
                                                                                'M_Sppp': M_Sppp,
                                                                                'Lpp': Lpp,
                                                                                'M_Lpp': M_Lpp,
                                                                                'Lppp': Lppp,
                                                                                'M_Lppp': M_Lppp,
                                                                                'Jp': Jp,
                                                                                'M_Jp': M_Jp,
                                                                                'Tp': Tp,
                                                                                'M_Tp': M_Tp
                                                                            }
                                                                                                
                                                                            combinations.append(d)
                                    
    return combinations


def quantum_number_array(quantum_numbers, file_name=None):
    """Returns all possible quantum numbers as an array where each row specifies
    a particular combination. If a file name is passed, then this function will
    save the array.
    """
    
    M = len(quantum_numbers)
    N = len(quantum_numbers[0])
    quantum_number_array = np.zeros((M, N))
    
    # Loop over sets of quantum numbers
    for i, d in enumerate(quantum_numbers):
        
        # Loop over quantum numbers of the set
        for j, key in enumerate(d):
            
            quantum_number_array[i, j] = d[key]
            
    if file_name is not None:
        
        np.savetxt('quantum_numbers/' + file_name, quantum_number_array,
                   fmt='%.1f')
        
    return quantum_number_array
            
        
        
    