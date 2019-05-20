# 07/11/18
# Calculates the phase shift as a function of energy from Magnus/SRG evolved
# potential (coupled-channel).
# Updated 09/25/18 to manually correct +/- shifts of pi in phase shifts

import numpy as np
import numpy.linalg as LA
from scipy.interpolate import RectBivariateSpline,interp2d

class phase_shifts(object):
    
    def __init__(self,vnn,gp,gw,conv='Stapp',spline=True):
        '''Define constants, interpolate vnn (vnn should be input in units fm)'''
        
        # Save convention
        self.conv = conv
        
        # Maximum momentum in fm^-1
        self.kmax = max(gp)

        # Units are fm^-1
        self.gp = gp # Momenta array
        self.gw = gw # Gaussian weights
        # Dimension
        N = len(gp)
        self.dim = N
        
        if spline:
            # Interpolate each sub-block with RectBivariateSpline
            self.V11_func = RectBivariateSpline(gp,gp,vnn[:N,:N])
            self.V12_func = RectBivariateSpline(gp,gp,vnn[:N,N:2*N])
            self.V21_func = RectBivariateSpline(gp,gp,vnn[N:2*N,:N])
            self.V22_func = RectBivariateSpline(gp,gp,vnn[N:2*N,N:2*N])
        else:
            # Interpolate each sub-block with interp2d
            self.V11_func = interp2d(gp,gp,vnn[:N,:N])
            self.V12_func = interp2d(gp,gp,vnn[:N,N:2*N])
            self.V21_func = interp2d(gp,gp,vnn[N:2*N,:N])
            self.V22_func = interp2d(gp,gp,vnn[N:2*N,N:2*N])
            
        self.spline = spline
        
    def delta(self,e):
        '''Returns phase shift as a function of lab energy, e [MeV]'''
        
        N = self.dim # Length of each sub-block
        
        k0 = np.sqrt(e/2.0/41.47) # fm^-1
        kmax = self.kmax # fm^-1
        
        # Build u_j vector
        k_vec = self.gp # fm^-1
        w_vec = self.gw # fm^-1
        # First N elements of u_vec
        u_vec = 2.0/np.pi * ( w_vec*k_vec**2 ) / (k_vec**2-k0**2)
        # N+1 element of u_vec
        u_last = -2.0/np.pi*k0**2*( np.sum(w_vec/(k_vec**2-k0**2)) + \
        np.log((kmax+k0)/(kmax-k0))/(2.0*k0) )
        u_vec = np.append(u_vec,u_last) # Length is N+1, units fm^-1
        
        k_full = np.append(self.gp,k0)
        
        if self.spline: # RectBivariateSpline
            col_mesh,row_mesh = np.meshgrid(k_full,k_full)
            #col_mesh2,row_mesh2 = np.meshgrid(k_vec,k_vec)
        
            # Append k0 to sub-blocks
            v11 = self.V11_func.ev(row_mesh,col_mesh) # Spline
            v12 = self.V12_func.ev(row_mesh,col_mesh)
            v21 = self.V21_func.ev(row_mesh,col_mesh)
            v22 = self.V22_func.ev(row_mesh,col_mesh)
            
        else: # interp2d
        
            v11 = np.zeros([N+1,N+1])
            v12 = np.zeros([N+1,N+1])
            v21 = np.zeros([N+1,N+1])
            v22 = np.zeros([N+1,N+1])
            
            for i in range(N+1):
                for j in range(N+1):
                    
                    v11[i,j] = self.V11_func(k_full[i],k_full[j])
                    v12[i,j] = self.V12_func(k_full[i],k_full[j])
                    v21[i,j] = self.V21_func(k_full[i],k_full[j])
                    v22[i,j] = self.V22_func(k_full[i],k_full[j])
            
        V_matrix = np.vstack((np.hstack((v11,v12)),np.hstack((v21,v22))))
        # Build A matrix, N+1 x N+1, unitless
        A_matrix = np.identity(2*(N+1)) + np.tile(u_vec,(2*(N+1),2))*V_matrix

        # Calculate R matrices and define extremes of R_matrix
        R_matrix = LA.solve(A_matrix,V_matrix) # Units fm

        R11 = R_matrix[N,N]
        R12 = R_matrix[N,2*N+1]
        # R21 = R12
        R22 = R_matrix[2*N+1,2*N+1]

        # Coupled-channel variables
        
        eps = 0.5*np.arctan(2.0*R12/(R11-R22))
        R_eps = (R11-R22)/np.cos(2.0*eps)
        delta_a = -np.arctan(0.5*k0*(R11+R22+R_eps))
        delta_b = -np.arctan(0.5*k0*(R11+R22-R_eps))
            
        # Restrict values on phases
        while delta_a-delta_b <= 0:
            delta_a += np.pi
        while delta_a-delta_b > np.pi/2.0:
            delta_b += np.pi
        
        if self.conv == 'Blatt':
            if delta_b > 0.0: # Manually fix +/- shifts in pi
                delta_b -= np.pi
            # Return phase shifts in degrees
            return np.array([delta_a,delta_b,eps])*180.0/np.pi
            
        else: # Stapp
        
            # Return eps_bar, delta_bar_a, delta_bar_b
            eps_bar = 0.5*np.arcsin(np.sin(2.0*eps)*np.sin(delta_a-delta_b))
            delta_bar_a = 0.5*(delta_a+delta_b+\
            np.arcsin(np.tan(2.0*eps_bar)/np.tan(2.0*eps)))
            delta_bar_b = 0.5*(delta_a+delta_b-\
            np.arcsin(np.tan(2.0*eps_bar)/np.tan(2.0*eps)))
            
            if delta_b > 0.0: # Manually fix +/- shifts in pi
                delta_bar_b -= np.pi
                eps_bar *= -1.0
                
            while delta_bar_a < -100.0*np.pi/180.0:
            #while delta_bar_a < 0.0:
                delta_bar_a += np.pi
            if e > 120.0:
                ang = 80.0*np.pi/180.0
            else:
                ang = np.pi
            while delta_bar_a > ang:
                delta_bar_a -= np.pi
                
            # Return phase shifts in degrees
            return np.array([delta_bar_a,delta_bar_b,eps_bar])*180.0/np.pi
            
    def delta_array(self,e_array):
        '''Calculates array of delta_bar_a phase shifts for given energies (in MeV)'''
    
        deltas = np.zeros(len(e_array))
        i = 0
        for ie in e_array:
        
            deltas[i] = self.delta(ie)[0]
            # Manually correct shifts in +/- pi
            if i > 0:
                if deltas[i] < (deltas[i-1]-100.0):
                    deltas[i] += 180.0
                if deltas[i] > (deltas[i-1]+100.0):
                    deltas[i] -= 180.0
            i += 1

        return deltas  