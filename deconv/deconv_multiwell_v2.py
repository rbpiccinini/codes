# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 08:11:03 2020

@author: bfv8
"""

import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.signal as sg

class classHist():
    """
    A class for rate and pressure history.
    """
    def __init__(self, t, q, p):
        self.t = t
        self.q = q
        self.p = p
        self.q_1 = np.zeros(len(self.q))
        self.q_1[:-1] = self.q[1:]


class classWell():
    """
    A class with well element in the deconvolution.
    Array indices are:
        i -> time,
        j -> modes,
        k -> wells.
    """        
    def __init__(self, name, well_type, hist, psi, eig, p0):
        self.name = name
        self.type = well_type
        self.hist = hist
        
        self.tau = self.hist.t - self.hist.t[0] 
        self.dtau = np.diff(self.tau)
        
        self.psi = psi
        self.set_Q(hist.q_1) # integrate with backward value
        self.set_eig(eig)
        self.p0 = p0

    
    def set_eig(self, eig):
        self.eig=eig
    
    def set_psi(self, psi):
        self.psi=psi
    
    def set_Q(self, q):
        """
        Builds the flow rate Toeplitz matrix.
        
        Qi1,i2 = q(ti1 -ti2) if i2<=i1, or
        Qmn = 0, otherwise.

        Parameters 
        ----------
        q : numpy.array(ni,)
            Flow rate array (length = ni).

        Returns
        -------
        Q : np.array(ni,ni)
            Flow rate Toeplitz matrix
        """
        
        # Builds Toepltiz matrix for flow rate       
        padding = np.zeros(q.shape[0] - 1, q.dtype)
        first_col = q
        first_row = np.r_[q[0], padding]
        self.Q = sp.linalg.toeplitz(first_col, first_row)



class classConv():
    """
    A class to evaluate multiwell pressure-rate convolution.

    Evaluates pressure-rate deconvolution for multiwell data using the
    integral transform technique by [1].

    References
    ----------
    [1] Mikhailov, Mikhail Dimitrov, and M. Necati Ozisik.
        "Unified analysis and solutions of heat and mass diffusion." (1984).

    Parameters
    ----------
    tk : numpy.array(k,)
        Time [hours].
    qk : numpy.array(k,)
        Flow rate [m^3/day].
    pk : numpy array(k,)
        Bottom-hole pressure in kgf/cm^2.
    p0 : float
        Initial pressure. Reservoir is assumed to be initially at rest.
    dt : float
        Time-step for resampling pressure and flow-rate data. An equal sampling
        frequency for pressure and flow-rate is required to compute
        the convolution integral.
    """

    def __init__(self, wells, t, ne, p0):
        self.wells = wells
        self.nw = len(wells)
        self.ne = ne
        self.nt = len(t)
        self.p0 = p0
        self.p = np.vstack([well.hist.p for well in wells]).T
        
        self.convolve_path = np.einsum('wtx,xe,ewv -> tv',
                                   self.Q(),
                                   self.D(eig, t),
                                   self.C(psis),
                                   optimize=True)
        np.einsum_path('ij,jk,kl->il', a, b, c, optimize='greedy')
        
    def D(self, eig, t):
        """
        Builds the exponential decay matrix
        D : np.array(ni,nj).
        
        D_{ij}=  \int _{t_{i-1}}^{t_{i}}e^{- \lambda _{j} \tau}d \tau

        Parameters
        ----------
        eig : numpy.array(ne,)
            Array of exponents including the null array (length = nj).
        t   : numpy.array(nt,)
            Array of exponents including the null array (length = nj).
        """
        SMALL = 1e-12
        
        D = np.zeros([len(t), len(eig)])
        D[1:,:] = (np.exp(-np.einsum('i,j->ij', t[:-1], eig))
                        -np.exp(-np.einsum('i,j->ij', t[1:], eig)))/(eig+SMALL)  
        D[:,0] = np.diff(t, prepend=0.)
        
        return D
    
    def Q(self):
        return np.array([well.Q for well in self.wells])
    
    def C(self, psis):
        psis = psis.reshape([self.nw, self.ne]).T
        r = []
        for w in psis.T:
            r.append(np.einsum('i,ij -> ij', w, psis))
        return np.einsum('ijk->jik', np.array(r))
    
    def convolve(self, t, eig, psis):
        return self.p0 - np.einsum('wtx,xe,ewv -> tv',
                                   self.Q(),
                                   self.D(eig, t),
                                   self.C(psis),
                                   optimize=True)
    
    def error(self, x):
        psi0 = x[0]
        eig = np.concatenate([[0.], x[1:self.ne-1]])
        psis = np.zeros([self.ne, self.nw])
        psis[1:,:] = x[self.ne-1:].reshape([self.nw, self.ne-1]).T
        psis = np.ravel(psis)
        
        return np.ravel(self.convolve(self.t, eig, psis) - self.p)
    
    def eig_wells_to_x0(self, eig, wells):
        psi0 = [wells[0].psi[0]]
        psis = np.ravel([well.psi[1:] for well in wells])
        return np.concatenate([psi0, eig[1:], psis])
        
        
    def deconvolve(self,
                   eig,
                   wells,
                   bounds=None,
                   method='trf',
                   xtol=1e-8,
                   ftol=1e-8,
                   gtol=1e-8):
        
        # Create inital guess for eigenvalues
        # x0 = [1/Vpct, eigs, cs]
        psis = np.ravel([well.psi for well in wells])
        x0 = self.eig_wells_to_x0(eig, wells)
        
        # bounds for eingevalues
        lb = [0.]+[0.01]*(self.ne-1)+[-5]*(self.ne-1)*self.nw
        ub = [np.inf]+[500.]*(self.ne-1)+[5]*(self.ne-1)*self.nw
        sqbounds = [lb,ub]
        print(np.array(sqbounds).shape)
        
        r = sp.optimize.least_squares(self.error,
                                        x0,
                                        bounds=sqbounds,
                                        method='trf',
                                        xtol=xtol,
                                        ftol=ftol,
                                        gtol=gtol,
                                        verbose=2)
        return r
        


