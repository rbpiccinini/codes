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
        self.set_D(eig)
    
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
    
    def Cs(self, wells):
        """
        Builds the eigenfunction vector.
        
        C_{j,k}=  \Psi_j(x_w)\Psi_j(x_k)

        Parameters
        ----------
        wells : list of classWell objects.
        
        Returns
        -------
        C : np.array(ni,)
            Coefficients vector
        """
        C = np.array([well.psi for well in wells])
        return np.einsum('i,ji->ij', self.psi, np.array(C))
    
    def Qs(self, wells):
        """
        Builds the multiwell flow rate Toeplitz matrix.
        
        Qi1,i2,k = q_w(ti1 -ti2) if i2<=i1, or
        Qi1,i2,k = 0, otherwise.

        Parameters 
        ----------
        wells : list of classWell objects.

        Returns
        -------
        Qs : np.array(ni, ni, nk)
            Flow rate multiwell Toeplitz matrix
        """
        Qs = np.array([well.Q for well in wells])
        return np.einsum('kij', Qs)

    def set_D(self, eig):
        """
        Builds the exponential decay matrix
        D : np.array(ni,nj).
        
        D_{ij}=  \int _{t_{i-1}}^{t_{i}}e^{- \lambda _{j} \tau}d \tau

        Parameters
        ----------
        eig : numpy.array(nj,)
            Array of exponents including the null array (length = nj).
        """
        
        # b = eig[1:]
        # ti = t[1:]
        # ti_1 = t[:-1]
        # self.D = np.zeros([len(t), len(eig)])
        # self.D[1:,1:] = (np.exp(-np.einsum('i,j->ij', ti_1, b))
        #                 -np.exp(-np.einsum('i,j->ij', ti  , b)))/b
        # self.D[1:,0] = np.diff(t)
        
        # self.D = self.D.T
        
        b = eig[1:]
        self.D = np.zeros([len(eig), len(self.tau)])
        self.D[1:,1:] = (np.exp(-np.einsum('i,j->ij', b, self.tau[:-1]))
                        -np.exp(-np.einsum('i,j->ij', b, self.tau[1:])))/b[:, np.newaxis]  
        self.D[0,1:] = self.dtau
        
        
        # self.D = np.diff(sp.integrate.cumtrapz(y=np.exp(np.einsum('i,j->ij',-eig,t)),
        #                                             x=t,
        #                                             initial=0.),
        #                   prepend=0.)
    
    def convolve(self, wells=None):
        """
        Computes the convolution integral with a list of wells.
        
        Parameters
        ----------
        wells : list of well objects
            A list of objects
        
        Returns
        -------
        dp : np.array
            Returns dp(t) = p0 - p(t)
        """
        if wells==None:
            wells = [self]
        # Qs = self.Qs(wells)
        # D  = self.D
        # C  = self.Cs(wells)
        # DC = np.einsum('ji,jk -> ik', D, C)
        # QDC = np.einsum('ijk,jk -> i', Qs, DC)
        
        QDC = np.einsum('ijk,jk -> i', self.Qs(wells), 
                        np.einsum('ji,jk -> ik', self.D, self.Cs(wells)),
                        optimize=True)
        return QDC
    
    def bourdet(self, t=None):
        """
        Computes the bourdet derivative.
        
        Parameters
        ----------
        t : numpy.array(n,)
            Times
        
        Returns
        -------
        tdp/dt : np.array
            Returns tdp(t)/dt
        """  
        if t==None:
            t=self.hist.t
        
        D = np.exp(-np.einsum('i,j->ij', t, self.eig))
        return t*np.einsum('ij,j -> i', D, self.psi**2)
    
    def PI(self, t):
        """
        Computes the well productivity index (PI).
        
        Parameters
        ----------
        t : numpy.array(n,)
            Times
        
        Returns
        -------
        PI : np.array
            Returns PI
        """
        dec = 1.0-np.exp(-np.einsum('i,j->ij', t, self.eig[1:]))
        den = np.einsum('ij,j->i', dec, self.psi[1:]**2/self.eig[1:]) 
        return 1.0/den


class classDeconv():
    """
    A class to evaluate multiwell pressure-rate deconvolution.

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

    def __init__(self, wells):
        self.wells = wells
        self.nw = len(wells)
        self.nexp = len(wells[0].eig)-1

    def find_psi(self, psi):
         
        return sp.optimize.lsq_linear(self.convMatrix(b), self.dpk,
                                      bounds=[0, np.inf],
                                      lsq_solver='exact',
                                      tol=1e-8,
                                      verbose=0).x
        # return sp.optimize.nnls(self.convMatrix(b)*self.W,
        #                         self.dpk*np.ravel(self.W))[0]
        
    def find_eig(self, eig):
        print(eig)
        print(eig.shape)
        eig = eig.reshape((self.nw, len(eig)//self.nw))
        
        for well, eigi in zip(self.wells, eig.tolist()):
            well.set_eig(np.concatenate(([0.], eigi)))
        
        x0 = np.concatenate([x.psi for x in self.wells ], axis=0)
        r = sp.optimize.least_squares(self.find_psi,
                                            x0,
                                            args=(wells),
                                            bounds=sqbounds,
                                            method='trf',
                                            xtol=xtol,
                                            ftol=ftol,
                                            gtol=gtol,
                                            verbose=2,
                                            jac=self.jac)

    def error_wells(self, wells):
        """
        doc
        """
        s = []
        for well in wells:
            s.append((well.p0 - well.convolve(wells) - well.hist.p)**2)
        return np.array(s)
    
    def error(self, x0):
        """
        doc
        """
        self.x0_to_wells(x0)
        s = []
        for well in self.wells:
            s.append((well.p0 - well.convolve(self.wells) - well.hist.p))
        return np.ravel(np.array(s))

    def jac(self, x0):
        b=sqb**2
        a = self.find_a(b)
        M = self.convMatrix(b)[:,1:]
        Mpinv = np.linalg.pinv(M, rcond=self.tol)
        Mjac = self.convMatrixJac(b)[:,1:]
        return -np.matmul(np.dot((np.eye(self.k)-np.dot(M, Mpinv)), Mjac),
                          np.diag(a[1:]))#*self.W
    
    def wells_to_x0(self):
        eigs = self.wells[0].eig[1:]
        psi0 = self.wells[0].psi[0]
        psis = np.ravel([x.psi[1:] for x in self.wells])
        x0 = np.concatenate(([psi0], eigs, psis), axis=0)
        return x0
    
 
    def x0_to_wells(self, x0):
        eigs = np.concatenate(([0.], x0[1:self.nexp+1]), axis=0)
        psi0 = x0[0]
        psis = list(x0[1+self.nexp:].reshape((self.nw, self.nexp )))
        for well,psi in zip(self.wells, psis):
            well.set_eig(eigs)
            well.set_psi(np.concatenate(([psi0], psi), axis=0))

    def deconv(self, bounds=None, tol=1e-5, method='trf',
               xtol=1e-8,
               ftol=1e-8,
               gtol=1e-8):
        
        # Create inital guess for eigenvalues
        # x0 = [1/Vpct, eigs, cs]
        x0 = self.wells_to_x0()
        
        # bounds for eingevalues
        lb = [0.]+[0.01]*(self.nexp)+[-5]*self.nexp*self.nw
        ub = [np.inf]+[500.]*(self.nexp)+[5]*self.nexp*self.nw
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
