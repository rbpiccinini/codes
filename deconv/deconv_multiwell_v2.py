# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 08:11:03 2020

@author: bfv8
"""
import numpy as np
import scipy as sp
import pandas as pd
import mpmath as mp

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

    """        
    def __init__(self, name, well_type, hist, psi, eig, p0):
        self.name = name
        self.type = well_type
        self.hist = hist
        
        self.tau = self.hist.t - self.hist.t[0] 
        self.dtau = np.diff(self.tau)
        self.nt = len(self.hist.t)
        
        self.psi = psi
        self.set_Qnu(hist.q_1) # integrate with backward value
        self.set_eig(eig)
        self.p0 = p0
        self.pconv=None

    
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
        
        
    def set_Qnu(self, q):
        """
        Builds the flow rate Toeplitz matrix for irregularly spaced data.
        
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
        q_interp = sp.interpolate.interp1d(self.hist.t, q, kind='linear')
        nt = len(self.hist.t)
        self.Q = np.zeros([nt, nt])
        
        dt = np.tril(np.tile(self.hist.t, (nt, 1)).T - np.tile(self.hist.t, (nt, 1)))
        self.Q = q_interp(dt)

    def bourdet(self, t=None, well=None):
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
        if t is None:
            t=self.hist.t
        if well is None:
            psib=self.psi
        else:
            psib=well.psi
        
        D = np.exp(-np.einsum('i,j->ij', t, self.eig))
        return t*np.einsum('ij,j -> i', D, self.psi*psib)
    
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
    """

    def __init__(self, wells, t, ne, p0):
        self.wells = wells
        self.nw = len(wells)
        self.ne = ne
        self.nt = len(t)
        self.p0 = p0
        self.p = np.vstack([well.hist.p for well in wells]).T
        self.t = t
        self.conv_path=True
        
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
        
        # Analytical calculation
        SMALL = 1e-12
        
        D = np.zeros([len(t), len(eig)])
        
        D[1:,:] = (np.exp(-np.einsum('i,j->ij', t[:-1], eig))
                   *(1.0-np.exp(-np.einsum('i,j->ij', t[1:]-t[:-1], eig))))/(eig+SMALL)  
        
        D[:,0] = np.diff(t, prepend=0.)
        
        # Numerical calculation
        # D = np.diff(sp.integrate.cumtrapz(y=np.exp(np.einsum('i,j->ij',-eig,t)),
        #                                             x=t,
        #                                             initial=0.),
        #                   prepend=0.).T
        return D
    
    def Q(self):
        return np.array([well.Q for well in self.wells])
    
    def C(self, psis):
        psis = psis.reshape([self.nw, self.ne]).T
        r = []
        for w in psis.T:
            r.append(np.einsum('i,ij -> ij', w, psis))
        return np.einsum('ijk->jik', np.array(r))
    
    def convolve(self, t, eig, psis, optimize=True):
        return self.p0 - np.einsum('wtx,xe,ewv -> tv',
                                   self.Q(),
                                   self.D(eig, t),
                                   self.C(psis),
                                   optimize=optimize)
    
    def error(self, x, optimize=True):
        eig, psis = self.x2conv(x)
        self.x = x
        return (np.ravel(self.convolve(self.t, eig, psis, optimize) - self.p)
                /np.ravel(self.p))*self.p0
    
    def x2conv(self, x):
        psi0 = x[0]
        eig = np.zeros(self.ne)
        eig[1:] = x[1:self.ne]
        psis = np.zeros([self.nw, self.ne])
        psis[:, 1:] = x[self.ne:].reshape([self.nw, self.ne-1])
        psis[:,0] = psi0
        
        return eig, np.ravel(psis)
        
    def conv2x(self, eig, psis):
        x = np.zeros(len(eig)+(self.ne-1)*self.nw)
        x[0] = psis[0]
        
        psis = psis.reshape([self.nw, self.ne])
        psis = np.ravel(psis[:, 1:])

        x[1:self.ne] = eig[1:]
        x[self.ne:] = psis
        
        return x
    
    def set_wells(self, eig, psis):
        psis = psis.reshape([self.nw, self.ne])
        pconvs = self.convolve(self.t, eig, psis).T
        for well, psi, pconv in zip(self.wells, psis, pconvs):
            well.psi = psi
            well.eig = eig
            well.pconv = pconv
    
    def get_df(self):
        dfs = []
        keys = []
        for well in self.wells:
            df = pd.DataFrame(columns=['eig', 'psi'])
            df['eig'] = well.eig
            df['psi'] = well.psi
            dfs.append(df)
            keys.append(well.name)
        return pd.concat(dfs, keys=keys, names=['well'])
        
    def deconvolve(self,
                   eig,
                   psis,
                   bounds=None,
                   method='trf',
                   xtol=1e-6,
                   ftol=1e-7,
                   gtol=1e-5):
        
        # Create inital guess for eigenvalues
        # x0 = [1/Vpct, eigs, cs]
        
        # Einsum optimal path not improving speed??
        # self.conv_path = np.einsum_path('wtx,xe,ewv -> tv',
        #                                 self.Q(),
        #                                 self.D(eig, self.t),
        #                                 self.C(psis),
        #                                 optimize='optimal')[0]
        self.conv_path=True
        x0 = self.conv2x(eig, psis)       
        r = sp.optimize.least_squares(self.error,
                                        x0,
                                        bounds=bounds,
                                        method='trf',
                                        jac='3-point',
                                        xtol=xtol,
                                        ftol=ftol,
                                        gtol=gtol,
                                        max_nfev=500,
                                        verbose=2,
                                        kwargs={'optimize':self.conv_path})
        return r
        



