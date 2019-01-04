import numpy as np
import pandas as pd
import scipy.optimize as sp

class classPVT:
    """Class for fluid PVT properties.

    Parameters
    ----------
    pvt_table : string
        Excel xlsx file with PVT properties table.
    Pb : float
        Bubble-point pressure (kgf/cm2).
    co : float
        Undersaturated oil compressibility (cm2/kgf).
    cw : float
        Water compressibility (cm2/kgf).
    T : float
        Reservoir fluid temperature (K).

    Attributes
    -------
    Bo : float
        Returns Bo for given oil phase pressure and bubble-point
        pressure (m3/sm3).
    Bg : float
        Returns Bg for given gas phase pressure (m3/sm3).
    Bw : float
        Returns Bw for given water phase pressure (m3/sm3).
    Rs : float
        Returns gas-oil solubility ratio (sm3/sm3).
    """

    def __init__(self, pvt_table, co, cw, T):
        # reading pvt table
        self.pvt_table = pd.read_excel(pvt_table)
        self.pe = self.pvt_table['p']
        self.Boe = self.pvt_table['Bo']
        self.Rse = self.pvt_table['Rs']
        self.Bge = self.pvt_table['Bg']
        self.Bwe = self.pvt_table['Bw']
        # reading other properties
        self.co = co
        self.cw = cw
        self.T = T
        # computing gas compressibility factor Z
        self.Bgi = 1.033/self.pe*(self.T)/288.
        self.Ze = self.Bge/self.Bgi

    def Bo(self, p, Pb):
        """Returns Bo for given bubble-point and oil phase pressure.

        Parameters
        ----------
        p : float
            Oil phase pressure (kgf/cm2).
        Pb : float
            Oil phase bubble-point pressure (kgf/cm2).

        Returns
        -------
        Bo : float
            Formation-volume factor of oil phase (m3/sm3).

        """

        if p >= Pb:
            return np.interp(Pb, self.pe[::-1], self.Boe[::-1]) \
                    * (1.-self.co*(p-Pb))
        else:
            return np.interp(p, self.pe[::-1], self.Boe[::-1])

    def Bg(self, p):
        """Returns Bg for given gas phase pressure.

        Parameters
        ----------
        p : float
            Gas phase pressure (kgf/cm2).

        Returns
        -------
        Bg : float
            Formation-volume factor of gas phase (m3/sm3).

        """
        Z=np.interp(p,self.pe[::-1],self.Ze[::-1])
        return 1.033/p*(self.T)/288.*Z

    def Bw(self,p):
        """Returns Bw for given gas phase pressure.

        Parameters
        ----------
        p : float
            Water phase pressure (kgf/cm2).

        Returns
        -------
        Bw : float
            Formation-volume factor of water phase (m3/sm3).

        """
        return np.interp(p,self.pe[::-1],self.Bwe[::-1])

    def Rs(self,p,Pb):
        """Returns Rs for given oil-phase pressure.

        Parameters
        ----------
        p : float
            Water phase pressure (kgf/cm2).
        Pb : float
            Oil phase bubble-point pressure (kgf/cm2).

        Returns
        -------
        Rs: float
            Gas-oil solubility ratio (m3/m3).

        """
        if p >= Pb:
            return np.interp(Pb,self.pe[::-1],self.Rse[::-1])
        else:
            return np.interp(p,self.pe[::-1],self.Rse[::-1])



class classAquifer:
    """Class for aquifer properties and modeling.

    Parameters
    ----------
    Waq : float
        Aquifer volume (m3).
    p0 : float
        Initial aquifer pressure (kgf/cm2).
    cr : float
        pore volume compressibility (cm2/kgf).
    pvt : classPVT object
        Fluid PVT properties.

    Attributes
    -------
    calcWe : float
        Returns aquifer influx for a given final pressure.

    """
    def __init__(self, Waq, p0, cr, pvt):
        self.pvt = pvt
        self.Waq = Waq
        self.p0 = p0
        self.cr = cr

    def calcWe(self,p):
        """Returns aquifer influx for a given final pressure.

        Returns aquifer influx for a given final pressure assuming the simple
        tank model.

        Parameters
        ----------
        p : float
            Final aquifer pressure (kgf/cm2).

        Returns
        -------
        We : float
            Volume of aquifer influx (m3).

        """
        ct=self.cr+self.pvt.cw
        We=ct*self.Waq*(self.p0-p)
        return We

class classRes:
    """Class for reservoir properties.

    Parameters
    ----------
    N : float
        Stock-tank original oil in place - STOOIP (sm3).
    p0 : float
        Initial reservoir pressure (kgf/cm2).
    Pb0 : float
        Initial bubble-point pressure of oil phase (kgf/cm2).
    cr : float
        Pore volume compressibility (cm2/kgf).
    Swcon : float
        Connate water saturation.
    m : float
        Ratio of initial gas-cap volume to initial oil zone volume.
    pvt : classPVT object
        Fluid PVT properties.

    Attributes
    -------
    Vp : float
        Returns reservoir pore volume.

    """
    def __init__(self, N, p0, Pb0, cr, Swcon, m, pvt):
        self.N = N
        self.p0 = p0
        self.Pb0 = Pb0
        self.cr = 0.
        self.Swcon = 0.
        self.m = 0.
        self.pvt = pvt
        # Pore volume of HC zone
        self.Vp0 = (1.+m)*N*pvt.Bo(p0, Pb0)/(1.-Swcon)

    def Vp(self, p):
        """Returns reservoir pore volume.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).

        Returns
        -------
        Vp : float
            Reservoir pore volume (m3).

        """
        return self.Vp0*(1.-self.cr*(self.p0-p))

class classHist:
    """Class for loading production, injection and pressure history.

    Parameters
    ----------
    hist : string
       Excel xlsx file with production, injection and pressure history.

    Attributes
    -------
    Vp : float
        Returns reservoir pore volume.

    """
    def __init__(self, data):
        self.data = pd.read_excel(data, index_col='Time')
        self.data[['Np', 'Gp', 'Wp', 'Wi', 'Gi']] = \
        self.data[['Np', 'Gp', 'Wp', 'Wi', 'Gi']].interpolate(method='index')
        self.time = self.data.index.values
        self.p = self.data['p']
        self.Np = self.data['Np']
        self.Gp = self.data['Gp']
        self.Wp = self.data['Wp']
        self.Wi = self.data['Wi']
        self.Gi = self.data['Gi']

class classMBE:
    """Class for defining Material Balance Equation (MBE).

    Parameters
    ----------
    hist : object of classHist
        Production, injection and pressure history.
    pvt : object of classPVT
        Fluid PVT properties.
    res : object of classRes
        Reservoir properties.
    aquif : object of classAquifer
        Aquifer properties.

    Attributes
    -------
    xx : xxxx
        xxxx.

    """
    def __init__(self, hist, pvt, res, aquif):
        self.hist = hist
        self.pvt = pvt
        self.res = res
        self.aquif = aquif
        self.out = self.hist.data

    def Eo(self, p, Pb):
        pvt = self.pvt
        res = self.res
        return pvt.Bo(p, Pb)-pvt.Bo(res.p0, res.Pb0) \
                + (pvt.Rs(res.p0, res.Pb0)\
                - pvt.Rs(p, Pb))*pvt.Bg(p)

    def Eg(self, p):
        pvt = self.pvt
        res = self.res
        return pvt.Bo(res.p0, res.Pb0)*(pvt.Bg(p)/pvt.Bg(res.p0)-1.)

    def Efw(self, p):
        pvt = self.pvt
        res = self.res
        m = res.m
        cr = res.cr
        cw = pvt.cw
        return -(1.+m)*pvt.Bo(res.p0,res.Pb0)*(cr+cw*res.Swcon)\
                *(p-res.p0)/(1.-res.Swcon)

    def Sw(self, p, Wi, Wp, We):
        res=self.res
        pvt=self.pvt
        return ((1.+res.m)*res.N*pvt.Bo(res.p0,res.Pb0)*res.Swcon \
                /(1.-res.Swcon)/pvt.Bw(res.p0)+(Wi+We-Wp))*pvt.Bw(p)/res.Vp(p)

    def So(self, p, Pb, Np):
        res=self.res
        pvt=self.pvt
        return (res.N-Np)/res.Vp(p)*pvt.Bo(p, Pb)

    def Sg(self, p, Pb, N, Np, Gp, Gi):
        res=self.res
        pvt=self.pvt
        return pvt.Bg(p)/res.Vp(p)*(N*((pvt.Rs(res.p0, res.Pb0)-pvt.Rs(p, Pb))\
                +res.m*pvt.Bo(res.p0, res.Pb0)/pvt.Bg(res.p0)) \
                +Gi-Np*(Gp/Np-pvt.Rs(p, Pb)))

    def residual(self, p, Pb, Np, Gp, Gi, Wi, Wp):
        res=self.res
        pvt=self.pvt
        We=self.aquif.calcWe(p)
        # computing residual:
        # rhs=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))+(Wi+We)*pvt.Bw(p)
        # lhs=Np*(pvt.Bo(p,pvt.Pb)+(Gp/Np-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp*pvt.Bw(p)
        E=res.N*(self.Eo(p, Pb)+res.m*self.Eg(p)+self.Efw(p))
        Fi=(Wi+We)*pvt.Bw(p)+Gi*pvt.Bg(p)
        Fp=Np*(pvt.Bo(p, Pb)+(Gp/Np-pvt.Rs(p, Pb))*pvt.Bg(p)) \
                +Wp*pvt.Bw(p)
        #return (E-Fp+Fi)/(abs(Fp)+abs(Fi)+abs(E))
        return (E-Fp+Fi)/res.Vp0


    def calcP(self):
        res = self.res
        pvt = self.pvt

        # Create dataframe to be returned
        r = pd.DataFrame(columns=['pCalc', 'Pb', 'So', 'Sg', 'Sw', 'Vp', \
               'Bo', 'Bg', 'Bw', 'Rs', 'Eo', 'Eg', 'Efw', 'We', 'iter','resd'])
        r = pd.concat([self.hist.data, r], axis=1)
        # Solve MBE for pressure
        for index, row in r.iterrows():
            Pb = self.updatePb(row['Np'], row['Gp'], row['Gi'])
            pbrentq,brentq_r=sp.brentq(self.residual, min(pvt.pe[:]), \
                    max(pvt.pe[:]), args=(Pb, row['Np'], row['Gp'], \
                    row['Gi'], row['Wi'], row['Wp']), xtol=1e-5, rtol=1e-8, \
                    maxiter=2000, full_output=True, disp=True)
            p = pbrentq
            We = self.aquif.calcWe(p)

            # Save output data
            r.loc[index, 'pCalc'] = pbrentq
            r.loc[index, 'Pb'] = Pb
            r.loc[index, 'resd'] = self.residual(p, Pb, row['Np'],
                    row['Gp'], row['Gi'], row['Wi'], row['Wp'])
            r.loc[index, 'iter'] = brentq_r.iterations
            r.loc[index, 'Vp'] = res.Vp(p)
            r.loc[index, 'We'] = self.aquif.calcWe(p)
            r.loc[index, 'Sw'] = self.Sw(p, row['Wi'], row['Wp'], We)
            r.loc[index, 'So'] = self.So(p, Pb, row['Np'])
            r.loc[index, 'Sg'] = self.Sg(p, Pb, res.N, row['Np'], \
                    row['Gp'], row['Gi'])
            r.loc[index, 'Rs'] = pvt.Rs(p, Pb)
            r.loc[index, 'Bo'] = pvt.Bo(p, Pb)
            r.loc[index, 'Eo'] = self.Eo(p, Pb)
            r.loc[index, 'Eg'] = self.Eg(p)
        return r

    def EqPb(self,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        pmax=max(pvt.pe[:])
        return (res.N*(pvt.Rs(res.p0,res.Pb0) \
                +res.m*pvt.Bo(res.p0,res.Pb0)/pvt.Bg(res.p0))-Gp+Gi)/(res.N-Np)

    def updatePb(self,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        RsPb=self.EqPb(Np,Gp,Gi)
        r=np.interp(RsPb,pvt.Rse[::-1],pvt.pe[::-1])
#        r=max(r,pvt.pe[-2])
        return r

#------------------------------------------------------------------------
# Miscelaneous functions
#------------------------------------------------------------------------


    def calcGp(self):
        res=self.res
        teste=self.res.teste
        Np=teste.Np
        Wi=teste.Wi
        Wp=teste.Wp
        pvt=self.pvt
        dt=self.dt
        out=zeros([len(teste.p),2])
        Gp=0.
        GpOld=0.
        for i in range(len(teste.p)):
            p=teste.p[i]
            t=np.interp(Np[i],self.Np,self.dt)
            We=self.aquifero.We(p,t,Np[i])
            F=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))+(Wi[i]+We)*pvt.Bw(p)
            Gp=pvt.Rs(p,pvt.Pb)*Np[i]-pvt.Bo(p,pvt.Pb)/pvt.Bg(p)*Np[i]-pvt.Bw(p)/pvt.Bg(p)*Wp[i]+F/pvt.Bg(p)
#            Gp=pvt.Rs(p,pvt.Pb)*Np[i]+(pvt.Rs(res.p0,res.Pb0)-pvt.Rs(p,pvt.Pb))*self.N
            Gp=max(Gp,GpOld)
            pvt.Pb=self.updatePb(Np[i],Gp)
            out[i,0]=Np[i]
            out[i,1]=Gp
            GpOld=Gp
#               Reset original reservoir Bubble Pressure
        pvt.Pb=res.Pb0
        return out

    def fw(self):
        p=self.out.p
        Sw=self.out.Sw
        pvt=self.pvt
        petroPhys=self.petroPhys
        M=(petroPhys.krw(Sw)/pvt.viscw(p))/(petroPhys.kro(Sw)/pvt.visco(p))
        return 1./(1.+1./M)

    def RGO(self,Sw,So,p):
        krg=-0.6*(Sw+So-0.5)+0.3
        Bo=self.pvt.Bo(p)
        Bg=self.pvt.Bg(p)
        muo=1000.
        mug=1.
        RGO=krg*muo/mug*Bo/Bg
        return RGO

    def We_residual(self,We,p,Np,Gp,Gi,Wi,Wp):
        res=self.res
        pvt=self.pvt
        # computing residual:
        # rhs=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))+(Wi+We)*pvt.Bw(p)
        # lhs=Np*(pvt.Bo(p,pvt.Pb)+(Gp/Np-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp*pvt.Bw(p)
        E=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))
        Fi=(Wi+We)*pvt.Bw(p)+Gi*pvt.Bg(p)
        Fp=Np*(pvt.Bo(p,pvt.Pb)+(Gp/Np-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp*pvt.Bw(p)
        #return (E-Fp+Fi)/(abs(Fp)+abs(Fi)+abs(E))
        return (E-Fp+Fi)/res.Vp0

    def calcWe(self):
        out=self.out
        out.reset()
        Np=self.Np
        Gp=self.Gp
        Gi=self.Gi
        Wi=self.Wi
        Wp=self.Wp
        We=self.We
        res=self.res
        pvt=self.pvt
        dt=self.dt
        p=self.p
        Wecalc=0.

        RGO=gradient(Gp)/gradient(Np)

##        if 'out' in locals():
##            del out
##            out=[]
        for i in range(len(Np)):
            pvt.Pb=self.updatePb(Np[i],Gp[i],Gi[i])
            # print 'Residual = ', self.We_residual(0.,p[i],Np[i],Gp[i],Gi[i],Wi[i],Wp[i])
#            We[i],brentq_r=sp.bisect(self.We_residual, 0. ,10.*(Wp[i]+Np[i]), args=(p[i],Np[i],Gp[i],Gi[i],Wi[i],Wp[i]),xtol=1e-8, rtol=1e-8, maxiter=2000, full_output=True, disp=True)
            We[i],brentq_r=sp.brentq(self.We_residual, 0. ,5.*(Wp[i]+Np[i]), args=(p[i],Np[i],Gp[i],Gi[i],Wi[i],Wp[i]),xtol=1e-8, rtol=1e-8, maxiter=2000, full_output=True, disp=True)
            fval = self.We_residual(We[i],p[i],Np[i],Gp[i],Gi[i],Wi[i],Wp[i])

            # Save output data
            out.p.append(p[i])
            out.Vp.append(res.Vp(p[i]))
            out.Sw.append(self.Sw(p[i],Wi[i],Wp[i],We[i]))
            out.So.append(self.So(p[i],Np[i]))
            out.Sg.append(self.Sg(p[i],self.N,Np[i],Gp[i],Gi[i]))
            out.residual.append(fval)
            out.Eo.append(self.Eo(p[i]))
            out.Eg.append(self.Eg(p[i]))
            out.Efw.append(self.Efw(p[i]))
            out.We.append(self.We[i])
            #out.fw.append(gradient(Wp)[i]/(gradient(Np)[i]+gradient(Wp)[i]))

            # Compute new bubble pressure
            out.Pb.append(pvt.Pb)
#             print '[iter = ',i,']'
#             print 'brentq iter, fcalls = ',brentq_r.iterations,brentq_r.function_calls
#             print 'p [kgf/cm2], residual = ',p[i],fval
#             print 'We [MM m3] = ', We[i]/1e6
#             print 'MBE oil = ',out.So[-1]*out.Vp[-1]/pvt.Bo(p[i],pvt.Pb)/(self.N-Np[i])
#             print 'Pb [kgf/cm2] = ',out.Pb[-1]
#             print 'RGO = ', RGO[i]
#             print 'FR Oil, FR Gas = ',Np[i]/self.N, Gp[i]/(self.N*pvt.Rs(res.p0,res.Pb0))
#             print 'So, Sg, Sw = ',out.So[-1],out.Sg[-1],out.Sw[-1]
#             print '--------------------------------------------'

        return out
