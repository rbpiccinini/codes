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

    def __init__(self, pvt_table, co, T):
        # reading pvt table
        self.pe = pvt['P']
        self.Boe = pvt['BO']
        self.Rse = pvt['RS']
        self.Bge = pvt['BG']
        self.Bwe = pvt['BW']

        # reading other properties
        self.co = co
        self.T = T

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
            return interp(Pb, self.pe[::-1], self.Boe[::-1]) \
                    * (1.-self.co*(p-Pb))
        else:
            return interp(p, self.pe[::-1], self.Boe[::-1])

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
        if type(self.Ze)==type(1.0):
            Bgi=1.033/self.pe*(self.T)/293.
            self.Ze=self.Bge/Bgi

        Z=interp(p,self.pe[::-1],self.Ze[::-1])
        return 1.033/p*(self.T)/293.*Z

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
        return interp(p,self.pe[::-1],self.Bwe[::-1])

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
            return interp(Pb,self.pe[::-1],self.Rse[::-1])
        else:
            return interp(p,self.pe[::-1],self.Rse[::-1])



class classAquifer:
    """Class for aquifer modeling.

    Parameters
    ----------
    WW : float
        Aquifer volume (m3).
    p0 : float
        Initial aquifer pressure (kgf/cm2).
    cw : float
        water compressbility (cm2/kgf).
    cr : float
        pore volume compressibility (cm2/kgf).

    Attributes
    -------
    calcWe : float
        Returns aquifer influx for a given final pressure.

    """
    def __init__(self, WW, p0, cw, cr):
        self.WW = WW
        self.p0 = p0
        self.cw = cw
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
        ct=self.cr+self.cw
        We=ct*self.WW*(self.p0-p)
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
        Ratio of initial gas-cap volume to initial oil zone volume

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
    def __init__(self, file):
        self.data = pd.read_excel(file)
        self.data =


class classMBE:
"""Class for defining material balance equation.

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
        Ratio of initial gas-cap volume to initial oil zone volume

    Attributes
    -------
    Vp : float
        Returns reservoir pore volume.

    """

    p=0.
    pold=0.

    N=0.
    Np=0.
    Gp=0.

    Wp=0.
    Wi=0.
    We=0.

    dt=0.

    aquifero=classAquifero()

    res=classRes()
    pvt=classPVT()
    out=classHist()

    petroPhys=classPetrophys()


    def Eo(self,p):
        pvt=self.pvt
        res=self.res
        return pvt.Bo(p,pvt.Pb)-pvt.Bo(res.p0,res.Pb0)+(pvt.Rs(res.p0,res.Pb0)-pvt.Rs(p,pvt.Pb))*pvt.Bg(p)

    def Eg(self,p):
        pvt=self.pvt
        res=self.res
        return pvt.Bo(res.p0,res.Pb0)*(pvt.Bg(p)/pvt.Bg(res.p0)-1.)

    def Efw(self,p):
        pvt=self.pvt
        res=self.res
        m=res.m
        cr=res.cr
        cw=pvt.cw(p,res.p0)
        return -(1.+m)*pvt.Bo(res.p0,res.Pb0)*(cr+cw*res.Sw0)*(p-res.p0)/(1.-res.Sw0)

    def Sw(self,p,Wi,Wp,We):
        res=self.res
        pvt=self.pvt
        return ((1.+res.m)*self.N*pvt.Bo(res.p0,res.Pb0)*res.Sw0/(1.-res.Sw0)/pvt.Bw(res.p0)+(Wi+We-Wp))*pvt.Bw(p)/res.Vp(p)

    def So(self,p,Np):
        res=self.res
        pvt=self.pvt
        return (self.N-Np)/res.Vp(p)*pvt.Bo(p,pvt.Pb)

    def Sg(self,p,N,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        return pvt.Bg(p)/res.Vp(p)*(N*((pvt.Rs(res.p0,res.Pb0)-pvt.Rs(p,pvt.Pb))+res.m*pvt.Bo(res.p0,res.Pb0)/pvt.Bg(res.p0))+Gi-Np*(Gp/Np-pvt.Rs(p,pvt.Pb)))

    def residual(self,p,Np,Gp,Gi,Wi,Wp,dt):
        res=self.res
        pvt=self.pvt
        We=self.aquifero.We(p,dt,Np)
        # computing residual:
        # rhs=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))+(Wi+We)*pvt.Bw(p)
        # lhs=Np*(pvt.Bo(p,pvt.Pb)+(Gp/Np-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp*pvt.Bw(p)
        E=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))
        Fi=(Wi+We)*pvt.Bw(p)+Gi*pvt.Bg(p)
        Fp=Np*(pvt.Bo(p,pvt.Pb)+(Gp/Np-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp*pvt.Bw(p)
        #return (E-Fp+Fi)/(abs(Fp)+abs(Fi)+abs(E))
        return (E-Fp+Fi)/res.Vp0

    def VOIP(self):
        p=self.res.teste.p
        Np=self.res.teste.Np
        Wp=self.res.teste.Wp
        Wi=self.res.teste.Wi
        We=zeros(len(Wi))
        Gi=self.res.teste.Gi
        Gp=self.res.teste.Gp
        res=self.res
        pvt=self.pvt
        N=[]
        for i in range(len(p)):
            E=self.Eo(p[i])+res.m*self.Eg(p[i])+self.Efw(p[i])
            Fi=(Wi[i]+We[i])*pvt.Bw(p[i])+Gi[i]*pvt.Bg(p[i])
            Fp=Np[i]*(pvt.Bo(p[i],pvt.Pb)+(Gp[i]/Np[i]-pvt.Rs(p[i],pvt.Pb))*pvt.Bg(p[i]))+Wp[i]*pvt.Bw(p[i])
            N.append((Fp-Fi)/E)

    def VOIPW(self):
        p=self.res.teste.p
        Np=self.res.teste.Np
        Wp=self.res.teste.Wp
        Wi=self.res.teste.Wi
        We=zeros(len(Wi))
        Gi=self.res.teste.Gi
        Gp=self.res.teste.Gp
        res=self.res
        pvt=self.pvt
        A=zeros(len(p))
        B=zeros(len(p))
        C=zeros(len(p))
        for i in range(len(p)):
            A[i]=self.Eo(p[i])+res.m*self.Eg(p[i])+self.Efw(p[i])
            B[i]=(self.res.cr+self.pvt.cw)*pvt.Bw(p[i])*(self.res.p0-p[i])
            C[i]=Np[i]*(pvt.Bo(p[i],pvt.Pb)+(Gp[i]/Np[i]-pvt.Rs(p[i],pvt.Pb))*pvt.Bg(p[i]))+Wp[i]*pvt.Bw(p[i])-Wi[i]*pvt.Bw(p[i])-Gi[i]*pvt.Bg(p[i])
        # least square regression of Ai*N+Bi*We-Ci=0

    def runHist(self,pmin,pmax):
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

        px=linspace(pmin, pmax, 100)
        out.px.append(px/res.Pb0)
        RGO=gradient(Gp)/gradient(Np)

##        if 'out' in locals():
##            del out
##            out=[]
        for i in range(len(Np)):
            pvt.Pb=self.updatePb(Np[i],Gp[i],Gi[i])
#            pbisect,bisect_r=sp.bisect(self.residual, min(pvt.pe[:]) ,0.99*res.p0, args=(Np[i],Gp[i],Wi[i],Wp[i],dt[i]),xtol=1e-5, rtol=1e-6, maxiter=1000, full_output=True, disp=True)
            pbrentq,brentq_r=sp.brentq(self.residual, min(pvt.pe[:]) ,2.*res.p0, args=(Np[i],Gp[i],Gi[i],Wi[i],Wp[i],dt[i]),xtol=1e-5, rtol=1e-8, maxiter=2000, full_output=True, disp=True)
#            pnewton=sp.newton(self.residual, pbisect, args=(Np[i],Gp[i],Wi[i],Wp[i],dt[i]),tol=1e-5, maxiter=1000)

            p=pbrentq
            pnewton=pbrentq
            fval = self.residual(p,Np[i],Gp[i],Gi[i],Wi[i],Wp[i],dt[i])


            if mod(i,2)==0:
                resd=zeros(len(px))
                figure(1)
                for j in range(len(px)):
                    resd[j]=self.residual(px[j],Np[i],Gp[i],Gi[i],Wi[i],Wp[i],We[i])
                plot(px/res.Pb0,resd,'k-',linewidth=0.5)
                plot(pnewton/res.Pb0,self.residual(pnewton,Np[i],Gp[i],Gi[i],Wi[i],Wp[i],We[i]),'sk',markersize=5)
                plot(pvt.Pb/res.Pb0,self.residual(pvt.Pb,Np[i],Gp[i],Gi[i],Wi[i],Wp[i],We[i]),'ok',markersize=5,markerfacecolor='#ffffff',markeredgewidth=2)
                xlabel('Pressure [units of $P_{b1}$]')
                ylabel('Residual Function')
                grid(1)
                ylim([-0.1,0.15])
                xlim([0,2.2])


            # Save output data
            out.px.append(resd)
            out.Eo_p1.append(pvt.Bo(p,pvt.Pb)-pvt.Bo(res.p0,res.Pb0))
            out.Eo_p2.append((pvt.Rs(res.p0,res.Pb0)-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))
            out.p.append(p)
            out.Vp.append(res.Vp(p))
            out.Sw.append(self.Sw(p,Wi[i],Wp[i],self.aquifero.We(p,dt[i],Np[i])))
            out.So.append(self.So(p,Np[i]))
            out.Sg.append(self.Sg(p,self.N,Np[i],Gp[i],Gi[i]))
            out.Rs.append(pvt.Rs(p,pvt.Pb))
            out.Bo.append(pvt.Bo(p,pvt.Pb))
            out.residual.append(fval)
            out.Eo.append(self.Eo(p))
            out.Eg.append(self.Eg(p))
            out.Efw.append(self.Efw(p))
            out.We.append(self.aquifero.We(p,dt[i],Np[i]))
            out.fw.append(gradient(Wp)[i]/(gradient(Np)[i]+gradient(Wp)[i]))
            out.Pa.append(self.aquifero.pa(out.We[-1]))
            out.Fp.append(Np[i]*(pvt.Bo(p,pvt.Pb)+(Gp[i]/Np[i]-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp[i]*pvt.Bw(p))
            out.iter.append(brentq_r.iterations)

            # Compute new bubble pressure
            out.Pb.append(pvt.Pb)
            print '[iter = ',i,']'
            print 'brentq iter, fcalls = ',brentq_r.iterations,brentq_r.function_calls
            print 'p [kgf/cm2], residual = ',p,fval
            print 'MBE oil = ',out.So[-1]*out.Vp[-1]/pvt.Bo(p,pvt.Pb)/(self.N-Np[i])
            print 'Pb [kgf/cm2] = ',out.Pb[-1]
            print 'RGO = ', RGO[i]
            print 'FR Oil, FR Gas = ',Np[i]/self.N, Gp[i]/(self.N*pvt.Rs(res.p0,res.Pb0))
            print 'So, Sg, Sw = ',out.So[-1],out.Sg[-1],out.Sw[-1]
            print '--------------------------------------------'

        plot(pnewton/res.Pb0,self.residual(pnewton,Np[i],Gp[i],Gi[i],Wi[i],Wp[i],We[i]),'sk',markersize=5,label='solution')
        plot(pvt.Pb/res.Pb0,self.residual(pvt.Pb,Np[i],Gp[i],Gi[i],Wi[i],Wp[i],We[i]),'ok',markersize=5,markerfacecolor='#ffffff',markeredgewidth=2,label='bubble-point pressure')
        legend(loc='upper right')
        subplots_adjust(bottom=0.15,left=0.15)
        return out

    def EqPb(self,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        pmax=max(pvt.pe[:])
        return (self.N*(pvt.Rs(res.p0,res.Pb0)+res.m*pvt.Bo(res.p0,res.Pb0)/pvt.Bg(res.p0))-Gp+Gi)/(self.N-Np)

    def updatePb(self,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        RsPb=self.EqPb(Np,Gp,Gi)
        r=interp(RsPb,pvt.Rse[::-1],pvt.pe[::-1])
        r=max(r,pvt.pe[-2])
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
            t=interp(Np[i],self.Np,self.dt)
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
            print 'Residual = ', self.We_residual(0.,p[i],Np[i],Gp[i],Gi[i],Wi[i],Wp[i])
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
            print '[iter = ',i,']'
            print 'brentq iter, fcalls = ',brentq_r.iterations,brentq_r.function_calls
            print 'p [kgf/cm2], residual = ',p[i],fval
            print 'We [MM m3] = ', We[i]/1e6
            print 'MBE oil = ',out.So[-1]*out.Vp[-1]/pvt.Bo(p[i],pvt.Pb)/(self.N-Np[i])
            print 'Pb [kgf/cm2] = ',out.Pb[-1]
            print 'RGO = ', RGO[i]
            print 'FR Oil, FR Gas = ',Np[i]/self.N, Gp[i]/(self.N*pvt.Rs(res.p0,res.Pb0))
            print 'So, Sg, Sw = ',out.So[-1],out.Sg[-1],out.Sw[-1]
            print '--------------------------------------------'

        return out
