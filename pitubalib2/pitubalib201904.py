import numpy as np
import pandas as pd
import scipy.optimize as sp
#from timeit import default_timer as timer


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
    Bwi : float
        Bw at reference pressure (sm3/m3)
    pw : float
        Reference pressure for Bwi (kgf/cm2)
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

    def __init__(self, pvt_table, co, cw, Bwi, pw, T):
        # reading pvt table
        self.pvt_table = pd.read_excel(pvt_table).sort_values(by='p',
                                      ascending=True)
        self.pe = self.pvt_table['p']
        self.Boe = self.pvt_table['Bo']
        self.Rse = self.pvt_table['Rs']
        self.Bge = self.pvt_table['Bg']
        self.Bwe = self.pvt_table['Bw']
        # reading other properties
        self.co = co
        self.cw = cw
        self.Bwi = Bwi
        self.pw = pw
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

        return np.where(p >= Pb,
                         np.interp(Pb, self.pe, self.Boe)*(1.-self.co*(p-Pb)),
                         np.interp(p, self.pe, self.Boe))

    def dBo(self, p, Pb):
        """Returns dBo/dp for given bubble-point and oil phase pressure.

        Parameters
        ----------
        p : float
            Oil phase pressure (kgf/cm2).
        Pb : float
            Oil phase bubble-point pressure (kgf/cm2).

        Returns
        -------
        dBo/dp : float
            Derivative of formation-volume factor of oil phase (cm^2/kgf).

        """
        return np.where(p >= Pb,
                        -self.co*self.Bo(p, Pb),
                        np.interp(p,
                                  self.pe,
                                  np.gradient(self.Boe)/np.gradient(self.pe)))


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
        Z=np.interp(p, self.pe, self.Ze)
        return 1.033/p*(self.T)/288.*Z

    def dBg(self, p):
        """Returns dBg/dp for given gas phase pressure.

        Parameters
        ----------
        p : float
            Gas phase pressure (kgf/cm2).

        Returns
        -------
        dBg/dp : float
            Derivative of formation-volume factor of gas phase (cm^2/kgf).
        """
        Z=np.interp(p, self.pe, self.Ze)
        dBg = -1.0/self.Bge**2*np.gradient(self.Bge)/np.gradient(self.pe)
        return np.interp(p, self.pe, dBg)

    def Bw(self, p):
        """Returns Bw for given water phase pressure.

        Parameters
        ----------
        p : float
            Water phase pressure (kgf/cm2).

        Returns
        -------
        Bw : float
            Formation-volume factor of water phase (m3/sm3).

        """
        return self.Bwi*(1.0-self.cw*(p-self.pw))

    def dBw(self, p):
        """Returns dBw/dp for given gas phase pressure.

        Parameters
        ----------
        p : float
            Water phase pressure (kgf/cm2).

        Returns
        -------
        dBw/dp : float
            Derivative of formation-volume factor of water phase (cm^2/kgf).

        """
        return -self.Bw(p)*self.cw

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
        dRs/dp: float
            Gas-oil solubility ratio (m3/m3).

        """
        return np.where(p >= Pb,
                        np.interp(Pb,self.pe,self.Rse),
                        np.interp(p,self.pe,self.Rse))

    def dRs(self,p,Pb):
        """Returns dRs/dp for given oil-phase pressure.

        Parameters
        ----------
        p : float
            Water phase pressure (kgf/cm2).
        Pb : float
            Oil phase bubble-point pressure (kgf/cm2).

        Returns
        -------
        dRs/dp: float
            Derivative of gas-oil solubility ratio (cm^2/kgf).

        """
        return np.where(p >= Pb,
                        0.,
                        np.interp(p,
                                  self.pe,
                                  np.gradient(self.Rse)/np.gradient(self.pe)))




class classAquifer:
    """Class for aquifer properties and modeling.

    Parameters
    ----------
    Waq : float
        Aquifer volume (m3).
    p0 : float
        Initial aquifer pressure (kgf/cm^2).
    cr : float
        pore volume compressibility (cm^2/kgf).
    pvt : classPVT object
        Fluid PVT properties.
    J : float
        Aquifer productivity (m^3/d)/(kgf/cm^2).

    Attributes
    -------
    dWe : float
        Increment of aquifer influx.

    """
    def __init__(self, Waq, paq0, cr, pvt, J):
        self.pvt = pvt
        self.Waq = Waq
        self.paq0 = paq0
        self.cr = cr
        self.J = J

    def dWe(self, paq, p, dt):
        """Returns the increment of aquifer influx.

        Parameters
        ----------
        dp : float
            Difference between final and initial aquifer pressure inside a
            time-step (kgf/cm2).
        dt : float
            Time interval (days).

        Returns
        -------
        dWe : float
            Increment of aquifer influx (m3).

        """
        ct=self.cr+self.pvt.cw
        dWe = ct*self.Waq*(paq-p)*(1.0-np.exp(-self.J/self.Waq/ct*dt))
        return dWe

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
        Returns reservoir pore volume (m3).

    """
    def __init__(self, N, p0, Pb0, cr, Swcon, m, pvt):
        self.N = N
        self.p0 = p0
        self.Pb0 = Pb0
        self.cr =cr
        self.Swcon = Swcon
        self.m = m
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
    xlsx : string
       Excel xlsx file with production, injection and pressure history.
    """
    def __init__(self, xlsx):
        self.data = pd.read_excel(xlsx, index_col='Dates')
        self.data[['Np', 'Gp', 'Wp', 'Wi', 'Gi']] = \
        self.data[['Np', 'Gp', 'Wp', 'Wi', 'Gi']].interpolate(method='index')
        self.data = self.data.groupby(by=self.data.index).mean()
        self.date = self.data.index.values
        self.time = (self.data.index - self.data.index[0]).days.values
        self.data['time'] = self.time
        self.data['dtime'] = np.gradient(self.time)
        self.p = self.data['p']
        self.Np = self.data['Np']
        self.Gp = self.data['Gp']
        self.Wi = self.data['Wi']
        self.Gi = self.data['Gi']

class classOutput:
    """Class for saving output data.

    Parameters
    ----------
    hist : dataframe
       Excel xlsx file with production, injection and pressure history.

    Attributes
    -------
    append : None
        Appends data to output dataframe.

    """
    def __init__(self, hist):
        # create column levels
        hist = hist.data
        self.pvt_rock = pd.DataFrame(index=hist.index,
                                columns=['Bo1', 'Bg1', 'Bw1', 'Rs1', 'Pb1',
                                         'Bo2', 'Bg2', 'Bw2', 'Rs2', 'Pb2',
                                         'co', 'cw', 'cr'])
        self.ooip = pd.DataFrame(index=hist.index, columns=['m', 'N'])
        self.aquifer = pd.DataFrame(index=hist.index,
                                    columns=['Waq', 'We', 'J',
                                             'paq1', 'paq2',])
        self.reservoir = pd.DataFrame(index=hist.index,
                                columns=['p1', 'p2', 'So', 'Sg', 'Sw'])
        self.drive = pd.DataFrame(index=hist.index,
                             columns=['oil and solution gas expansion',
                                      'gas cap expansion',
                                      'pore volume compaction',
                                      'aquifer influx and water injection'])

        # concatenate all dfs in a single multilevel df
        keys = ['reservoir history',
                'fluid and rock properties',
                'original HC volumes',
                'aquifer data',
                'reservoir data',
                'drive mechanisms']
        self.data = pd.concat([hist,
                               self.pvt_rock,
                               self.ooip,
                               self.aquifer,
                               self.reservoir,
                               self.drive], axis=1)
    def append(self, index,
               Bo1, Bg1, Bw1, Rs1, Pb1, Bo2, Bg2, Bw2, Rs2, Pb2, co, cw, cr,
               m, N,
               Waq, We, J, paq1, paq2,
               p1, p2, So, Sg, Sw,
               drive_og, drive_gascap, drive_pv, drive_wiwe):

        self.data.loc[index, self.pvt_rock.columns] = \
                      [Bo1, Bg1, Bw1, Rs1, Pb1, Bo2, Bg2, Bw2, Rs2, Pb2,
                       co, cw, cr]
        self.data.loc[index, self.ooip.columns] = [m, N]
        self.data.loc[index, self.aquifer.columns] = [Waq, We, J, paq1, paq2]
        self.data.loc[index, self.reservoir.columns] = [p1, p2, So, Sg, Sw]
        self.data.loc[index, self.drive.columns] = \
                                 [drive_og, drive_gascap, drive_pv, drive_wiwe]



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
    output : object of ClassOutput
        Template for output data.

    Attributes
    -------
    xx : xxxx
        xxxx.

    """
    def __init__(self, hist, pvt, res, aquif, output):
        self.hist = hist
        self.pvt = pvt
        self.res = res
        self.aquif = aquif
        self.output = output

    def Eo(self, p, Pb):
        """Returns the expansion of oil and solution gas for final bubble-point
        and reservoir pressures.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Pb : float
            Oil-phase bubble-point pressure (kgf/cm2).

        Returns
        -------
        Eo : float
            Expansion of oil and solution gas (rm3/sm3).
        """
        pvt = self.pvt
        res = self.res
        return pvt.Bo(p, Pb)-pvt.Bo(res.p0, res.Pb0) \
                + (pvt.Rs(res.p0, res.Pb0)\
                - pvt.Rs(p, Pb))*pvt.Bg(p)

    def Eg(self, p):
        """Returns the expansion of free gas for a final reservoir pressure.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).

        Returns
        -------
        Eg : float
            Expansion of reservoir free gas (rm3/sm3).
        """
        pvt = self.pvt
        res = self.res
        return pvt.Bo(res.p0, res.Pb0)*(pvt.Bg(p)/pvt.Bg(res.p0)-1.)

    def Efw(self, p):
        """Returns the expansion of connate water and the compaction of pore
        volume for a final reservoir pressure.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).

        Returns
        -------
        Efw : float
            Expansion of connate water and compaction of pore
            volume for a final reservoir pressure (rm3/sm3).
        """
        pvt = self.pvt
        res = self.res
        m = res.m
        cr = res.cr
        cw = pvt.cw
        return -(1.+m)*pvt.Bo(res.p0,res.Pb0)*(cr+cw*res.Swcon)\
                *(p-res.p0)/(1.-res.Swcon)

    def Sw(self, p, Wi, Wp, We):
        """Returns final water saturation.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Wi : float
            Cumulative injection of water (sm3)
        Wp : float
            Cumulative production of water (sm3)
        We : float
            Cumnulative aquifer influx (sm3)

        Returns
        -------
        Sw : float
            Final water-phase saturation.
        """
        res=self.res
        pvt=self.pvt
        return ((1.+res.m)*res.N*pvt.Bo(res.p0,res.Pb0)*res.Swcon \
                /(1.-res.Swcon)/pvt.Bw(res.p0)+(Wi+We-Wp))*pvt.Bw(p)/res.Vp(p)

    def So(self, p, Pb, Np):
        """Returns final oil saturation.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Pb : float
            Oil-phase bubble-point pressure (kgf/cm2).
        Np : float
            Cumulative production of oil (sm3)

        Returns
        -------
        So : float
            Final oil-phase saturation.
        """
        res=self.res
        pvt=self.pvt
        return (res.N-Np)/res.Vp(p)*pvt.Bo(p, Pb)

    def Sg(self, p, Pb, N, Np, Gp, Gi):
        """Returns final gas saturation.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Pb : float
            Oil-phase bubble-point pressure (kgf/cm2).
        N : float
            Original volume of oil in place (sm3)
        Np : float
            Cumulative production of oil (sm3)
        Gp : float
            Cumulative production of gas (sm3)
        Gi : float
            Cumulative injection of gas (sm3)

        Returns
        -------
        Sg : float
            Final gas-phase saturation.
        """
        res=self.res
        pvt=self.pvt
        return pvt.Bg(p)/res.Vp(p)*(N*((pvt.Rs(res.p0, res.Pb0)-pvt.Rs(p, Pb))\
                +res.m*pvt.Bo(res.p0, res.Pb0)/pvt.Bg(res.p0)) \
                +Gi-Gp+Np*pvt.Rs(p, Pb))

    def Fp(self, p, Pb, Np, Gp, Wp):
        """Returns the produced volume in reservoir conditions.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Pb : float
            Oil-phase bubble-point pressure (kgf/cm2).
        N : float
            Original volume of oil in place (sm3).
        Np : float
            Cumulative production of oil (sm3).
        Gp : float
            Cumulative production of gas (sm3).
        Wp : float
            Cumulative production of water (sm3).

        Returns
        -------
        Fp : float
            Produced volume of oi, gas and water in reservoir conditions.
        """
        res=self.res
        pvt=self.pvt
        return Np*pvt.Bo(p, Pb)+(Gp-Np*pvt.Rs(p, Pb))*pvt.Bg(p)+Wp*pvt.Bw(p)

    def Fpog(self, p, Pb, Np, Gp):
        """Returns the produced volume of oil and gas in reservoir conditions.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Pb : float
            Oil-phase bubble-point pressure (kgf/cm2).
        N : float
            Original volume of oil in place (sm3).
        Np : float
            Cumulative production of oil (sm3).
        Gp : float
            Cumulative production of gas (sm3).
        Wp : float
            Cumulative production of water (sm3).

        Returns
        -------
        Fp : float
            Produced volume of oi, gas and water in reservoir conditions.
        """
        res=self.res
        pvt=self.pvt
        return Np*pvt.Bo(p, Pb)+(Gp-Np*pvt.Rs(p, Pb))*pvt.Bg(p)

    def Fi(self, p, Wi, We, Gi):
        """Returns the injected volume of water and gas in reservoir conditions.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Wi : float
            Cumulative injection of water (sm3)
        We : float
            Cumulative influx of water (sm3)
        Gi : float
            Cumulative injection of gas (sm3)

        Returns
        -------
        Fi : float
            Injected volume of water and gas in reservoir conditions.
        """
        pvt=self.pvt
        return pvt.Bw(p)*(We+Wi)+Gi*pvt.Bg(p)

    def Fiw(self, p, Wi, We, Wp):
        """Returns the net injected volume of water in reservoir conditions.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Wi : float
            Cumulative injection of water (sm3)
        We : float
            Cumulative influx of water (sm3)
        Wp : float
            Cumulative production of water (sm3)

        Returns
        -------
        Fiw : float
            Injected volume of water and gas in reservoir conditions.
        """
        pvt=self.pvt
        return pvt.Bw(p)*(We+Wi-Wp)

    def Fig(self, p, Gi):
        """Returns the injected volume of gas in reservoir conditions.

        Parameters
        ----------
        p : float
            Reservoir pressure (kgf/cm2).
        Gi : float
            Cumulative injection of gas (sm3)

        Returns
        -------
        Fig : float
            Injected volume of gas in reservoir conditions.
        """
        pvt=self.pvt
        return pvt.Bg(p)*Gi

    def residual(self, p, paq, dt, Pb, Np, Gp, Gi, Wi, Wp, Weold):
        res=self.res
        pvt=self.pvt
        We=Weold+self.aquif.dWe(paq, p, dt)
        # computing residual:
        # rhs=self.N*(self.Eo(p)+res.m*self.Eg(p)+self.Efw(p))+(Wi+We)*pvt.Bw(p)
        # lhs=Np*(pvt.Bo(p,pvt.Pb)+(Gp/Np-pvt.Rs(p,pvt.Pb))*pvt.Bg(p))+Wp*pvt.Bw(p)
        E=res.N*(self.Eo(p, Pb)+res.m*self.Eg(p)+self.Efw(p))
        Fiw = self.Fiw(p, Wi, We, Wp)
        Fig = self.Fig(p, Gi)
        Fpog = self.Fpog(p, Pb, Np, Gp)

        return (E+Fiw+Fig)/Fpog - 1


    def calcP(self, info=True):
        hist = self.hist
        res = self.res
        pvt = self.pvt
        aquif = self.aquif
        output = self.output


        # Initialize p, paq and We values
        pold = res.p0
        paq = aquif.paq0
        Weold = 0.
        We = 0.

        # Solve MBE for pressure
        for index, row in hist.data.iterrows():
            Pb = self.updatePb(row['Np'], row['Gp'], row['Gi'])
            dt = row['dtime']

            # avoid errors at initial time step.
            if row['Np']+row['Gp'] <1e-15:
                continue
            else:
                pbrentq,brentq_r=sp.brentq(self.residual, pvt.pe.min(), \
                        pvt.pe.max(), args=(paq, dt, Pb, row['Np'], row['Gp'], \
                        row['Gi'], row['Wi'], row['Wp'], Weold),\
                        xtol=1e-2, rtol=1e-6, \
                        maxiter=2000, full_output=True, disp=True)
            p = pbrentq

            Fpog = self.Fpog(p, Pb, row['Np'], row['Gp'])
            residual = 100*self.residual(p, paq, dt, Pb, row['Np'], row['Gp'], \
                        row['Gi'], row['Wi'], row['Wp'], Weold)

            We = Weold+aquif.dWe(paq, p, dt)
            paq = paq - aquif.dWe(paq, p, dt)/(aquif.Waq*(pvt.cw+aquif.cr))
            Eo = self.Eo(p, Pb)
            Eg = self.Eg(p)
            Efw = self.Efw(p)
            Fiw = self.Fiw(p, row['Wi'], We, row['Wp'])
            if info:
                print('{:%Y-%m-%d %H:%M} -- Matbal error (% of Fpog) = {:-5.2f}'.format(index, residual))

#            if row['Wp'] > (We + row['Wi']):
#                raise Exception('Wp > We + Wi: Produced water is greater than the \
#                                sum of aquifer influx and injected water. \
#                                Please revise input data.')

            # Save output data
            output.append(index = index,
                               Bo1 = pvt.Bo(res.p0, res.Pb0),
                               Bg1 = pvt.Bg(res.p0),
                               Bw1 = pvt.Bw(res.p0),
                               Rs1 = pvt.Rs(res.p0, res.Pb0),
                               Pb1 = res.Pb0,
                               Bo2 = pvt.Bo(p,Pb),
                               Bg2 = pvt.Bg(p),
                               Bw2 = pvt.Bw(p),
                               Rs2 = pvt.Rs(p, Pb),
                               Pb2 = Pb,
                               co = pvt.co,
                               cw = pvt.cw,
                               cr = res.cr,
                               m = res.m,
                               N = res.N,
                               Waq = aquif.Waq,
                               We = We,
                               J = aquif.J,
                               paq1 = aquif.paq0,
                               paq2 = paq,
                               p1 = res.p0,
                               p2 = p,
                               So = self.So(p, Pb, row['Np']),
                               Sg = self.Sg(p, Pb, res.N, row['Np'],
                                            row['Gp'], row['Gi']),
                               Sw = self.Sw(p, row['Wi'], row['Wp'], We),
                               drive_og = Eo*res.N/Fpog,
                               drive_gascap = res.m*Eg*res.N/Fpog,
                               drive_pv = Efw*res.N/Fpog,
                               drive_wiwe = Fiw/Fpog)
            # Save old p and We values
            pold = p
            Weold = We
        return output.data

    def calcWe(self):
        res = self.res
        pvt = self.pvt
        hist = self.hist
        output = self.output

        # Solve MBE for aquifer influx
        for index, row in hist.data.iterrows():
            Pb = self.updatePb(row['Np'], row['Gp'], row['Gi'])
            p = row['p']
            E = res.N*(self.Eo(p, Pb)+res.m*self.Eg(p)+self.Efw(p))
            Fi = (row['Wi'])*pvt.Bw(p)+row['Gi']*pvt.Bg(p)  # Fi - Bw2*We
            Fp = row['Np']*(pvt.Bo(p, Pb)+\
                    -pvt.Rs(p, Pb)*pvt.Bg(p)) +row['Gp']*pvt.Bg(p) \
                 +row['Wp']*pvt.Bw(p)
            We = (Fp-E-Fi)/pvt.Bw(p)

            Fp = self.Fp(p, Pb, row['Np'], row['Gp'])
            Eo = self.Eo(p, Pb)
            Eg = self.Eg(p)
            Efw = self.Efw(p)
            Fi = self.Fi(p, row['Wp'], row['Wi'], We)

            # Save output data
            output.append(index = index,
                               Bo1 = pvt.Bo(res.p0, res.Pb0),
                               Bg1 = pvt.Bg(res.p0),
                               Bw1 = pvt.Bw(res.p0),
                               Rs1 = pvt.Rs(res.p0, res.Pb0),
                               Pb1 = res.Pb0,
                               Bo2 = pvt.Bo(p,Pb),
                               Bg2 = pvt.Bg(p),
                               Bw2 = pvt.Bw(p),
                               Rs2 = pvt.Rs(p, Pb),
                               Pb2 = Pb,
                               co = pvt.co,
                               cw = pvt.cw,
                               cr = res.cr,
                               m = res.m,
                               N = res.N,
                               Waq = np.nan,
                               We = We,
                               J = np.nan,
                               paq1 = np.nan,
                               paq2 = np.nan,
                               p1 = res.p0,
                               p2 = p,
                               So = self.So(p, Pb, row['Np']),
                               Sg = self.Sg(p, Pb, res.N, row['Np'],
                                            row['Gp'], row['Gi']),
                               Sw = self.Sw(p, row['Wi'], row['Wp'], We),
                               drive_og = Eo*res.N/Fp,
                               drive_gascap = res.m*Eg*res.N/Fp,
                               drive_pv = Efw*res.N/Fp,
                               drive_wiwe = Fi/Fp)
        return output.data


    def EqPb(self,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        return (res.N*(pvt.Rs(res.p0,res.Pb0) \
                +res.m*pvt.Bo(res.p0,res.Pb0)/pvt.Bg(res.p0))-Gp+Gi)/(res.N-Np)

    def updatePb(self,Np,Gp,Gi):
        res=self.res
        pvt=self.pvt
        RsPb=self.EqPb(Np,Gp,Gi)
        r=np.interp(RsPb,pvt.Rse,pvt.pe)
#        r=max(r,pvt.pe[-2])
        return r
