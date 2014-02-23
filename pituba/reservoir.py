from pylab import *
import scipy.optimize as sp
from pitubalib import *

rc('font',family='serif')
rc('text', usetex=False)
rc('font', size=16)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
rc('legend', fontsize=16)
rc('figure',figsize=[8,6])

fbottom=0.15
fleft=0.15

#---------------------------------------------------------------------------
# PVT data from reservoir fluid
#---------------------------------------------------------------------------

## PVT
fluido=classPVT()
fluido.pe=array([200.,170.,140.,110.,80.,60.,48.,35.,25.,15.,1.033])
fluido.Boe=array([1.165,1.149,1.133,1.116,1.100,1.089,1.083,1.076,1.069,1.062,1.0000])
fluido.Rse=array([66.9,57.8,48.6,39.4,30.2,24.1,20.4,16.4,13.1,9.4,0.0000])
fluido.Bge=array([0.0047,0.055,0.0067,0.0085,0.0117,0.0156,0.0195,0.0310,0.0435,0.0720,1.0000])
fluido.Bwe=1.0*ones(len(fluido.pe))
fluido.viscoe=array([5.81000,7.23000,8.05000,9.05000,10.23000,11.58000,12.17000,12.80000,14.12000,14.83000,15.57000,16.33000,17.12000,17.96000,18])
fluido.viscwe=1.0*ones(len(fluido.pe))
fluido.co=142e-6
fluido.Pb=60.
fluido.T=273.+57.

#---------------------------------------------------------------------------
# COMPUTING EQUATION COEFFICIENTS
# Np*(Bo2+(Rp-Rso2)*Bg2)+Wp*Bw2=N*(Eo+m*Eg+Efw)+(Wi+We)*Bw2+Gi*Bg2
#---------------------------------------------------------------------------
Eq=classEBM()
prod=genfromtxt('production.dat',delimiter='\t',skiprows=1,usecols = (0,1,2,3,4), converters = {4: lambda s: float(s or 0)})
Eq.Np=prod[:,0]
Eq.Gp=prod[:,1]
Eq.Wp=prod[:,2]
Eq.Wi=prod[:,3]
Eq.We=prod[:,4]
Eq.Gi=zeros(len(Eq.Np))
Eq.dt=zeros(len(Eq.Np))

Eq.N=6e6

Res=classRes()
Res.p0= 100.
Res.Pb0=fluido.Pb

Res.phi=0.25
Res.Sw0=0.20
Res.Swi=0.20
Res.Sor=0.25
Res.Vp0=Eq.N*fluido.Bo(Res.p0,Res.Pb0)/(1.-Res.Sw0-Res.Sg0)


# Petrophysics
Res.m=0.
Res.cr=50e-6
petroPhys=classPetrophys()
Swe=array([0.13,0.23,0.34,0.44,0.54,0.65,0.75])
petroPhys.Swe=Swe
petroPhys.kroe=0.7*((1.-Res.Sor-Swe)/(1.-Res.Sw0-Res.Sor))**3
petroPhys.krwe=0.075*((Swe-Res.Sw0)/(1.-Res.Sw0-Res.Sor))**6

# Aquifer
aquifero=classAquifero()
aquifero.WW=10*Eq.N # 400e6*0.27 # 35e6
aquifero.p0=Res.p0
aquifero.cr=Res.cr
aquifero.cw=1e-5
aquifero.pvt=fluido
aquifero.k=0.35
aquifero.L=1000.
aquifero.model='user_Np'
aquifero.Npe=Eq.Np.copy()
aquifero.Wee=Eq.We.copy()

Eq.res=Res
Eq.aquifero=aquifero
Eq.pvt=fluido
Eq.petroPhys=petroPhys

#-----------------------------------------------------------------------
# COMPUTE GAS PRODUCTION FROM EXP PRESSURE
#-----------------------------------------------------------------------
out=classHist()
out=Eq.runHist(1.,Res.p0)
savetxt('Fp.txt',out.Fp)
mksize=8

figure()
plot(Eq.Np/Eq.N,array(out.p)/Res.Pb0,'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
plot(Eq.Np/Eq.N,array(out.Pb)/Res.Pb0,'ko-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
legend(('reservoir pressure','bubble-point pressure'),loc='upper center')
xlabel('Oil Recovery Factor ($N_p/N$)')
ylabel('Pressure [units of $P_{b1}$]')
grid(1)
subplots_adjust(bottom=fbottom,left=fleft)
savefig('./article/python/matbal_p.pdf',dpi=600)

figure()
plot(Eq.Np/Eq.N,Eq.Wp/Eq.N,'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
plot(Eq.Np/Eq.N,Eq.Wi/Eq.N,'ko-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
plot(Eq.Np/Eq.N,Eq.We/Eq.N,'kx-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
xlabel('Recovery Factor ($N_p/N$)')
ylabel('Water Volume [units of $N$]')
legend(('production - $W_p$','injected - $W_i$','aquifer influx - $W_e$'),loc='upper left')
grid(1)
subplots_adjust(bottom=fbottom,left=fleft,right=0.85)
savefig('./article/python/matbal_water.pdf',dpi=600)

figure()
plot(Eq.Np/Eq.N,Eq.Gp/Eq.N/fluido.Rs(Res.p0,Res.Pb0),'ko-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
#legend(loc='upper center')
xlabel('Oil Recovery Factor ($N_p/N$)')
ylabel('Gas Recovery Factor ($G_p/N/R_s$)')
grid(1)
subplots_adjust(bottom=fbottom,left=fleft,right=0.85)
savefig('./article/python/matbal_gas.pdf',dpi=600)

#figure()
#plot(Eq.Np*1e-6,array(Eq.Wi)*1e-6,'bx-')
#plot(Eq.Np*1e-6,array(out.We)*1e-6,'gx-')
#legend(('Wi','We'),'upper left')
#xlabel('Cumulative oil production (1e6 std m3)')
#ylabel('Water Volume (1e6 std m3)')

figure()
plot(Eq.Np/Eq.N,out.So,'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
plot(Eq.Np/Eq.N,out.Sg,'kx-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
plot(Eq.Np/Eq.N,out.Sw,'ko-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
xlabel('Oil Recovery Factor ($N_p/N$)')
ylabel('Fluid Saturation')
ylim([0,1])
grid(1)
legend(('oil','gas','water'),loc='upper right')
subplots_adjust(bottom=fbottom,left=fleft)
savefig('./article/python/matbal_S.pdf',dpi=600)

figure()
plot(range(len(out.iter)),out.iter,'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
xlabel('Time steps')
ylabel('Newtonian Iterations')
ylim([0,16])
xlim([0,25])
grid(1)
#legend(('oil','gas','water'),loc='upper right')
subplots_adjust(bottom=fbottom,left=fleft)
savefig('./article/python/matbal_iter.pdf',dpi=600)

figure()
semilogy(range(len(out.iter)),abs(array(out.residual)),'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
xlabel('Time steps')
ylabel('Residual [units of porous volume]')
ylim([1e-15,1e-5])
xlim([0,25])
grid(1)
#legend(('oil','gas','water'),loc='upper right')
subplots_adjust(bottom=fbottom,left=fleft)
savefig('./article/python/matbal_residual.pdf',dpi=600)


#-----------------------------------------------------------------------
show()
