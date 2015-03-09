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
fluido.pe=array([250.,170.,139.8,109.5,79.3,49.,41.,33.,25.,17.,9.,1.])
fluido.Boe=array([1.15,1.13,1.12,1.11,1.098,1.083,1.079,1.074,1.068,1.062,1.052,1.035])
fluido.Rse=array([60.,47.,41.,35.,28.,20.4,18.,15.5,13.,10.,5.,0.21])
fluido.Bge=array([0.004,0.0057,0.007,0.0093,0.0136,0.022,0.028,0.036,0.049,0.0738,0.1433,1.2443])
fluido.Bwe=1.02379*(1-3.45e-5*(fluido.pe - 97.45))
fluido.viscoe=array([5.81000,7.23000,8.05000,9.05000,10.23000,11.58000,12.17000,12.80000,14.12000,14.83000,15.57000,16.33000,17.12000,17.96000,18])
fluido.viscwe=1.0*ones(len(fluido.pe))
fluido.co=1.05e-4
fluido.Pb=48.
fluido.T=273.+57.

#---------------------------------------------------------------------------
# COMPUTING EQUATION COEFFICIENTS
# Np*(Bo2+(Rp-Rso2)*Bg2)+Wp*Bw2=N*(Eo+m*Eg+Efw)+(Wi+We)*Bw2+Gi*Bg2
#---------------------------------------------------------------------------
Eq=classEBM()
## prod=genfromtxt('production.dat',delimiter='\t',skiprows=2,usecols = (0,1,2,3,4), converters = {4: lambda s: float(s or 0)})
## prod=genfromtxt('production_datum_850.dat',delimiter='\t',skiprows=3)
prod=genfromtxt('production.dat',delimiter='\t',skiprows=2)

Eq.dt=prod[:,0]
Eq.p=prod[:,1]
Eq.Np=prod[:,2]
Eq.Gp=prod[:,3]
Eq.Wp=prod[:,4]
Eq.Wi=prod[:,5]
Eq.Gi=prod[:,6]
Eq.We=prod[:,7]

Eq.N=4.89055e6

Res=classRes()
Res.p0= 97.22 # p0 = 97.46 kgf/cm2 at datum z=-850 m
Res.Pb0=fluido.Pb

Res.phi=0.2092435
Res.Sw0=0.35
Res.Vp0=Eq.N*fluido.Bo(Res.p0,Res.Pb0)/(1.-Res.Sw0-Res.Sg0)


# Petrophysics
Res.m=0.
Res.cr=21.76e-6

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

#-----------------------------------------------------------------------
# COMPUTE GAS PRODUCTION FROM EXP PRESSURE
#-----------------------------------------------------------------------
out=classHist()
out=Eq.runHist(1.,Res.p0)
savetxt('Fp.txt',out.Fp)
mksize=8

print len(Eq.p), len(Eq.Np)
figure()
plot(Eq.Np/Eq.N,array(out.p)/Res.Pb0,'ks',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
plot(Eq.Np/Eq.N,Eq.p/Res.Pb0,'r-')
plot(Eq.Np/Eq.N,array(out.Pb)/Res.Pb0,'ko',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
legend(('reservoir pressure','3D model','bubble-point pressure'),loc='upper center')
xlabel('Oil recovery factor ($N_p/N$)')
ylabel('Pressure [units of $P_{b1}$]')
grid(1)
subplots_adjust(bottom=fbottom,left=fleft)
ylim([0.,2.0])
savefig('./article/python/matbal_p.pdf',dpi=600)

figure()
plot(Eq.Np/Eq.N,Eq.Wp/Eq.N,'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
plot(Eq.Np/Eq.N,Eq.Wi/Eq.N,'ko-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
plot(Eq.Np/Eq.N,Eq.We/Eq.N,'kx-',linewidth=0.5,markersize=mksize,markerfacecolor='k',markeredgewidth=2)
xlabel('Recovery factor ($N_p/N$)')
ylabel('Water volume [units of $N$]')
legend(('production - $W_p$','injected - $W_i$','aquifer influx - $W_e$'),loc='upper left')
grid(1)
subplots_adjust(bottom=fbottom,left=fleft,right=0.85)
savefig('./article/python/matbal_water.pdf',dpi=600)

figure()
plot(Eq.Np/Eq.N,Eq.Gp/Eq.N/fluido.Rs(Res.p0,Res.Pb0),'ko-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
#legend(loc='upper center')
xlabel('Oil recovery factor ($N_p/N$)')
ylabel('Gas recovery factor ($G_p/N/R_s$)')
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
xlabel('Oil recovery factor ($N_p/N$)')
ylabel('Phase saturation')
ylim([0,1])
grid(1)
legend(('oil','gas','water'),loc='upper right')
subplots_adjust(bottom=fbottom,left=fleft)
savefig('./article/python/matbal_S.pdf',dpi=600)

figure()
plot(range(len(out.iter)),out.iter,'ks-',linewidth=0.5,markersize=mksize,markerfacecolor='#ffffff',markeredgewidth=2)
xlabel('Time steps')
ylabel('Number of iterations')
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
