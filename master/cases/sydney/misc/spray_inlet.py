#!/usr/bin/env python

from pylab import *
import numpy as np

def pv(T):
	return 101325.*np.exp(-3613.6637*(1./T-3.03951e-3))

rhol=780. #kg/m3
L=518e3 # J/kg
mpfuel=7. *1e-3/60.
mpair=135. *1e-3/60. # kg/s 
mpac= 5. *1e-3/60.

mpspray= mpfuel - mpac

cpair=1005. # J/kg/K
cpl=1.6733*cpair
kappal=0.161
Ti=298.- mpac*L/(mpair*cpair+mpac*cpl)

Yair=mpair/(mpair+mpac)
Yac=1-Yair

YO2air=0.233
YN2air=1.-YO2air

MMair=0.029 # kg/mol
MMac=0.058 # kg/mol

MMmix=1./(Yair/MMair+Yac/MMac)

R=8.314
Rmix=R/MMmix
Pi=1e5 # Pa
rhomix=Pi/Rmix/Ti

D=0.0098
mug=1.874e-4 *1e-3*1e2 #1.78e-5
mul=0.4013e-3
Reg=(mpair+mpac)*4/D/pi/mug
Rel=(mpfuel)*4/D/pi/mul

Umean=(mpair+mpac)/(pi*D**2)/rhomix
tnew=0.5/Umean

px=linspace(0,0.0049,50)
py=-467.61*px+2.347
for i in range(len(py)):
	py[i]=min(py[i],1.0)
py=py*31.

exp=loadtxt('exp/bc_x_0p5/Xac.csv',skiprows=1,delimiter=',')
y=exp[:,0]*D
Yac2=exp[:,1]*MMac/MMmix/100.
fluxInlet=trapz(Yac2*y,x=y)*2.*pi*26.4*rhomix
print 'fluxInlet [g/min] = ', fluxInlet*1e3*60.

print 'Ti = ',Ti
print 'mass_nozzle = ', mpair+mpspray
print 'mpspray = ', mpspray
print 'Yac = ',Yac
print 'YO2 = ', (1.-Yac)*YO2air
print 'YN2 = ', (1.-Yac)*YN2air
print 'Reg = ', Reg
print 'Rel = ', Rel
print 'Umean = ', Umean
print 'tnew = ', tnew

fig_width = 8.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 13,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
          'figure.figsize': fig_size}
rcParams.update(params)
plot(px/D,py)
xlabel('Radial position - y/D')
ylabel('Droplet injection velocity (m/s)')

# checking vapor pressure

figure()
Tv=280.#linspace(280.,330.,100)
Tb=329.
#print L*MMac/R,1./Tb
pv=np.exp(-L*MMac/R*(1./Tv-1./Tb))
#pv2=133.32*pow(10.0,(7.02447-1161.0/(224.0+Tv-273.0)))
#pv3=101325.*exp(-3613.6637*(1./Tv-3.03951e-3))
plot(Tv,pv,'b')
#plot(Tv,pv3,'r')
xlabel('T')
ylabel('pv')


figure()
Yf=0.035
Ys=pv
mp=5e-6*3.34e-5*5./360.

DD=1.24e-05
SMD=13.7e-6

Dd=60e-6
dmdt=-2.*pi*Dd*DD*rhomix*log((1.-Yf)/(1.-Ys))
md=rhol*pi/6.*Dd**3
tevap=-md/dmdt
N=2*rhomix*Dd*cpl/kappal*log((1.-Yf)/(1.-Ys))
print 'tevap,N = ',tevap,N
plot(Tv,tevap)
xlabel('T')
ylabel('time')
#show()
