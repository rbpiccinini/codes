from pylab import *
from scipy.integrate import odeint
from scipy.integrate import cumtrapz


fig_width = 7.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 18,
          'text.fontsize': 18,
          'legend.fontsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
          'figure.figsize': fig_size}
rcParams.update(params)

def calcCd(U,D):
	rho=1.2
	rhod=780.
	mu=1.8e-5
	Re=rho*U*D/mu
	
	return 24./Re*(1.+1./6.*Re**0.6667)
	
def dudt(U,t,D):
	rho=1.2
	rhod=780.
	mu=1.8e-5
	Re=rho*U*D/mu
	Cd=max(calcCd(U,D),0.424)
	tu=4./3.*rhod*D/(rho*U*Cd)
	return -U/tu
	
#D=array([5.,10.,20.,30])*1e-6
D=linspace(5,40,50)*1e-6
Urel=100.
rho=1.2
rhod=780.
mu=1.8e-5

Re=rho*Urel*D/mu
print max(Re)
Cd=24./Re*(1.+1./6.*Re**0.6667)
dudt1=calcCd(5.,D)
dudt2=calcCd(10.,D)
dudt3=calcCd(30.,D)


fig = figure()
# plot data
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
ax.plot(D*1e6,dudt1/max(dudt1),'k-')
ax.plot(D*1e6,dudt2/max(dudt2),'k--')
ax.plot(D*1e6,dudt3/max(dudt3),'k:')
#legend(('')).draw_frame(False)
xlabel('Droplet Diameter ($\mu m$)')
ylabel('Normalized Droplet Acceleration')

time=linspace(0,0.01,100)


U2 = odeint(dudt,10.,time,args=((10e-6,)))
U3 = odeint(dudt,10.,time,args=((20e-6,)))
U4 = odeint(dudt,10.,time,args=((30e-6,)))
U5 = odeint(dudt,10.,time,args=((40e-6,)))
print len(U2),len(time)
x2=cumtrapz(U2,x=time,axis=0)
x3=cumtrapz(U3,x=time,axis=0)
x4=cumtrapz(U4,x=time,axis=0)
x5=cumtrapz(U5,x=time,axis=0)
print len(x5)

fig=figure()
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
ax.plot(1e3*time,U2,'k--',label='$10\mu m$')
ax.plot(1e3*time,U3,'k:',label='$20\mu m$')
ax.plot(1e3*time,U4,'k-.',label='$30\mu m$')
ax.plot(1e3*time,U5,'k-',label='$40\mu m$')
ylim([-0.5,10])
xlabel('time (ms)')
ylabel('Velocity (m/s)')
legend(loc='upper right').draw_frame(False)

savetxt('drag.dat',vstack([1e3*time,U2.T,U3.T,U4.T,U5.T]).T,delimiter='\t')

fig=figure()
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
ax.plot(1e3*time[:-1],x2,'k--',label='$10\mu m$')
ax.plot(1e3*time[:-1],x3,'k:',label='$20\mu m$')
ax.plot(1e3*time[:-1],x4,'k-.',label='$30\mu m$')
ax.plot(1e3*time[:-1],x5,'k-',label='$40\mu m$')
xlabel('time (ms)')
ylabel('Position (m)')
#legend(loc='upper right').draw_frame(False)
show()
