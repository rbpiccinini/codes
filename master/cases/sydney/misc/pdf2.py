
#!/usr/bin/env python

from pylab import *
import scipy.special as sp
import numpy as np
import scipy.integrate as spi
import scipy.optimize as opt

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

def erfInv(y):
	a=0.147
	k=2./pi/a+0.5*(log(1.0-y**2))
	h=log(1.0 - y*y)/a
	x=sqrt(-k + sqrt(k*k - h))
	
	if (y < 0.0):
		x = -1.0*x
	return x
	

def sample(mu,sd):
	minValue=1.e-6
	maxValue=60.e-6
	
	a = sp.erf((minValue - mu)/sd)
	b = sp.erf((maxValue - mu)/sd)
	y=random()
	x = erfInv(y)*sd + mu
	x = min(max(exp(x), exp(minValue)), exp(maxValue))
	return log(x);
	
def estimatelogn(a,draw):
	de=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,0]*1e-6
	pe=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,1]*1e6
	
	E=a[0]
	V=a[1]**2
	
	mu=log(E)-0.5*log(1+V/E**2)
	sd=sqrt(log(1.+V/E**2))
	
	dn=linspace(0,max(de),200)

	pn=1./(dn*sd*sqrt(2*pi))*exp(-(log(dn)-mu)**2/(2.*sd**2))
	#clogn=0.5*sp.erfc(-(log(dn)-mu)/sqrt(2)/sd)
	
	#DD=12.
	#q=3.
	#logn=q*xp**(q-1.)/DD**q*exp(-(xp/DD)**q)
	#clogn=1.-exp(-(xp/DD)**q)
	
	#pn=logn
	cpn=spi.cumtrapz(pn,x=dn)
	z=array([0])
	cpn=hstack((z,cpn))
	
	cpe=spi.cumtrapz(pe,x=de)
	cpe=hstack((z,cpe))
	
	smd=sum(dot(diff(cpn),diff(dn))*dn[1:]**3)/sum(dot(diff(cpn),diff(dn))*dn[1:]**2)
	smde=sum(dot(diff(cpe),diff(de))*de[1:]**3)/sum(dot(diff(cpe),diff(de))*de[1:]**2)
	
	print ' --------------------'
	print 'sum(pn) = ', cpn[-1]
	print 'sum(pe) = ', cpe[-1]
	print 'E = ',E
	print 'V = ',V
	print 'SD = ', sqrt(V)
	print 'mu = ',mu
	print 'var = ',sd**2
	print 'sd = ',sd
	print 'smd  = ',smd
	print 'smde = ', smde
	print '--------------------'
	
	delta=interp(dn,de,pe)-pn
	
	if draw==1:
		fig = figure()
		ax = fig.add_subplot(1,1,1)
		subplots_adjust( left=0.2, bottom=0.15 )
		ax.plot(dn*1e6,pn,'-k')
		#ax.plot(de*1e6,pe,'s-r')
		ylabel('lognormal pdf ($1/m$)')
		xlabel('Droplet diameter ($\mu m$)')
		#legend(('fitted', 'experiment'))
		grid(0)
	
	return sum(delta**2)/dn[-1]**2 # (smd-smde)**2/smd**2


#mu=6.5e-6
#sd=sqrt(2.5e-11)

mu=10e-6
sd=sqrt(5e-11)

#figure()
#out=opt.fmin_cg(estimatelogn,[mu,sd],args=[0],gtol=1e-07,maxiter=100)
#mu=out[0]
#sd=out[1]
#print out

estimatelogn([mu,sd],1)

d=linspace(1e-6,1,1000)
for i in range(len(d)):
	d[i]=sample(mu,sd)*1e6

figure()
hist(d,50)
xlim([1,60])
print 'SMD =',sum(d**3)/sum(d**2)
print 'Mean = ',average(d)

show()
