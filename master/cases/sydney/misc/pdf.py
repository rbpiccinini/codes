#!/usr/bin/env python

from pylab import *
import scipy.special as sp
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


def estimate_plot(a):
	de=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,0]*1e-6
	pe=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,1]*1e6
	
	#E=a[0]
	#V=a[1]
	
	DD=a[0]
	q=a[1]
	
	#mu=log(E)-0.5*log(1+V/E**2)
	#sd=sqrt(log(1.+V/E**2))
	
	dn=linspace(min(de),max(de),200)
	#logn=1./(dn*sd*sqrt(2*pi))*exp(-(log(dn)-mu)**2/(2.*sd**2))
	#clogn=0.5*sp.erfc(-(log(dn)-mu)/sqrt(2)/sd)
	rl=q*dn**(q-1.)/DD**q*exp(-(dn/DD)**q)
	
	#DD=12.
	#q=3.
	#logn=q*xp**(q-1.)/DD**q*exp(-(xp/DD)**q)
	#clogn=1.-exp(-(xp/DD)**q)
	
	#pn=logn
	pn=rl
	cpn=spi.cumtrapz(pn,x=dn)
	z=array([0])
	cpn=hstack((z,cpn))
	
	cpe=spi.cumtrapz(pe,x=de)
	cpe=hstack((z,cpe))
	
	
	print ' --------------------'
	print 'sum(pn) = ', cpn[-1]
	print 'sum(pe) = ', cpe[-1]
#	print 'E = ',E
#	print 'V = ',V
#	print 'mu = ',mu
#	print 'sd = ',sd
	print 'DD = ',DD
	print 'q  = ',q
	print 'smd  = ',sum(pn*dn**3)/sum(pn*dn**2)
	print 'smde = ',sum(pe*de**3)/sum(pe*de**2)
	print '--------------------'
	
	figure()
	plot(dn*1e6,pn,'-k')
	plot(de*1e6,pe,'sr')
	ylabel('Rosin-Rammler pdf')
	xlabel('Droplet diameter ($\mu m$)')
	legend(('fitted', 'experiment'))
	grid(1)
	


def estimate(a):
	de=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,0]*1e-6
	pe=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,1]*1e6
	
	#E=a[0]
	#V=a[1]
	
	DD=a[0]
	q=a[1]
	
	#mu=log(E)-0.5*log(1+V/E**2)
	#sd=sqrt(log(1.+V/E**2))
	
	dn=linspace(8e-6,50e-6,200)

	#logn=1./(dn*sd*sqrt(2*pi))*exp(-(log(dn)-mu)**2/(2.*sd**2))
	#clogn=0.5*sp.erfc(-(log(dn)-mu)/sqrt(2)/sd)
	rl=q*dn**(q-1.)/DD**q*exp(-(dn/DD)**q)
	
	#DD=12.
	#q=3.
	#logn=q*xp**(q-1.)/DD**q*exp(-(xp/DD)**q)
	#clogn=1.-exp(-(xp/DD)**q)
	
	#pn=logn
	pn=rl
	cpn=spi.cumtrapz(pn,x=dn)
	z=array([0])
	cpn=hstack((z,cpn))
	
	cpe=spi.cumtrapz(pe,x=de)
	cpe=hstack((z,cpe))
	
	smd=sum(pn*dn**3)/sum(pn*dn**2)
	smde=sum(pe*de**3)/sum(pe*de**2)
	
	print ' --------------------'
	print 'sum(pn) = ', cpn[-1]
	print 'sum(pe) = ', cpe[-1]
#	print 'E = ',E
#	print 'V = ',V
#	print 'mu = ',mu
#	print 'sd = ',sd
	print 'DD = ',DD
	print 'q  = ',q
	print 'smd  = ',smd
	print 'smde = ', smde
	print '--------------------'
	
	delta=interp(dn,de,pe)-pn
	
	return (smd-smde)**2 + 1e7*sum(delta**2)/max(dn)**2

def estimatelogn(a):
	de=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,0]*1e-6
	pe=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,1]*1e6
	
	E=a[0]
	V=a[1]
	
	mu=log(E)-0.5*log(1+V/E**2)
	sd=sqrt(log(1.+V/E**2))
	
	dn=linspace(8e-6,50e-6,200)

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
	pn=pn/cpn[-1]
	cpn=cpn/cpn[-1]
	
	cpe=spi.cumtrapz(pe,x=de)
	cpe=hstack((z,cpe))
	pe=pe/cpe[-1]
	cpe=cpe/cpe[-1]
	
	smd=sum(pn*dn**3)/sum(pn*dn**2)
	smde=sum(pe*de**3)/sum(pe*de**2)
	
	print ' --------------------'
	print 'sum(pn) = ', cpn[-1]
	print 'sum(pe) = ', cpe[-1]
	print 'E = ',E
	print 'V = ',V
	print 'mu = ',mu
	print 'sd = ',sd
	print 'smd  = ',smd
	print 'smde = ', smde
	print '--------------------'
	
	delta=interp(dn,de,pe)-pn
	
	return (smd-smde)**2
	
def estimatelogn_plot(a):
	de=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,0]*1e-6
	pe=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,1]*1e6
	
	E=a[0]
	V=a[1]
	
	mu=log(E)-0.5*log(1+V/E**2)
	sd=sqrt(log(1.+V/E**2))
	
	dn=linspace(8e-6,50e-6,200)

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
	pn=pn/cpn[-1]
	cpn=cpn/cpn[-1]
	
	cpe=spi.cumtrapz(pe,x=de)
	cpe=hstack((z,cpe))
	pe=pe/cpe[-1]
	cpe=cpe/cpe[-1]
	
	smd=sum(pn*dn**3)/sum(pn*dn**2)
	smde=sum(pe*de**3)/sum(pe*de**2)
	
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
	
	figure()
	plot(dn*1e6,pn,'-k')
	ylabel('lognormal pdf')
	xlabel('Droplet diameter ($\mu m$)')
	legend(('fitted', 'experiment'))
	grid(1)
	
	return (smd-smde)**2

#out=opt.fmin_cg(estimate,[7e-6,2.],gtol=1e-07,maxiter=100)
#estimate_plot(out)
#estimate_plot([7.5e-6,1.36])
#show()

de=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,0]*1e-6
pe=loadtxt('exp/bc_x_0p5/pdf.csv',delimiter=',')[:,1]*1e6

E=15e-6
V=(2.5e-6)**2

print E, V
mu=log(E)-0.5*log(1.+V/E**2)
var=sqrt(log(1.+V/E**2))

print mu,var
#out=opt.fmin_cg(estimatelogn,[E,V],gtol=1e-07,maxiter=100)
estimatelogn_plot([E,V])
show()

