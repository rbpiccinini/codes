##
## Solves 1D transient heterogeneous diffusion in a line
## using the Generalized Integral Transform Technique (GITT)
##
## Author: Rodrigo Piccinini
## Reference: Carolina P. Naveira-Cotta. Eigenfunction expansions for transient diffusion in heterogeneous media

##
# HEADER
##

import scipy.special
import scipy.stats
import numpy as np
import scipy as sp
from pylab import *

params = {'backend': 'ps',
#'text.latex.preamble': ['\usepackage[latin1]{inputenc}'],
'axes.labelsize': 14, # fontsize for x and y labels (was 10)
'axes.titlesize': 14,
'text.fontsize': 14, # was 10
'legend.fontsize': 14, # was 10
'xtick.labelsize': 14,
'ytick.labelsize': 14,
'text.usetex': True,
'figure.figsize': [6.9,8.0],
'font.family': 'serif',
'lines.linewidth': 1
}

rcParams.update(params)

class classGITT():
	def __init__(self,x,w,k,N,f,phi,alpha,beta):
		self.x=x
		self.dx=np.gradient(x)
		self.w=w
		self.k=k
		self.N=N
		self.f=f
		self.phi=phi
		self.alpha=alpha
		self.beta=beta

	def omega(self,i):
		if i==0:
			return np.ones(len(self.x))
		else:
			return np.sqrt(2)*np.cos(i*pi*self.x)

	def calcA(self):
		self.A=np.zeros([self.N,self.N])
		for n in range(self.N):
			omegan=self.omega(n)
			for m in range(self.N):
				omegam=self.omega(m)
				domegan=np.gradient(self.k*np.gradient(omegan,self.dx),self.dx)
				self.A[n,m]=np.trapz(omegam*domegan,x=self.x)
	
	def calcB(self):
		self.B=np.zeros([self.N,self.N])
		for n in range(self.N):
			for m in range(self.N):
				self.B[n,m] = -1.*np.trapz(self.omega(n)*self.w*self.omega(m),x=self.x)

	def solveEigenProb(self):
		self.calcA()
		self.calcB()
		self.eigenvalues, self.eigenvectors = sp.linalg.eig(self.A,b=self.B)
		self.eigenvalues=self.eigenvalues.real + 1e-10 # avoid close to zero negative values
		self.ind=np.argsort(self.eigenvalues)

	def Psi(self,i):
		Psi=0.0
		for n in range(self.N):
			Psi = Psi + self.omega(n)*self.eigenvectors[n,i]
		return Psi
	
	def gi(self,i):
		r1=self.phi[-1]*(self.Psi(i)[-1]-self.k[-1]*gradient(self.Psi(i),self.dx)[-1]/(self.alpha[-1]+self.beta[-1]))
		r0=self.phi[0]*(self.Psi(i)[0]-self.k[0]*gradient(self.Psi(i),self.dx)[0]/(self.alpha[0]+self.beta[0]))
		return r1-r0

	def Ni(self,i):
		return np.trapz(self.w*self.Psi(i)**2,x=self.x)
	
	def fi(self,i):
		return np.trapz(self.w*self.Psi(i)*self.f/self.Ni(i),x=self.x)

	def p(self,t):
		r=0.
		for i in range(self.N):
			r=r+self.Psi(i)*self.fi(i)*exp(-self.eigenvalues[i]*t)+self.gi(i)/self.eigenvalues[i]*(1.0-np.exp(-self.eigenvalues[i]*t))
		return r

##
# BODY
##

x=np.linspace(0.,1.,400)

#beta=3.0
#k0=1.0
#w0=10.0

#k=k0*exp(2.0*beta*x)
#w=w0*exp(2.0*beta*x)
#f=(1.0-np.exp(2*beta*(1.-x)))/(1.-np.exp(2.0*beta))


N=100

f=1.0
kho=0.01
khe1=kho*(1.0 + 0.2*sin(4*pi*x))
khe2=kho*(1.0 + 0.8*sin(4*pi*x))
w0=1.0

print('kho = ', 1.0/trapz(1.0/khe2,x))

homog=classGITT(x=x,w=w0,k=kho,N=N,f=f,phi=array([1.,0.]), alpha=array([1.,1.]), beta=array([0.,0.]))
homog.solveEigenProb()

heterog1=classGITT(x=x,w=w0,k=khe1,N=N,f=f,phi=array([1.,0.]), alpha=array([1.,1.]), beta=array([0.,0.]))
heterog1.solveEigenProb()

heterog2=classGITT(x=x,w=w0,k=khe2,N=N,f=f,phi=array([1.,0.]), alpha=array([1.,1.]), beta=array([0.,0.]))
heterog2.solveEigenProb()

eigho=np.sqrt(homog.eigenvalues)[homog.ind]
eighe1=np.sqrt(heterog1.eigenvalues)[heterog1.ind]
eighe2=np.sqrt(heterog2.eigenvalues)[heterog2.ind]
print(eigho)
print(eighe1)
print(eighe2)

figure()

subplot(211)
plot(homog.eigenvalues[homog.ind[1:30]],'kx:')
plot(heterog1.eigenvalues[heterog1.ind[1:30]],'s',markerfacecolor='None',markeredgewidth=1.25)
plot(heterog2.eigenvalues[heterog2.ind[1:30]],'o',markerfacecolor='None',markeredgewidth=1.25)
ylim([0.,100.])
xlabel('Eigenvalue number')
ylabel('Eiganvalue magnitude [1/s]')
legend(('homogeneous', 'low heterogeneity', 'high heterogeneity'),loc='upper left')

subplot(212)
plot(x,kho*ones(len(x)),'k-')
plot(x,khe1,'k--',markerfacecolor='None',markeredgewidth=1.25, markevery=10)
plot(x,khe2,'k:',markerfacecolor='None',markeredgewidth=1.25, markevery=10)
ylim([0.,0.02])
xlabel('Position [m]')
ylabel('diffusion cofficient [m^2/s]')
legend(('homogeneous', 'low heterogeneity', 'high heterogeneity'))
tight_layout()
show()

#figure()
#plot(x,homog.p(0.1),'b')
#show()
