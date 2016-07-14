##
## Solves 1D transient heterogeneous diffusion using GITT
##
import scipy.special
#from sympy import *
import numpy as np
import scipy as sp
from pylab import *

class classGITT():
	def __init__(self,x,w,k,N):
		self.x=x
		self.dx=np.gradient(x)
		self.w=w
		self.k=k
		self.N=N

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
				self.B[m,n] = -1.*np.trapz(self.omega(n)*self.w*self.omega(m),x=self.x)

	def solveEigenProb(self):
		self.calcA()
		self.calcB()
		self.eigenvalues, self.eigenvectors = sp.linalg.eig(self.A,b=self.B)
		self.eigenvalues=self.eigenvalues.real
		ind=np.argsort(self.eigenvalues)

		self.eigenvalues=self.eigenvalues[ind]
		self.eigenvectors=self.eigenvectors[ind]


	def Psi(i):
		Psi=0.0
		for n in range(self.N):
			Psi = Psi + self.omega(n)*self.eigenvectors[i][n] 
		return Psi

x=np.linspace(0.,1.,500)
beta=1.0
k0=1.0
w0=10.0

k=k0*exp(2.0*beta*x)
w=w0*exp(2.0*beta*x)
khe=k0*(5.0+4.0*sin(2*pi*x))
N=20

homog=classGITT(x=x,w=w,k=k,N=N)
homog.solveEigenProb()

print(np.sqrt(homog.eigenvalues)[1:11:1])

