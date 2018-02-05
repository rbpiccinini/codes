from pylab import *
from numpy.linalg import eig
from scipy.sparse.linalg import eigsh
import numpy as np
from scipy import stats
from scipy import integrate

def householder(A):
	n = len(A)
	for k in range(n-2):
		alpha = np.sign(A[k+1,k])*np.sqrt(sum(A[k+1:,k]**2))
		r = np.sqrt(0.5*(alpha**2-A[k+1,k]*alpha))
		v = zeros(n)
		v[k+1] = (A[k+1,k]-alpha)/(2.0*r)
		v[k+2:] = A[k+2:,k]/(2.0*r)
		P = np.identity(n)-2.0*np.outer(v,v)
		A = dot(P,dot(A,P))

	a = zeros(n-1)
	b = zeros(n)
	c = zeros(n-1)
	
	# a \ b \ c
	for i in range(1,n-1):
		a[i-1] = A[i,i-1]
		b[i] = A[i,i]
		c[i] = A[i, i+1]
	
	b[0] = A[0,0]
	b[n-1] = A[n-1,n-1]

	c[0] = A[0,1]
	a[n-2] = A[n-1,n-2]
	
	A = np.diag(a,-1)+np.diag(b, 0)+np.diag(c,1)
	
	return A
	 

def eigCalc(n, L, k):
	x = linspace(0.,L,n)
	dx = L/n
	
	# wn=1
	# k=(1.0 + 0.7*sin(wn*pi*x))*1e-2

	kph = zeros(n-1)
	for i in range(n-1):
		kph[i]=2.0*k[i]*k[i+1]/(k[i]+k[i+1])
	
	kmh = zeros(n)
	for i in range(1,n):
		kmh[i]=2.0*k[i]*k[i-1]/(k[i]+k[i-1])
		
	A=zeros([n,n])
	a = zeros(n-1)
	b = zeros(n)
	c = zeros(n-1)
	
	# a \ b \ c
	for i in range(1,n-1):
		a[i-1] = kmh[i] 
		b[i] = -kph[i]-kmh[i]
		c[i] = kph[i]

#		A[i,i-1] = kmh[i]
#		A[i,i] = -kph[i]-kmh[i]
#		A[i,i+1] = kph[i] 

		b[0] = -kph[0]
		c[0] = kph[0]

		a[n-2] = kmh[n-1]
		b[n-1] =-kmh[n-1] #-A[n-1,n-2]


	a = -a/dx**2
	b = -b/dx**2
	c = -c/dx**2
	
	
	A = np.diag(a,-1)+np.diag(b, 0)+np.diag(c,1)
	
	nest=100
	w, v = eigsh(A, k=nest, which='SM')
	k = argsort(w)
	
	
#     reconstructing k(x)
#    s = zeros(n)
#    for i in range(nest):
#        s = s + w[i]*v[:,i]
#    s = integrate.cumtrapz(s, dx=dx, initial=0.)
#    kest = -s*dx/gradient(sum(v[:,1:],axis=1))
	kest = recK(A, dx)
	
	return w, v, kest, A

def recK(A, dx):
	#     reconstructing k(x)
	n = len(A)
	
	k=zeros(n)
	kph = zeros(n)
	kmh = zeros(n)

	for i in range(1,n-1):
		k[i] = 0.5*A[i,i]
	k[0] = A[0,0]
	k[n-1] = A[n-1,n-1] #A[n-1,n-1]
#	for i in range(n-1):
#		kph[i]= A[i, i+1]
#		kmh[i+1]= A[i+1, i]

#	kph[n-1] = kmh[n-1]
#	kmh[0] = kph[0]
#
#	k = -2.0*kph*kmh/(kph+kmh)*dx**2

	return k*dx**2
	

def calcvho(x,n):
	L=max(x)
	vho=zeros([len(x),n])
	for i in range(n):
		vho[:,i] = cos(i*pi/L*x)
		vho[:,i] = vho[:,i]/sqrt(dot(vho[:,i],vho[:,i]))
	return -vho
				
#print(eigs[1:10:1]/(pi*arange(1,10,1)/L)**2)

#for n in [5,10,50,100,250,500,1000,2000]:
L = 1.
wn = 3
f=1.0

n = 500
m = 200
x = linspace(0.,L,n) 
khe=(1.0 + 0.7*sin(wn*pi*x))*1e-2 #np.random.random(n)
kho = L/trapz(1.0/khe,x)
w, v, kest, A = eigCalc(n,L=L, k=khe)
who, vho, kestho, Aho = eigCalc(n,L=L, k=kho*ones(n))

Arec = Aho #zeros([n,n])
nest = 100
for i in range(nest):
	Arec = Arec + (w[i]-who[i])*outer(vho[:,i],vho[:,i])
Arec = householder(Arec)
wrec, vrec = eigsh(Arec, k=10, which='SM')
krec = recK(Arec, L/n)

intveta=ones(len(w)-1)
for i in range(1,len(w)):
	intveta[i-1] = trapz(khe*gradient(v[:,i])**2/i**2.,x=x)/kho
	
keff=zeros(len(w))

keff[1:]= w[1:]*L**2/(pi*arange(1,len(w)))**2
for i in range(len(keff)):
	print('{:2d} \t keff = {:8.5f} \t w = {:8.5e}'.format(i+1, keff[i], w[i]))
#    print(i+1,'k = ',keff[i], '\t\tw = ',w[i])



