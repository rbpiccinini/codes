import scipy.special
#from sympy import *
import numpy as np
import scipy as sp
from pylab import *

def integrandoA(x,m,n,d,f,k):
	dx=np.gradient(x)
	if n==0:
		omegan=np.ones(len(x))
	else:
		omegan=np.sqrt(2.)*np.cos(n*np.pi*x)
	if m==0:
		omegam=np.ones(len(x))
	else:
		omegam=np.sqrt(2.)*np.cos(m*np.pi*x)
	
	domegan=np.gradient(k*np.gradient(omegan,dx),dx)
#	r2=k[-1]*omegam[-1]*gradient(omegan,dx)[-1]-k[0]*omegam[0]*gradient(omegan,dx)[0]
	r=np.trapz(omegam*domegan,x=x)
	return r

def integrandoB(x,m,n,phi,ct):
	if n==0:
		v=np.ones(len(x))
	else:
		v=np.sqrt(2.)*np.cos(n*np.pi*x)
	if m==0:
		u=np.ones(len(x))
	else:
		u=np.sqrt(2.)*np.cos(m*np.pi*x)
	r=-np.trapz(u*phi*ct*v,x=x)
	return r

def eigenfunction(x,v,i):
	Psi=0.0
	v=array(v[i])
	for n in range(len(v)):
		if n==0:
			omegan=np.ones(len(x))
		else:
			omegan=np.sqrt(2.)*np.cos(n*np.pi*x)
		Psi = Psi + omegan*v[n] 
	return Psi
#
#omegan=sqrt(2)*cos(n*pi*x)
#omegam=sqrt(2)*cos(m*pi*x)
#k=kk*(1+d*cos(x*pi*f))

x=np.linspace(0.,1.,500)

beta=1.0
k0=1.0
w0=10.0

#k=k0*exp(2.0*beta*x)
#w=w0*exp(2.0*beta*x)
khe=k0*(5.0+4.0*sin(2*pi*x))
k=khe
w=1.0

N=4
mu=1.0
k=k/mu
phi=1.0
ct=w
d=0.0
f=1

A=np.zeros([N,N])
for m in range(N):
	for n in range(N):
		A[n,m]=integrandoA(x,m,n,d,f,k)

B=np.zeros([N,N])
for m in range(N):
	for n in range(N):
		B[n,m]=integrandoB(x,m,n,phi,ct)

whe, v = sp.linalg.eig(A,b=B)
whe = np.sort(whe.real)

print('w=\n',np.sqrt(whe))
#print('v=\n',v)

#print('A=\n',A)
#print('B=\n',B)



##############################

keff=len(k)/(np.sum(1.0/k))
w=1.0

print('\nkeff=\n',keff)
#
figure()
plot(x,k)

k=keff

#N=30 
mu=1.0
k=k/mu
phi=1.0
ct=w
d=0.0
f=1

A=np.zeros([N,N])
for m in range(N):
	for n in range(N):
		A[n,m]=integrandoA(x,m,n,d,f,k)

B=np.zeros([N,N])
for m in range(N):
	for n in range(N):
		B[n,m]=integrandoB(x,m,n,phi,ct)

who, v = sp.linalg.eig(A,b=B)
who = np.sort(who.real)

print('\nw=\n',np.sqrt(who))

print('\ndelta =',(whe/who))
print('\nv =',v[0])

figure()
plot(range(len(whe)),np.sqrt(whe),'bs')
plot(range(len(who)),np.sqrt(who),'rs')

#figure()
#plot((whe)/who)

figure()
plot(x,eigenfunction(x,v,0))
plot(x,eigenfunction(x,v,1))
plot(x,eigenfunction(x,v,2))
plot(x,eigenfunction(x,v,3))
xlabel('distance')
ylabel('eigenfunction')

show()