from pylab import *
from numpy.linalg import eig
from scipy.sparse.linalg import eigsh
from scipy.linalg import solve
import numpy as np
from scipy import stats
from scipy import integrate

def Ln(n, L, k, p):
    x = linspace(0.,L,n)
    dx = L/n
    
    # wn=1
    # k=(1.0 + 0.7*sin(wn*pi*x))*1e-2
    
    kph = zeros(n-1)
    for i in range(n-1):
        kph[i]=2.0*k[i]*k[i+1]/(k[i]+k[i+1])/dx**2
    
    kmh = zeros(n)
    for i in range(1,n):
        kmh[i]=2.0*k[i]*k[i-1]/(k[i]+k[i-1])/dx**2
    
    A=zeros([n,n])
    for i in range(1,n-1):
        A[i,i-1] = kmh[i]
        A[i,i] = (-kph[i]-kmh[i])
        A[i,i+1] = kph[i]
    
    A[0,0] = -kph[0]
    A[0,1] = kph[0]

    A[n-1,n-2] = kmh[n-1]
    A[n-1,n-1] =-kmh[n-1] #-A[n-1,n-2]
    
    return A

def Ld(n, L, k, p):
    x = linspace(0.,L,n)
    dx = L/n
    
    # wn=1
    # k=(1.0 + 0.7*sin(wn*pi*x))*1e-2
    
    kph = zeros(n-1)
    for i in range(n-1):
        kph[i]=2.0*k[i]*k[i+1]/(k[i]+k[i+1])/dx**2
    
    kmh = zeros(n)
    for i in range(1,n):
        kmh[i]=2.0*k[i]*k[i-1]/(k[i]+k[i-1])/dx**2
    
    A=zeros([n,n])
    for i in range(1,n-1):
        A[i,i-1] = kmh[i]
        A[i,i] = (-kph[i]-kmh[i])
        A[i,i+1] = kph[i]
       
    return A
				
#print(eigs[1:10:1]/(pi*arange(1,10,1)/L)**2)

#for n in [5,10,50,100,250,500,1000,2000]:
L = 1.
wn = 3
f=1.0
n=100
x = linspace(0.,L,n) 
khe=(1.0 + 0.7*sin(wn*pi*x))*1e-2 #np.random.random(n)
kho = L/trapz(1.0/khe,x)*ones(n)


p = linspace(0.,1.,n)
dt=0.1
nt=1000

pho0=[p[0]]
for i in range(nt):
	Lap = Ld(n, L=L, k=kho, p=p)
	a = np.identity(n) - Lap*dt
	a[0,0] = 1.
	a[0,1] = 0.
	a[n-1,n-2] = 0.
	a[n-1,n-1] = 1.
	p = solve(a, p)
	pho0.append(p[0])

pho0=[p[0]]
for i in range(nt):
	Lap = Ln(n, L=L, k=kho, p=p)
	a = np.identity(n) - Lap*dt
	p = solve(a, p)
	pho0.append(p[0])

p = linspace(0.,1.,n)
phe0=[p[0]]
for i in range(nt):
	Lap = Ld(n, L=L, k=khe, p=p)
	a = np.identity(n) - Lap*dt
	a[0,0] = 1.
	a[0,1] = 0.
	a[n-1,n-2] = 0.
	a[n-1,n-1] = 1.
	p = solve(a, p)
	phe0.append(p[0])

phe0=[p[0]]
for i in range(nt):
	Lap = Ln(n, L=L, k=khe, p=p)
	a = np.identity(n) - Lap*dt
	p = solve(a, p)
	phe0.append(p[0])
	
figure(1)
plot(range(nt+1), array(pho0), 'b')
plot(range(nt+1), array(phe0), 'r')

figure(2)
plot(x, p)
show()
