from pylab import *
from numpy.linalg import eig
import numpy as np
from scipy import stats
from scipy import integrate

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
    for i in range(1,n-1):
        A[i,i-1] = kmh[i]
        A[i,i] = -kph[i]-kmh[i]
        A[i,i+1] = kph[i] 
    
    A[0,0] =-kph[0]-k[0]
    A[0,1] = kph[0]

    A[n-1,n-2] = kmh[n-1]
    A[n-1,n-1] =  -k[n-1]-kmh[n-1] #-A[n-1,n-2]
    
    A = -A/dx**2
    w, v = eigh(A)
    k = argsort(w)
    
#     reconstructing k(x)
    nest=10
    s = zeros(n)
    for i in range(1,nest):
        s = s - w[i]*v[:,i]
    s = integrate.cumtrapz(s, dx=dx, initial=0.)
    kest = s*dx/gradient(sum(v[:,1:],axis=1))
    
    return w[0:20:1], -v[:,0:20:1], kest

#print(eigs[1:10:1]/(pi*arange(1,10,1)/L)**2)

#for n in [5,10,50,100,250,500,1000,2000]:
L = 1.0
for n in [1000]:
    x = linspace(0.,L,n) 
    khe = 1.0*ones(n)
    kho = L/trapz(1.0/khe,x)
    w, v, kest = eigCalc(n,L=L, k=khe)


keff = w*L**2/(pi*arange(1,len(w)+1))**2
for i in range(len(keff)):
    print(i+1,'k = ',keff[i], '\t\tw = ',w[i])

print('Analytical = ', kho)
print('kmin = ',min(khe))
print('kmax = ',max(khe))
print('kest = ',mean(kest))

figure()
for i in range(3):
    plot(x,v[:,i])
      
figure()
plot(x,kest)
show()