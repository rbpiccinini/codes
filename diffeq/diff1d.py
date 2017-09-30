from pylab import *
from numpy.linalg import eigvals

n = 200
L = 1.
dx = L/n
k = 1.0*ones(n)

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
A[0,1] = -A[0,0] 
A[n-1,n-2] = kmh[n-1]+k[n-1]
A[n-1,n-1] = -A[n-1,n-2]
A = -A/dx**2

print(sort(eigvals(A))[:10])
