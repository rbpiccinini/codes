#-------------------------------------------------------------------------------
# Name:        Solve diffusion Eq. in 1D (linear solution)
#              SI units are used
# Purpose:
#
# Author:      Rodrigo, bfv8
#
# Created:     15/02/2013
#-------------------------------------------------------------------------------
#!/usr/bin/env python

from pylab import *

# Define domain size
np=250

# Define diffusion coefficient
k=100e-15
phi=0.2
ct=1e-2/1e5
mu=1e-3
eta=k/mu/phi/ct
print 'eta = ', eta


# Define reservoir limit
R2=5000.
R1=0.

# Define points
#r=logspace(log10(R1),log10(R2),np)
r=linspace(R1,R2,np)
dr=r[2]-r[1]

# Define time
dt=3600.*24.*30 # 1 month
nt=12

# Define Initial conditions
p=zeros(np)
p0=1e6
p[0]=p0

# Build linear system
C=eta*dt/dr**2
a=zeros([np,np])
for i in range(1,np-1):
	a[i,i]=-(1.+2*C)
	a[i,i-1]=C
	a[i,i+1]=C

# Apply boundary conditions

# Internal BC of fixed value
a[0,0]=-1. ##C-1. ##C/2-1.
a[0,1]=0. ##-2*C ##-C
a[0,2]=0. ## C    ##C/2

# External BC of zero gradient
a[np-1,np-2]=2.*C
a[np-1,np-1]=-(1.+2.*C)

for tt in range(nt):
	p=dot(inv(a),-1.*p)

	figure(1)
	plot(r,p/1e5, 'k-')
#	ylim([1e-6,1.])

grid(1)
show()
