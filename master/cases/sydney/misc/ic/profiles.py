#!/usr/bin/env python

from pylab import *

D=0.0098
print '(x/D)max = ', 0.5/D
print '(y/D)max = ', 0.075/D


path='fineGrid/'
rh=[]
locs=arange(15,50,5) #[5,20,40]
for n in locs:
	file=path+'sets/5000/Profiles'+str(n)+'_U.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	u=loadtxt(file,delimiter=' ')[:,1]
	ub=(u-3.)/(u[0]-3.)
	ih=argmin(abs(u-0.5*u[0]))
	rh.append(y[ih])
	print '(x/D,rh,error) =', n,rh[-1],min(abs(u-0.5*u[0]))/u[0]
	
	figure(1)
	plot(y/rh[-1],ub,'-',label='x/D = '+str(n))
	xlabel('y/D')
	ylabel('Axial Velocity - U (m/s)')
	title('Axial Position - x/D = '+str(n))
	legend()
	
	file=path+'sets/5000/Profiles'+str(n)+'_k.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	k=loadtxt(file,delimiter=' ')[:,1]
	figure(2)
	plot(y/rh[-1],2./3.*k/u[0]**2,'-',label='x/D = '+str(n))
	xlabel('y0.5/D')
	ylabel('Turbulent Kinetic Energy - k (m2/s2)')
	legend()
	title('Axial Position - x/D = '+str(n))
xlim([0,2.2])

figure(3)
S=(0.594595-1.39189)/(8.04545-18.9545)
B=sum(array(rh)+S*locs)/len(locs)

plot(locs,array(rh)/D,'x-')
plot(locs,S*locs+B,'r-')
show()
