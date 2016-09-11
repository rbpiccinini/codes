#!/usr/bin/env python

from pylab import *

D=0.0098
print '(x/D)max = ', 0.5/D
print '(y/D)max = ', 0.075/D


path=['coarseGrid/','fineGrid/']
locs=arange(0,30,5) #[5,20,40]
for n in locs:
	file=path[0]+'sets/5000/Profiles'+str(n)+'_U.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	u=loadtxt(file,delimiter=' ')[:,1]
	
	figure(n)
	plot(y/D,u,'-',label='x/D = '+str(n))
	xlabel('y/D')
	ylabel('Axial Velocity - U (m/s)')
	title('Axial Position - x/D = '+str(n))

	
	file=path[1]+'sets/5000/Profiles'+str(n)+'_U.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	u=loadtxt(file,delimiter=' ')[:,1]
	plot(y/D,u,'-',label='x/D = '+str(n))
	
	legend(('coarse','fine'))

show()
