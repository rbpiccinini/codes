#!/usr/bin/env python

from pylab import *

pathCFD='/home/piccinini/master/run/spraybench/2D/gasflow/simple/RANS.lowRE/sets/4150/'
#pathCFD='/home/rodrigo/master/run/spraybench/3D/RANS.kepsilon/sets/4550/'
#pathCFD='/home/rodrigo/OpenFOAM/run/spraybench/3D/rans/sets/3700/'
pathEXP='../exp/gasflow/'


locs=[2,25,50,100,200,300,400]
figure(1)
for i in locs:
	x=loadtxt(pathEXP+'z'+str(i)+'mm.dat',skiprows=2)
	x[:,0]=1e-3*x[:,0]
	y=loadtxt(pathCFD+'z'+str(i)+'mm_U.xy')
	figure()
	plot(x[:,0],x[:,1],'ro',y[:,0],y[:,3],'b-',-y[:,0],y[:,3],'b-')
	title('z = '+str(i)+'mm')
	savefig(pathCFD+'../../imgs/'+'z'+str(i)+'mm.png')

figure(2)
for i in locs:
	x=loadtxt(pathEXP+'z'+str(i)+'mm.dat',skiprows=2)
	x[:,0]=1e-3*x[:,0]
	y=loadtxt(pathCFD+'z'+str(i)+'mm_k.xy')
	figure()
	plot(x[:,0],x[:,2],'rx-',x[:,0],x[:,4],'gx-',y[:,0],sqrt(2./3.*y[:,1]),'k-')
	xlim([0,0.1])
	title('z = '+str(i)+'mm')
	legend(('Urms','Vrms','sqrt(2/3k)'))
	savefig(pathCFD+'../../imgs/'+'z'+str(i)+'mm_k.png')

#show()
