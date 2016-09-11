#!/usr/bin/env python

from pylab import *

def dpdf(x,y):
	r=zeros(len(x))
	r[0]=x[0]*y[0]
	for i in range(1,len(x)):
		r[i]=0.5*(y[i]+y[i-1])/(x[i]-x[i-1])
	
	return r

mul=[]
sdl=[]
for i in range(10):
	x=load('z'+str(i)+'.dat',skiprows=2)
	x[:,1]=x[:,1]*0.01
	E=dot(x[:,0],x[:,1])
	V=dot(x[:,1],(x[:,0]-E)**2)
	mu=log(E)-0.5*log(1+V/E**2)
	sd=sqrt(log(1.+V/E**2))
	xp=linspace(x[0,0],x[-1,0],200)
	logn=1./(xp*sd*sqrt(2*pi))*exp(-(log(xp)-mu)**2/(2.*sd**2))
#	logn=dpdf(xp,logn)
	x[:,1]=dpdf(x[:,0],x[:,1])
	print str(i)+' --------------------'
	print 'sum(P) = ',sum(x[:,1])
	print 'sum(logn) = ', sum(logn)
	print 'E = ',E
	print 'V = ',V
	print 'mu = ',mu
	print 'sd = ',sd
	print 'k = ',max(x[:,1])
	print '--------------------'
	figure(i+1)
	plot(x[:,0],x[:,1],'o-r',xp,logn,'-b')
	savefig('img/z'+str(i)+'.png')
	mul.append(mu)
	sdl.append(sd)

fid=open('fitting.dat','w')

fid.write('mu\n(')
for x in mul:
	fid.write('\n'+str(x))

fid.write('\n);\n\nsd\n(')
for x in sdl:
	fid.write('\n'+str(x))
fid.write(');')
fid.close()
show()
