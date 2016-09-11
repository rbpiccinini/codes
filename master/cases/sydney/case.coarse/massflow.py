#!/usr/bin/env python
from pylab import *

x=loadtxt('patchMassFlows_massFlow/0.232/massFlow',skiprows=1)
figure()
plot(x[:,0],x[:,1],label='nozzle')
legend()

figure()
plot(x[:,0],x[:,2],label='coflow')
legend()

figure()
plot(x[:,0],x[:,3],label='outlet')
legend()

figure()
plot(x[:,0],x[:,4],label='walls')
legend()

figure()
plot(x[:,0],sum(x[:,1:],axis=1), label='mass bal')
legend()
show()
