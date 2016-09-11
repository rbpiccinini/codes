#!/usr/bin/env python

from pylab import *

y=loadtxt('parcels.dat',delimiter='|',usecols=[1])
x=linspace(0,len(y)*2e-6,len(y))
plot(x,y)
ylabel('No of parcels in the domain')
xlabel('Time (s)')
show()
