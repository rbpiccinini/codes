#!/usr/bin/env python

from pylab import *

def mu(T,Tc,pc):
	fp=1.141824
	Tr=T/Tc
	muEps=0.807*Tr-0.357*exp(-0.449*Tr)+0.34*exp(-4.058*Tr)+0.018
	Eps=(Tc*8314)**(1./6.)/sqrt(60.1)/pc**0.666*6.02e26**0.33
	return muEps*fp/Eps
	



print mu(273+40,508.3,48.2e5)
