#!/usr/bin/env python
from pylab import *

def main(outputdir,n):
	b=loadtxt(outputdir+'../b.dat')
	modeLz=zeros(b.shape[1])
	Lz=zeros([b.shape[0],n])
	for i in range(n):
		modeLz[i]=loadtxt(outputdir+'Lz.'+str(i)+'.csv',skiprows=1,usecols=[0],delimiter=',')
	
	for i in range(n):
		Lz[:,i]=modeLz[i]*b[:,i]
	
	savetxt(outputdir+'../Lz.dat',Lz)

outputdir='/home/piccinini/pod_encit/trunk/python/sygma/pods/Lz/'
main(outputdir,10)

outputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/duas_abertas/pods/Lz/'
main(outputdir,10)

outputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/duas_meio_abertas/pods/Lz/'
main(outputdir,10)
