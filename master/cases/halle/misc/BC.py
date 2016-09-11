#!/usr/bin/env python

from pylab import *

path='../3D/rans/'
nTheta=36
z=1.75
x=load('channel.csv',delimiter=',',skiprows=1)
r=x[::5,0]
k=x[::5,2]
eps=x[::5,2]
U=x[::5,3]
th=linspace(0,2*pi,nTheta+1)
#th=th[:-1]
points=[]

w=0.032/(pi*(0.032**2-0.02**2))
meps=mean(eps)
mk=mean(k)
#-----------------------------------------------------------------------
# POINTS
#-----------------------------------------------------------------------
f=open(path+'constant/boundaryData/inlet/points','w')
f.write('/*---------------------------------------------------------------------------*\ \n')
f.write('| =========                 |                                                 | \n')
f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
f.write('|  \\    /   O peration     | Version:  1.6                                   | \n')
f.write('|   \\  /    A nd           | Web:      http://www.openfoam.org               | \n')
f.write('|    \\/     M anipulation  |                                                 | \n')
f.write('\*---------------------------------------------------------------------------*/ \n')
f.write('\n')
f.write('FoamFile \n')
f.write('{')
f.write('    version         2.0; \n')
f.write('    format          ascii; \n')
f.write('\n')
f.write('    class           vectorField; \n')
f.write('    object          points;\n')
f.write('}\n')
f.write('\n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
f.write('\n')
f.write('\n(')
for i in range(nTheta):
	for j in range(len(r)):
		f.write('\n ('+str(r[j]*cos(th[i]))+' '+str(r[j]*sin(th[i])) +' '+str(z)+')')
		points.append([r[j]*cos(th[i]),r[j]*sin(th[i])])
f.write('\n)')

#-----------------------------------------------------------------------
# U
#-----------------------------------------------------------------------
f=open(path+'constant/boundaryData/inlet/0/U','w')
f.write('/*---------------------------------------------------------------------------*\ \n')
f.write('| =========                 |                                                 | \n')
f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
f.write('|  \\    /   O peration     | Version:  1.6                                   | \n')
f.write('|   \\  /    A nd           | Web:      http://www.openfoam.org               | \n')
f.write('|    \\/     M anipulation  |                                                 | \n')
f.write('\*---------------------------------------------------------------------------*/ \n')
f.write('\n')
f.write('FoamFile \n')
f.write('{')
f.write('    version         2.0; \n')
f.write('    format          ascii; \n')
f.write('\n')
f.write('    class           vectorAverageField; \n')
f.write('    object          values;\n')
f.write('}\n')
f.write('\n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
f.write('\n')
f.write('(0 0 '+str(-w)+')')
f.write('\n')
f.write(str(nTheta*len(r)))
f.write('\n(')
for i in range(nTheta):
	for j in range(len(r)):
		f.write('\n (0 0 '+' '+str(-U[j])+')')
f.write('\n)')

#-----------------------------------------------------------------------
# K
#-----------------------------------------------------------------------
f=open(path+'constant/boundaryData/inlet/0/k','w')
f.write('/*---------------------------------------------------------------------------*\ \n')
f.write('| =========                 |                                                 | \n')
f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
f.write('|  \\    /   O peration     | Version:  1.6                                   | \n')
f.write('|   \\  /    A nd           | Web:      http://www.openfoam.org               | \n')
f.write('|    \\/     M anipulation  |                                                 | \n')
f.write('\*---------------------------------------------------------------------------*/ \n')
f.write('\n')
f.write('FoamFile \n')
f.write('{')
f.write('    version         2.0; \n')
f.write('    format          ascii; \n')
f.write('\n')
f.write('    class           scalarAverageField; \n')
f.write('    object          values;\n')
f.write('}\n')
f.write('\n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
f.write('\n')
f.write('\n'+str(mk))
f.write('\n')
f.write(str(nTheta*len(r)))
f.write('\n(')
for i in range(nTheta):
	for j in range(len(r)):
		f.write('\n'+str(k[j]))
f.write('\n)')

#-----------------------------------------------------------------------
# Epsilon
#-----------------------------------------------------------------------
f=open(path+'constant/boundaryData/inlet/0/epsilon','w')
f.write('/*---------------------------------------------------------------------------*\ \n')
f.write('| =========                 |                                                 | \n')
f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
f.write('|  \\    /   O peration     | Version:  1.6                                   | \n')
f.write('|   \\  /    A nd           | Web:      http://www.openfoam.org               | \n')
f.write('|    \\/     M anipulation  |                                                 | \n')
f.write('\*---------------------------------------------------------------------------*/ \n')
f.write('\n')
f.write('FoamFile \n')
f.write('{')
f.write('    version         2.0; \n')
f.write('    format          ascii; \n')
f.write('\n')
f.write('    class           scalarAverageField; \n')
f.write('    object          values;\n')
f.write('}\n')
f.write('\n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
f.write('\n')
f.write('\n'+str(meps))
f.write('\n')
f.write(str(nTheta*len(r)))
f.write('\n(')
for i in range(nTheta):
	for j in range(len(r)):
		f.write('\n'+str(eps[j]))
f.write('\n)')

points=array(points)
plot(points[:,0],points[:,1],'o')
#show()
