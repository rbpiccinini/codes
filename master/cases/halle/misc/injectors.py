#!/usr/bin/env python

from pylab import *

temperature=273. +35.
massflow=array([1.6483E-05, 8.0020E-05, 1.8868E-04, 2.6372E-04, 1.6256E-04, 7.3488E-05, 2.4049E-05, 8.4273E-06, 3.1542E-06, 1.4191E-06])
print sum(massflow)
t0=0.
tf=5.

diameter=0.001

x=array([4,5,6,7,8,9,10,11,12,13])*1e-3
ninj=len(x)
vel=array([12.56, 15.66, 17.39, 15.72, 12.75, 10.79, 9.89, 9.49, 9.29, 9.16])

interval=tf-t0

oneDrop=1000.*4./3.*pi*(30e-6)**3
nParcels=1 #int(mass/oneDrop)

f=open('injectorProperties','w')
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
f.write('    root            "";\n')
f.write('    case            "";\n')
f.write('    instance        "";\n')
f.write('    local           "";\n')
f.write('\n')
f.write('    class           dictionary; \n')
f.write('    object          injectorProperties;\n')
f.write('}\n')
f.write('\n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')

f.write('(')

for i in range(ninj):
	f.write('\n\n// z = '+str(x[i])+' ------------------------------------------------------------------------- \n')
	f.write('\n    {')
	f.write('\n        injectorType        definedInjector;')
	f.write('\n        definedInjectorProps')
	f.write('\n        {')
	f.write('\n            position        ('+str(x[i])+' 0 1e-4);')
	f.write('\n            direction       (0 0 1);')
	f.write('\n            diameter        '+str(diameter)+';')
	f.write('\n            Cd              1.0;')
	f.write('\n            mass            '+str(massflow[i])+';')
	f.write('\n            nParcels        '+str(nParcels)+';')
	f.write('\n            temperature     '+str(temperature)+';')
	f.write('\n            X ( 1.0 );')
	f.write('\n            ')
	f.write('\n            massFlowRateProfile ')
	f.write('\n            (')
	f.write('\n                ('+str(t0)+'     1)')
	f.write('\n                ('+str(tf)+'     1)')
	f.write('\n           );')
	f.write('\n            temperatureProfile ')
	f.write('\n            (')
	f.write('\n                ('+str(t0)+'     '+str(temperature)+')')
	f.write('\n                ('+str(tf)+'     '+str(temperature)+')')
	f.write('\n           );')
	f.write('\n            velocityProfile ')
	f.write('\n            (')
	f.write('\n                ('+str(t0)+'     '+str(vel[i])+')')
	f.write('\n                ('+str(tf)+'     '+str(vel[i])+')')
	f.write('\n           );')
	f.write('\n        }')
	f.write('\n    }')

f.write('\n// ------------------------------------------------------------------------- \n')
f.write(')')
f.close()
