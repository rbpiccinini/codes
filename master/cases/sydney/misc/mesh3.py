#!/usr/bin/env python
from pylab import *
import os

f=open('blockMeshDict','w')
f.write('/*---------------------------------------------------------------------------*\ \n')
f.write('| =========                 |                                                 | \n')
f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
f.write('|  \\    /   O peration     | Version:  1.4                                   | \n')
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
f.write('    object          blockMeshDict;\n')
f.write('}\n')
f.write('\n')
f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')

f.write('\n')
f.write('convertToMeters 1;\n')
f.write('\n')
f.write('vertices (')

wedge=pi/180.*5.
x0=0.5

# fine mesh
nx=200
ny=50
nyNozzle=24
nyWall=15
xGrad=10
yGrad=15
yGradWall=0.1

r1=0.0049
r2=0.065
r3=0.075

vertices=[]
#central point
vertices.append([0,0,0])#0
vertices.append([0,r1*cos(wedge/2.),r1*sin(-wedge/2.)])#1
vertices.append([x0,r1*cos(wedge/2.),r1*sin(-wedge/2.)])#2
vertices.append([x0,0,0])#3
vertices.append([0,r1*cos(wedge/2.),r1*sin(wedge/2.)])#4
vertices.append([x0,r1*cos(wedge/2.),r1*sin(wedge/2.)])#5

vertices.append([0,r2*cos(wedge/2.),r2*sin(-wedge/2.)])#6
vertices.append([x0,r2*cos(wedge/2.),r2*sin(-wedge/2.)])#7
vertices.append([0,r2*cos(wedge/2.),r2*sin(wedge/2.)])#8
vertices.append([x0,r2*cos(wedge/2.),r2*sin(wedge/2.)])#9

vertices.append([0,r3*cos(wedge/2.),r3*sin(-wedge/2.)])#10
vertices.append([x0,r3*cos(wedge/2.),r3*sin(-wedge/2.)])#11
vertices.append([0,r3*cos(wedge/2.),r3*sin(wedge/2.)])#12
vertices.append([x0,r3*cos(wedge/2.),r3*sin(wedge/2.)])#13

for v in vertices:
	f.write('\n\t('+str(v[0])+' '+str(v[1])+' '+str(v[2])+')')


f.write('\n);')
f.write('\n\n blocks (')
#inlet
f.write('\n\t hex (0 3 2 1 0 3 5 4) ('+str(nx)+' '+str(nyNozzle)+' 1) simpleGrading ('+str(xGrad)+' 1 1)')
f.write('\n\t hex (1 2 7 6 4 5 9 8) ('+str(nx)+' '+str(ny)+' 1) simpleGrading ('+str(xGrad)+' '+str(yGrad)+' 1)')
f.write('\n\t hex (6 7 11 10 8 9 13 12) ('+str(nx)+' '+str(nyWall)+' 1) simpleGrading ('+str(xGrad)+' '+str(yGradWall)+' 1)')

f.write('\n);')

f.write('\n\nedges (')

f.write('\n\t arc 1 4 (0 '+str(r1)+' 0)')
f.write('\n\t arc 2 5 ('+str(x0)+' '+str(r1)+' 0)')
f.write('\n\t arc 6 8 (0 '+str(r2)+' 0)')
f.write('\n\t arc 7 9 ('+str(x0)+' '+str(r2)+' 0)')
f.write('\n\t arc 10 12 (0 '+str(r3)+' 0)')
f.write('\n\t arc 11 13 ('+str(x0)+' '+str(r3)+' 0)')

f.write('\n);')

# patches
f.write('\n\npatches\n(')
# inlet
f.write('\n\t patch nozzle (')
f.write('\n\t\t (0 4 1 0)')
f.write('\n\t )')
# coflow
f.write('\n\t patch coflow (')
f.write('\n\t\t (1 4 8 6)')
f.write('\n\t\t (6 8 12 10)')
f.write('\n\t )')
# walls
f.write('\n\t wall walls (')
f.write('\n\t\t (10 12 13 11)')
f.write('\n\t )')

# wedge front and back
f.write('\n\t wedge back (')
f.write('\n\t\t (0 1 2 3)')
f.write('\n\t\t (1 6 7 2)')
f.write('\n\t\t (6 10 11 7)')
f.write('\n\t )')
f.write('\n\t wedge front (')
f.write('\n\t\t (0 3 5 4)')
f.write('\n\t\t (4 5 9 8)')
f.write('\n\t\t (8 9 13 12)')
f.write('\n\t )')
# outlet
f.write('\n\t patch outlet (')
f.write('\n\t\t (3 2 5 3)')
f.write('\n\t\t (5 2 7 9)')
f.write('\n\t\t (9 7 11 13)')
f.write('\n\t )')
f.write('\n);')
# --

f.write('\n\n mergePatchPairs\n(\n);')

f.write('\n\n// ************************************************************************* // ')
f.close()
#os.system('xterm -hold -e \'source ~/.cshrc ;  blockMesh\'')












