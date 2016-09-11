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
z0,z1,z2,z3=-0.5,0.0,0.2,1.5

nr=['55','33','85','30']
nt='1'
nz=['75','100']

rGrad=['1','1','2.5','0.2']
zGrad=['0.05','7','6']

rad=[]
rad.append(0.010)
rad.append(0.032)
rad.append(0.080)
rad.append(0.100)

vertices=[]
#central point
vertices.append([0,0,0])

#
vertices.append([rad[0]*cos(-wedge/2.),rad[0]*sin(-wedge/2.),-0.5])
vertices.append([rad[0]*cos(wedge/2.),rad[0]*sin(wedge/2.),z0])

vertices.append([rad[1]*cos(-wedge/2.),rad[1]*sin(-wedge/2.),-0.5])
vertices.append([rad[1]*cos(wedge/2.),rad[1]*sin(wedge/2.),z0])

for r in rad:
	vertices.append([r*cos(-wedge/2.),r*sin(-wedge/2.),0])
	vertices.append([r*cos(wedge/2.),r*sin(wedge/2.),z1])

# wall
vertices.append([0,0,0.2])
for r in rad:
	vertices.append([r*cos(-wedge/2.),r*sin(-wedge/2.),0.2])
	vertices.append([r*cos( wedge/2),r*sin( wedge/2),z2])

# outlet
vertices.append([0,0,1.5])
for r in rad:
	vertices.append([r*cos(-wedge/2.),r*sin(-wedge/2.),1.5])
	vertices.append([r*cos( wedge/2.),r*sin( wedge/2.),z3])


for v in vertices:
	f.write('\n\t('+str(v[0])+' '+str(v[1])+' '+str(v[2])+')')


f.write('\n);')
f.write('\n\n blocks (')
#inlet
f.write('\n\t hex (1 3 4 2 5 7 8 6) ('+nr[1]+' '+nt+' '+nz[0]+') simpleGrading ('+rGrad[1] +' 1 '+zGrad[0]+')')
#1
f.write('\n\t hex (0 5 6 0 13 14 15 13) ('+nr[0]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[0] +' 1 '+zGrad[1]+')')
f.write('\n\t hex (5 7 8 6 14 16 17 15) ('+nr[1]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[1] +' 1 '+zGrad[1]+')')
f.write('\n\t hex (13 14 15 13 22 23 24 22) ('+nr[0]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[0] +' 1 '+zGrad[2]+')')
f.write('\n\t hex (14 16 17 15 23 25 26 24) ('+nr[1]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[1] +' 1 '+zGrad[2]+')')
f.write('\n\t hex (7 9 10 8 16 18 19 17) ('+nr[2]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[2] +' 1 '+zGrad[1]+')')
f.write('\n\t hex (16 18 19 17 25 27 28 26) ('+nr[2]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[2] +' 1 '+zGrad[2]+')')
f.write('\n\t hex (9 11 12 10 18 20 21 19) ('+nr[3]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[3] +' 1 '+zGrad[1]+')')
f.write('\n\t hex (18 20 21 19 27 29 30 28) ('+nr[3]+' '+nt+' '+nz[1]+') simpleGrading ('+rGrad[3] +' 1 '+zGrad[2]+')')

f.write('\n);')

f.write('\n\nedges (')

v=array([5,14,23])
z=array([z1,z2,z3])

f.write('\n\t arc 1 2 ('+str(rad[0])+' 0 '+str(z0)+')')
f.write('\n\t arc 3 4 ('+str(rad[1])+' 0 '+str(z0)+')')

for i in range(len(v)):
	f.write('\n\t arc '+str(v[i])+' '+str(v[i]+1)+' ('+str(rad[0])+' 0 '+str(z[i])+')')
	f.write('\n\t arc '+str(v[i]+2)+' '+str(v[i]+3)+' ('+str(rad[1])+' 0 '+str(z[i])+')')
	f.write('\n\t arc '+str(v[i]+4)+' '+str(v[i]+5)+' ('+str(rad[2])+' 0 '+str(z[i])+')')
	f.write('\n\t arc '+str(v[i]+6)+' '+str(v[i]+7)+' ('+str(rad[3])+' 0 '+str(z[i])+')')

j=1
for i in range(4,16,4):
	r=rad[j]
	#~ f.write('\n\t arc '+str(i+16)+' '+str(i+17)+' (0 '+str(-r*k)+' '+str(z2)+')')
	j=j+1

f.write('\n);')

# patches
f.write('\n\npatches\n(')
# inlet
f.write('\n\t patch inlet (')
f.write('\n\t\t (1 2 4 3)')
f.write('\n\t )')
# walls
f.write('\n\t wall walls (')
f.write('\n\t\t (1 5 6 2)')
f.write('\n\t\t (3 7 8 4)')
f.write('\n\t\t (0 6 5 0)')
f.write('\n\t\t (7 8 10 9)')
f.write('\n\t\t (9 10 12 11)')
f.write('\n\t\t (20 11 12 21)')
f.write('\n\t\t (29 20 21 30)')
f.write('\n\t )')

# wedge back
f.write('\n\t wedge back (')
f.write('\n\t\t (1 3 7 5)')
f.write('\n\t\t (0 5 14 13)')
f.write('\n\t\t (5 7 16 14)')
f.write('\n\t\t (7 9 18 16)')
f.write('\n\t\t (9 11 20 18)')
f.write('\n\t\t (13 14 23 22)')
f.write('\n\t\t (14 16 25 23)')
f.write('\n\t\t (16 18 27 25)')
f.write('\n\t\t (18 20 29 27)')
f.write('\n\t )')

# wedge front
f.write('\n\t wedge front (')
f.write('\n\t\t (2 6 8 4)')
f.write('\n\t\t (0 13 15 6)')
f.write('\n\t\t (6 15 17 8)')
f.write('\n\t\t (8 17 19 10)')
f.write('\n\t\t (10 19 21 12)')
f.write('\n\t\t (13 22 24 15)')
f.write('\n\t\t (15 24 26 17)')
f.write('\n\t\t (17 26 28 19)')
f.write('\n\t\t (19 28 30 21)')
f.write('\n\t )')
# outlet
f.write('\n\t patch outlet (')
f.write('\n\t\t (22 23 24 22)')
f.write('\n\t\t (23 25 26 24)')
f.write('\n\t\t (25 27 28 26)')
f.write('\n\t\t (27 29 30 28)')
f.write('\n\t )')
f.write('\n);')
# --

f.write('\n\n mergePatchPairs\n(\n);')

f.write('\n\n// ************************************************************************* // ')
f.close()
#os.system('xterm -hold -e \'source ~/.cshrc ;  blockMesh\'')












