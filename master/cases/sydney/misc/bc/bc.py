#!/usr/bin/env python

from pylab import *

def points(y,theta):
	f=open('points','w')
	f.write('/*---------------------------------------------------------------------------*\ \n')
	f.write('| =========                 |                                                 | \n')
	f.write('| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
	f.write('|  \    /   O peration     | Version:  1.6                                   | \n')
	f.write('|   \  /    A nd           | Web:      http://www.openfoam.org               | \n')
	f.write('|    \/     M anipulation  |                                                 | \n')
	f.write('\*---------------------------------------------------------------------------*/ \n')
	f.write('\n')
	f.write('FoamFile \n')
	f.write('{    version         2.0; \n')
	f.write('	format          ascii; \n')
	f.write('\n')
	f.write('    class           vectorField; \n')
	f.write('    object          points;\n')
	f.write('}\n')
	f.write(' ')
	f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
	
	
	f.write('(')
	
	x=0.
	z=0.
	for i in y:
		f.write('\n\t\t ('+str(x)+' '+str(i*cos(theta))+' '+str(-1.)+')')
	for i in y[1:]:
		f.write('\n\t\t ('+str(x)+' '+str(i*cos(theta))+' '+str(1.)+')')
	
	f.write('\n)')	
	f.close()
	
def velocity(Ux):
	f=open('U','w')
	f.write('/*---------------------------------------------------------------------------*\ \n')
	f.write('| =========                 |                                                 | \n')
	f.write('| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
	f.write('|  \    /   O peration     | Version:  1.6                                   | \n')
	f.write('|   \  /    A nd           | Web:      http://www.openfoam.org               | \n')
	f.write('|    \/     M anipulation  |                                                 | \n')
	f.write('\*---------------------------------------------------------------------------*/ \n')
	f.write('\n')
	f.write('FoamFile \n')
	f.write('{    version         2.0; \n')
	f.write('	format          ascii; \n')
	f.write('\n')
	f.write('    class           vectorAverageField; \n')
	f.write('    object          values;\n')
	f.write('}\n')
	f.write(' ')
	f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
	f.write('\n\n('+str(mean(Ux))+' 0 0)')
	f.write('\n'+str(len(Ux)*2-1)+'\n')
	f.write('(')
	Uy=0.
	Uz=0.
	for i in Ux:
		f.write('\n\t\t ('+str(i)+' '+str(Uy)+' '+str(Uz)+')')
	for i in Ux[1:]:
		f.write('\n\t\t ('+str(i)+' '+str(Uy)+' '+str(Uz)+')')
	
	f.write('\n)')	
	f.close()
	
def turbulentK(k):
	f=open('k','w')
	f.write('/*---------------------------------------------------------------------------*\ \n')
	f.write('| =========                 |                                                 | \n')
	f.write('| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
	f.write('|  \    /   O peration     | Version:  1.6                                   | \n')
	f.write('|   \  /    A nd           | Web:      http://www.openfoam.org               | \n')
	f.write('|    \/     M anipulation  |                                                 | \n')
	f.write('\*---------------------------------------------------------------------------*/ \n')
	f.write('\n')
	f.write('FoamFile \n')
	f.write('{    version         2.0; \n')
	f.write('	format          ascii; \n')
	f.write('\n')
	f.write('    class           scalarAverageField; \n')
	f.write('    object          values;\n')
	f.write('}\n')
	f.write(' ')
	f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
	f.write('\n\n'+str(mean(k)))
	f.write('\n'+str(len(Ux)*2-1)+'\n')
	f.write('(')
	Uy=0.
	Uz=0.
	for i in k:
		f.write('\n\t\t '+str(i))
	for i in k[1:]:
		f.write('\n\t\t '+str(i))
	
	f.write('\n)')	
	f.close()
	
def epsilon(eps):
	f=open('epsilon','w')
	f.write('/*---------------------------------------------------------------------------*\ \n')
	f.write('| =========                 |                                                 | \n')
	f.write('| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
	f.write('|  \    /   O peration     | Version:  1.6                                   | \n')
	f.write('|   \  /    A nd           | Web:      http://www.openfoam.org               | \n')
	f.write('|    \/     M anipulation  |                                                 | \n')
	f.write('\*---------------------------------------------------------------------------*/ \n')
	f.write('\n')
	f.write('FoamFile \n')
	f.write('{    version         2.0; \n')
	f.write('	format          ascii; \n')
	f.write('\n')
	f.write('    class           scalarAverageField; \n')
	f.write('    object          values;\n')
	f.write('}\n')
	f.write(' ')
	f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n')
	f.write('\n\n'+str(mean(eps)))
	f.write('\n'+str(len(eps)*2-1)+'\n')
	f.write('(')
	
	for i in eps:
		f.write('\n\t\t '+str(i))
	for i in eps[1:]:
		f.write('\n\t\t '+str(i))
	
	f.write('\n)')	
	f.close()

#flow field
dt=loadtxt('Ubc.csv',delimiter=',',skiprows=1)
y=dt[:,8]
Ux=dt[:,4]/max(dt[:,4])*31.
k=dt[:,1]
eps=dt[:,3]

D=0.0098

#droplets
yp=array([0.0431603,0.11217,0.175433,0.242521,0.305758,0.370879,0.43393,0.498768,0.563718])*.0098
Ud=array([30.9176,30.9176,31.0588,30.9176,30.2118,28.5176,22.0235,11.4353,4.37647])/30.
p=polyfit(yp[-4:],Ud[-4:],1)
print p

#plot(yp,Ud,'or')
#plot(yp,polyval(p,yp),'-b')

file='../exp/bc_x_0p5/up.csv'
ye=loadtxt(file,delimiter=',',skiprows=1)[:,0]*D
ue=loadtxt(file,delimiter=',',skiprows=1)[:,1]
file='../exp/bc_x_0p5/vp.csv'
ve=loadtxt(file,delimiter=',',skiprows=1)[:,1]
ke=0.5*(ue**2+2.*ve**2)

# setting plot style

fig_width = 8.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 13,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
          'figure.figsize': fig_size}
rcParams.update(params)

D=0.0098

# plot data
figure()
plot(y/D,Ux,'k-')
plot(yp/D,Ud*30.,':sr')
ylabel('Axial velocity (m/s)')
xlabel('Radial position (y/D)')
xlim([0,0.5])
legend(('simulated BC','exp'),'lower right')
savetxt('bc_Ux.dat',vstack([yp/D,Ud*30.]).T,delimiter='\t')
savefig('bc_Ux.png')

figure()
plot(y/D,k,'k-')
xlabel('Radial Position (y/D)')
ylabel('Turbulent KE (m2/s2)')
savetxt('bc_k.dat',vstack([y/D,k]).T,delimiter='\t')
savefig('bc_k.png')

figure()
plot(y/D,eps*1e-4,'k-')
xlabel('Radial Position (y/D)')
ylabel('Dissipation Rate (1e4 m2/s3)')
savetxt('bc_eps.dat',vstack([y/D,eps*1e-4]).T,delimiter='\t')
savefig('bc_eps.png')
show()

if len(y) == len(Ux):
	print 'Dimensions ok.'
	points(y,5./180.*pi)
	velocity(Ux)
	turbulentK(k)
	epsilon(eps)
