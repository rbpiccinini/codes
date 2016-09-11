#!/usr/bin/env python

from pylab import *
import numpy as np
import xlwt

D=0.0098

print '(x/D)max = ', 0.5/D
print '(y/D)max = ', 0.075/D

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
          'axes.color_cycle' : ['#348ABD', '#7A68A6', '#A60628', '#467821', '#CF4457', '#188487', '#E24A33'],
          'figure.figsize': fig_size}
rcParams.update(params)

#'o' 	circle
#'v' 	triangle_down
#'^' 	triangle_up
#'<' 	triangle_left
#'>' 	triangle_right
#'1' 	tri_down
#'2' 	tri_up
#'3' 	tri_left
#'4' 	tri_right
#'s' 	square
#'p' 	pentagon
#'*' 	star
#'h' 	hexagon1
#'H' 	hexagon2
#'+' 	plus


stle=['o','x','v','^','>','s','p','*','+']
rh=[]
locs=range(5,30,5) #[5,20,40]
for n in range(len(locs)):
	
	file='data/fine/x'+str(locs[n])+'_U.xy'
	yf=loadtxt(file,delimiter=' ')[:,0]
	uf=loadtxt(file,delimiter=' ')[:,1]
	
	file='data/coarse/x'+str(locs[n])+'_U.xy'
	yc=loadtxt(file,delimiter=' ')[:,0]
	uc=loadtxt(file,delimiter=' ')[:,1]
	
	fig = figure()
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ax.plot(yf/D,uf,'k-',label='fine')
	ax.plot(yc/D,uc,'k+',label='coarse',markeredgewidth=1.5)
	#ax.plot(ys/D,us,'^',label='spray')
	
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	xlabel('y / D')
	ylabel('Axial Velocity (m/s)')
	title('Axial Position - x/D = '+str(locs[n]))
	xlim([0,4])
	#ylim([0,1])
	savefig('grid_'+str(n)+'_U.png')
	
	
	file='data/fine/x'+str(locs[n])+'_k.xy'
	yf=loadtxt(file,delimiter=' ')[:,0]
	kf=loadtxt(file,delimiter=' ')[:,1]
	
	file='data/coarse/x'+str(locs[n])+'_k.xy'
	yc=loadtxt(file,delimiter=' ')[:,0]
	kc=loadtxt(file,delimiter=' ')[:,1]
	
	fig = figure()
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ax.plot(yf/D,kf,'k-',label='fine')
	ax.plot(yc/D,kc,'k+',label='coarse',markeredgewidth=1.5)
	
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	xlabel('y / D')
	ylabel('Turbulent KE (J/kg)')
	title('Axial Position - x/D = '+str(locs[n]))
	xlim([0,4])
	#ylim([0,1])

	savefig('grid_'+str(n)+'_k.png')
show()

