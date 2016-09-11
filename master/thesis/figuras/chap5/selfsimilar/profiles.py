#!/usr/bin/env python

from pylab import *
import numpy as np

D=0.0098

print '(x/D)max = ', 0.5/D
print '(y/D)max = ', 0.075/D

fig_width = 8.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 20,
          'text.fontsize': 20,
          'legend.fontsize': 15,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
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
locs=range(5,50,5) #[5,20,40]
for n in range(len(locs)):
	file='../kEpsilon/sets/0.33/Profiles'+str(locs[n])+'_U.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	u=loadtxt(file,delimiter=' ')[:,1]
	ub=(u-3.)/(u[0]-3.)
	ih=argmin(abs((u-3.)-0.5*(u[0]-3.)))
	rh.append(y[ih])
	print '(x/D,rh,error) =', locs[n],rh[-1],min(abs(u-0.5*u[0]))/u[0]
	fig = figure(1)
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ax.plot(y/rh[-1],ub,'k:'+stle[n],label='x/D = '+str(locs[n]))
	
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	xlabel('$\hat{y}$')
	ylabel('$\hat{U}$')
	xlim([0,2.6])
	ylim([0,1])
	savefig('selfsimilar_U.png')
	savetxt('data_U/U_x'+str(locs[n])+'.dat',vstack([y/rh[-1],ub]).T,delimiter='\t')
	
	file='../kEpsilon/sets/0.252/Profiles'+str(locs[n])+'_k.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	k=loadtxt(file,delimiter=' ')[:,1]
	fig = figure(2)
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ax.plot(y/rh[-1],k/(u[0]-3.)**2,'k:'+stle[n],label='x/D = '+str(locs[n]))
	xlabel('$\hat{y}$')
	ylabel('$\hat{k}$')
	
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	xlim([0,2.6])
	savefig('selfsimilar_k.png')
	savetxt('data_k/k_x'+str(locs[n])+'.dat',vstack([y/rh[-1],k/(u[0]-3.)**2]).T,delimiter='\t')
	
	file='../kEpsilon/sets/0.252/Profiles'+str(locs[n])+'_mut.xy'
	y=loadtxt(file,delimiter=' ')[:,0]
	mut=loadtxt(file,delimiter=' ')[:,1]
	
	file='../kEpsilon/sets/0.252/Profiles'+str(locs[n])+'_rho.xy'
	rho=loadtxt(file,delimiter=' ')[:,1]
	fig = figure(3)
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ax.plot(y/rh[-1],mut/rh[-1]/(u[0]-3.)/rho[0],'k:'+stle[n],label='x/D = '+str(locs[n]),markevery=10)
	xlabel('$\hat{y}$')
	ylabel('\nu_T')
	xlim([0,2.6])
	
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	savefig('selfsimilar_nut.png')
	savetxt('data_nut/nut_x'+str(locs[n])+'.dat',vstack([y/rh[-1],mut/rh[-1]/(u[0]-3.)/rho[0]]).T,delimiter='\t')
	


fig=figure(4)
exp=array([[8.04545,0.594595],[18.9545,1.39189]])
pexp=polyfit(exp[:,0],exp[:,1],1)
Sexp=pexp[0]
x0exp=pexp[1]

pnum=polyfit(locs[:5],array(rh[:5])/D,1)
Snum=pnum[0]
x0num=pnum[1]
print 'Experimental spread rate = ',Sexp
print 'Num. S and x0 = ',Snum

error=(array(rh[:5])/D-polyval(pexp,locs[:5]))/polyval(pexp,locs[:5])

ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
ax.plot(locs[:5],polyval(pnum,locs[:5]),'ok',markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5)
ax.plot(locs[:5],polyval(pexp,locs[:5]),'sr',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5)
xlabel('Axial position - x/D')
ylabel('Half-Radius - $y_{1/2}/D$')
legend(('numerical','experimental'),2).draw_frame(False)
xlim([0,30])
ylim([0,2])
savefig('selfsimilar_spread.png')
print 'spread rate error = ', average(error)
savetxt('spread.dat',vstack([locs[:5],polyval(pnum,locs[:5]),polyval(pexp,locs[:5])]).T,delimiter='\t')

fig=figure(5)
x=25
file='../kEpsilon/sets/0.252/Profiles'+str(x)+'_U.xy'
y=loadtxt(file,delimiter=' ')[:,0]
ye=linspace(0,1.0,500)
u=loadtxt(file,delimiter=' ')[:,1]
ub=(u-3.)/(u[0]-3.)
uexp=np.exp(-log(2.)*(ye/x/D/0.0665)**2)
ih=argmin(abs((u-3.)-0.5*(u[0]-3.)))
rh=y[ih]
print '(x/D,rh,error) =', x,rh,min(abs(u-0.5*u[0]))/u[0]
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
ax.plot(y/x/D,ub,'k-',label='num')
ax.plot(ye/x/D,uexp,'rs',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5,label='exp')
legend(('numerical','experimental'),'upper right').draw_frame(False)
xlabel('$\hat{y}$')
ylabel('$\hat{U}$')
xlim([0,0.25])
ylim([-0.1,1.1])
savefig('selfsimilar_num_exp.png')
savetxt('U_num.dat',vstack([y/x/D,ub]).T,delimiter='\t')
savetxt('U_exp.dat',vstack([ye/x/D,uexp]).T,delimiter='\t')

###### Centerline velocities

file='../kEpsilon/sets/0.252/centerline_U.xy'
yg=loadtxt(file)[:,0]
ug=loadtxt(file)[:,1]
#ug=ug/ug[0]

ye=array([5,10,15,20,25])
ue=[]
for n in [5,10,15,20,25]:
	file='../exp/x'+str(n)+'/d5/Ux.csv'
	ue.append(loadtxt(file,delimiter=',')[0,1])
	
ue=array(ue)

print x0num,x0exp
fig=figure(6)
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
ax.plot(yg/D,ug[0]/ug,'k-',label='num')
ax.plot(ye,ue[0]/ue,'rs',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5,label='exp')
legend()

show()
