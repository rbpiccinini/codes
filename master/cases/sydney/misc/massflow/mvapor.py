#!/usr/bin/env python

from pylab import *

def adjust_spines(ax,spines):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


fig_width = 6.  # Get this from LaTeX using \showthe\columnwidth
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

# 0 "origId"
# 1 "origProcId"
# 2 "m"
# 3 "d"
# 4 "yDot"
# 5 "injector"
# 6 "ct"
# 7 "ms"
# 8 "aC3H6O"
# 9 "T"
# 10 "tTurb"
# 11 "liquidCore"
# 12 "y"
# 13 "U:0"
# 14 "U:1"
# 15 "U:2"
# 16 "n:0"
# 17 "n:1"
# 18 "n:2"
# 19 "Uturb:0"
# 20 "Uturb:1"
# 21 "Uturb:2"
# 22 "Points:0"
# 23 "Points:1"
# 24 "Points:2"

D=0.0098
S=0.25*pi*D**2
print '(x/D)max = ', 0.5/D
print '(y/D)max = ', 0.075/D

xl=[1,2,3,3]
nticks=[5,5,5,5,5]

locations=[5,10,15,20]

minlet=5./S/1e4

for n in range(len(locations)):
	file='../kEpsilon/mflux/coarse/mvapor_'+str(locations[n])+'.csv'
	yn=loadtxt(file,delimiter=',',skiprows=1)[:-1,-2]
	mn=loadtxt(file,delimiter=',',skiprows=1)[:-1,-7]
	
	file='../exp/x'+str(locations[n])+'/mvapor.csv'
	ye=loadtxt(file,delimiter=',',skiprows=1)[:,0]
	me=loadtxt(file,delimiter=',',skiprows=1)[:,1]
	
	minu=trapz(mn*yn,x=yn)*2*pi*1e3*60.
	mexp=(0.5*me[0]*(ye[0]*D)**2+trapz(me*ye*D,x=ye*D))*2*pi*1e4
	print 'vapour mass flow [inlet, num, exp] = ', minlet, minu, mexp
	print 'ratio = ', mexp/minu
	# convert to exp units
	yn=yn/D
	mn=mn*6.
	
	#print min(ls.x),max(ls.x)
	#figure(1)
	#plot(n,sum(ls.m[:]),'bo')
	#title('Liquid Mass [kg]')
	#xlabel('x/D')
	
	#figure()
	#plot(ls.yc,ls.uc,'-b',linewidth=2)
	#plot(ls.y,ls.u,'.k')
	#plot(ye,ue,'o-r')
	#xlim([0,max(ye)])
	#ylim([0,max(ue)])
	#xlabel('y/D')
	#ylabel('Axial Velocity - U (m/s)')
	#title('Axial Position - x/D = '+str(n))
	
	
	fig = figure()
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
	ax.plot(yn,mn,'k-')
	ax.plot(ye,me,'sr',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5)
	
	
	legend(('numerical','experimental')).draw_frame(False)
	xlim([0,2])
	ylim([0,1.1*max(mn)])
	xlabel('y/D')
	ylabel('Vapor Mass Flux $(g/cm2/min)$')
	title('Axial Position - x/D = '+str(locations[n]))
	
	savefig('mvapor'+str(locations[n])+'.png')
	savetxt('mvapor_x'+str(locations[n])+'_num.dat',vstack([yn,mn]).T,delimiter='\t')
	savetxt('mvapor_x'+str(locations[n])+'_exp.dat',vstack([ye,me]).T,delimiter='\t')


file='../Sevap_vol.csv'
y=loadtxt(file,delimiter=',',skiprows=1)[:,-2]
Sevap=loadtxt(file,delimiter=',',skiprows=1)[:,11]

figure()
plot(y/D,Sevap/max(Sevap))
xlabel('y/D')
show()
