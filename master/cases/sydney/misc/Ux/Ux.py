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

def reyn(x,dx,xp,yp):
	w=[]
	for i in range(len(xp)):
		if xp[i] > x-dx and xp[i] < x+dx:
			w.append(yp[i])
	w=array(w)
	if len(w) == 1:
		return 0.
	return average(w)

class parcel:
	def __init__(self,d,u,v,uu,vv,x,y,T,m):
		self.d=d
		self.u=u
		self.v=v
		self.uu=uu
		self.vv=vv
		self.x=x
		self.y=y
		self.T=T
		self.m=m
	
	def reynolds(self,y,dy):
		self.yc=array(y)
		self.uc=zeros(len(y))
		self.uuc=zeros(len(y))
		self.vc=zeros(len(y))
		self.vvc=zeros(len(y))
		
		for j in range(len(y)):
			ut=[]
			uut=[]
			vt=[]
			vvt=[]
			mt=[]
			for i in range(len(self.y)):
				if self.y[i] > y[j]-dy and self.y[i] < y[j]+dy:
					#mt.append(self.m[i]/(pi/6.*self.d[i]**3))
					mt.append(self.m[i])
					ut.append(self.u[i])
					uut.append(self.uu[i])
					vt.append(self.v[i])
					vvt.append(self.vv[i])
			ut=array(ut)
			uut=array(uut)
			vt=array(vt)
			vvt=array(vvt)
			if not mt:
				break
			self.uc[j]=average(ut,weights=mt)
			self.uuc[j]=sqrt(average(uut**2,weights=mt))
			self.vc[j]=average(vt,weights=mt)
			self.vvc[j]=sqrt(average(vvt**2,weights=mt)) 


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

data=loadtxt('../drops.csv',skiprows=1,delimiter=',')
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
print '(x/D)max = ', 0.5/D
print '(y/D)max = ', 0.075/D

centerlineU=[]
for n in [5,10,15,20,25]:
	
	file='../exp/x'+str(n)+'/d5/Ux.csv'
	ye=loadtxt(file,delimiter=',')[:,0]
	ue=loadtxt(file,delimiter=',')[:,1]
	
	centerlineU.append([n,ue[0]])
	
	file='../exp/x'+str(n)+'/d15/Ux.csv'
	yed=loadtxt(file,delimiter=',')[:,0]
	ued=loadtxt(file,delimiter=',')[:,1]
	
	file='../kEpsilon/sets/0.252/Profiles'+str(n)+'_U.xy'
	yg=loadtxt(file)[:,0]
	ug=loadtxt(file)[:,1]
	
	#ye=array([0.0431603,0.11217,0.175433,0.242521,0.305758,0.370879,0.43393,0.498768,0.563718])
	#ue=array([30.9176,30.9176,31.0588,30.9176,30.2118,28.5176,22.0235,11.4353,4.37647])
	
	ls=[]
	for i in range(len(data)):
		drop=data[i,:]
		if drop[22]/D > n-.25 and drop[22]/D < n+.25:
			ls.append([drop[3],drop[13],drop[14],drop[19],drop[20],drop[22]/D,drop[23]/D,drop[9],drop[2]])
	
	ls=array(ls)
	ls=parcel(ls[:,0],ls[:,1],ls[:,2],ls[:,3],ls[:,4],ls[:,5],ls[:,6],ls[:,7],ls[:,8])
	
	#yc=linspace(0.,max(ls.y),30)
	ls.reynolds(ye,0.2)
	
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
	
### Figure 1	
	fig = figure()
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
	ax.plot(ls.y,ls.u,'k.')
	ax.plot(yed,ued,'sr',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5)

	print average(ls.u)/average(ue)
	
	legend(('numerical','experimental')).draw_frame(False)
	ylim(0,max(ued)+2.)
	xlim(0,max(ye))
	
	xlabel('y/D')
	ylabel('Axial Velocity - U_x (m/s)')
	title('Axial Position - x/D = '+str(n))
	
	savetxt('Ux_drops/Ux_x'+str(n)+'_num.dat',vstack([ls.y,ls.u]).T,delimiter='\t')
	savetxt('Ux_drops/Ux_x'+str(n)+'_exp.dat',vstack([yed,ued]).T,delimiter='\t')
	savefig('Ux_drops/Ux_x'+str(n)+'.png')


### Figure 2
	fig = figure()
	# plot data
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
	ax.plot(yg/D,ug,'k-')
	ax.plot(ye,ue,'sr',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5)
	
	print average(ls.u)/average(ue)
	
	legend(('numerical','experimental')).draw_frame(False)
	ylim(0,max(ued)+2.)
	xlim(0,max(ye))
	
	xlabel('y/D')
	ylabel('Axial Velocity - U_x (m/s)')
	title('Axial Position - x/D = '+str(n))
	
	savetxt('Ux_gas/Ux_gas_x'+str(n)+'_num.dat',vstack([yg/D,ug]).T,delimiter='\t')
	savetxt('Ux_gas/Ux_gas_x'+str(n)+'_exp.dat',vstack([ye,ue]).T,delimiter='\t')
	savefig('Ux_gas/Ux_gas'+str(n)+'.png')
	

### Figure 3
centerlineU=array(centerlineU)

file='../kEpsilon/sets/0.252/centerline_U.xy'
xn=loadtxt(file)[:,0]
un=loadtxt(file)[:,1]

fig_width = 7.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
rcParams.update({'figure.figsize': fig_size})

fig = figure()
# plot data
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
ax.plot(xn/D,un/max(un),'k-')
ax.plot(centerlineU[:,0],centerlineU[:,1]/max(un),'sr',markerfacecolor='None',markeredgecolor='r',markeredgewidth=1.5)
ylim([0,1.1])
xlim([0,50])

print average(ls.u)/average(ue)

legend(('numerical','experimental')).draw_frame(False)

xlabel('x/D')
ylabel('Centerline U_x (m/s)')

savetxt('centerline/centerline_U_num.dat',vstack([xn/D,un/max(un)]).T,delimiter='\t')
savetxt('centerline/centerline_U_exp.dat',vstack([centerlineU[:,0],centerlineU[:,1]/max(un)]).T,delimiter='\t')
savefig('centerline/centerline_U.png')

show()
