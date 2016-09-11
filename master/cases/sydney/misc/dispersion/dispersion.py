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
	
	def smd(self,y,dy):
		self.yc=array(y)
		self.smd=zeros(len(y))
		for j in range(len(y)):
			num=[]
			den=[]
			for i in range(len(self.y)):
				if self.y[i] > y[j]-dy and self.y[i] < y[j]+dy:
					den.append(self.d[i]**2)
					num.append(self.d[i]**3)
			num=array(num)
			den=array(den)
			self.smd[j]=average(num/den)
		
	
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



fig_width = 8.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 18,
          'text.fontsize': 18,
          'legend.fontsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
          'figure.figsize': fig_size,
          'axes.color_cycle' : ['#348ABD', '#7A68A6', '#A60628', '#467821', '#CF4457', '#188487', '#E24A33']
                      # E24A33 : orange
                      # 7A68A6 : purple
                      # 348ABD : blue
                      # 188487 : turquoise
                      # A60628 : red
                      # CF4457 : pink
                      # 467821 : green
           }
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

locs=[5,10,25]
for n in locs:
	ls=[]
	for i in range(len(data)):
		drop=data[i,:]
		if drop[22]/D > n-1. and drop[22]/D < n+1.:
			ls.append([drop[3],drop[13],drop[14],drop[19],drop[20],drop[22]/D,drop[23]/D,drop[9],drop[2]])
	
	ls=array(ls)
	ls=parcel(ls[:,0],ls[:,1],ls[:,2],ls[:,3],ls[:,4],ls[:,5],ls[:,6],ls[:,7],ls[:,8])
	print 'len(ls) = ',len(ls.x)
	
	
	classe=[]
	d0010=[]
	d1020=[]
	d2030=[]
	d3060=[]
	for i in range(len(ls.d)):
		k=[]
		if ls.d[i] <= 10e-6:
			d0010.append(i)
		elif ls.d[i] > 10e-6 and ls.d[i] <= 20e-6:
			d1020.append(i)
		elif ls.d[i] > 20e-6 and ls.d[i] <= 30e-6:
			d2030.append(i)
		elif ls.d[i] > 30e-6:
			d3060.append(i)
	
	classe.append(d0010)
	classe.append(d1020)
	classe.append(d2030)
	classe.append(d3060)
	
	
	
	color=['s','o','^','v']
	
	#fig = figure()
	#ax = fig.add_subplot(1,1,1)
	#subplots_adjust( left=0.15, bottom=0.15 )
	#ax.plot(ls.u[classe[0]],abs(ls.v[classe[0]]),color[0])
	#ax.plot(-1.*ls.u[classe[1]],abs(ls.v[classe[1]]),color[1])
	#ax.plot(-ls.u[classe[2]],-abs(ls.v[classe[2]]),color[2])
	#ax.plot(ls.u[classe[3]],-abs(ls.v[classe[3]]),color[3])
	#xlim([-max((ls.u[:])),max((ls.u[:]))])
	#ylim(-4.,4.)
	
	#title('Axial Position - x/D = '+str(n))
	#xlabel('Axial Velocity (m/s)')
	#ylabel('Radial Velocity (m/s)')
	
	#box = ax.get_position()
	#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

	## Put a legend below current axis
	#ax.legend(('< 10 mu','10-20 mu',' 20-30 mu','> 30 mu'),loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=False, shadow=False, ncol=4).draw_frame(True)
	#savefig('jointUV_'+str(n)+'.png')
	
	fig = figure()
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ymax=max(ls.y[classe[0]])
	yn=linspace(0,ymax,10)
	p=polyfit(ls.y[classe[0]],ls.v[classe[0]],1)
	ax.plot(ls.y[classe[0]],ls.v[classe[0]],'ko',markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5)
	ax.plot(yn, polyval(p,yn),'k-')
	legend(('$D<10\ \mu m$', 'trend line'),'lower right')
	title('Axial coordinate - x/D = '+str(n))
	ylabel('Velocity (m/s)')
	xlabel('y/D')
	ylim([-2,2])
	
	savetxt('class0_x'+str(n)+'.dat',vstack([ls.y[classe[0]],ls.v[classe[0]]]).T,delimiter='\t')
	savefig('class0_x'+str(n)+'.png')
	
	#fig = figure()
	#ax = fig.add_subplot(1,1,1)
	#subplots_adjust( left=0.15, bottom=0.15 )
	#ymax=max(ls.y[classe[1]])
	#yn=linspace(0,ymax,10)
	#p=polyfit(ls.y[classe[1]],ls.v[classe[1]],1)
	#ax.plot(ls.y[classe[1]],ls.v[classe[1]],'ko',markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5)
	#ax.plot(yn, polyval(p,yn),'k-')
	#ylim([-2,2])

	#fig = figure()
	#ax = fig.add_subplot(1,1,1)
	#subplots_adjust( left=0.15, bottom=0.15 )
	#ymax=max(ls.y[classe[2]])
	#yn=linspace(0,ymax,10)
	#p=polyfit(ls.y[classe[2]],ls.v[classe[2]],1)
	#ax.plot(ls.y[classe[2]],ls.v[classe[2]],'ko',markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5)
	#ax.plot(yn, polyval(p,yn),'k-')
	#legend(('$20 \mu m<D<30\ \mu m$', 'trend line'),'lower right')
	#title('Axial coordinate - x/D = '+str(n))
	#ylabel('Velocity (m/s)')
	#xlabel('y/D')
	#ylim([-2,2])

	fig = figure()
	ax = fig.add_subplot(1,1,1)
	subplots_adjust( left=0.15, bottom=0.15 )
	ymax=max(ls.y[classe[3]])
	yn=linspace(0,ymax,10)
	p=polyfit(ls.y[classe[3]],ls.v[classe[3]],1)
	ax.plot(ls.y[classe[3]],ls.v[classe[3]],'ko',markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5)
	ax.plot(yn, polyval(p,yn),'k-')
	legend(('$D>30\ \mu m$', 'trend line'),'lower right')
	title('Axial coordinate - x/D = '+str(n))
	ylabel('Velocity (m/s)')
	xlabel('y/D')
	ylim([-2,2])
	
	savetxt('class4_x'+str(n)+'.dat',vstack([ls.y[classe[3]],ls.v[classe[3]]]).T,delimiter='\t')
	savefig('class4_x'+str(n)+'.png')

show()
