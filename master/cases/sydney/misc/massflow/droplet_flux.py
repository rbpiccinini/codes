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


fig_width = 7.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 15,
          'text.fontsize': 15,
          'legend.fontsize': 13,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
          'figure.figsize': fig_size}
rcParams.update(params)

times=10
data=[]
for i in range(times):
	#data.append(loadtxt('drops.csv',skiprows=1,delimiter=','))
	data.append(loadtxt('../drops_average/drops.'+str(i)+'.csv',skiprows=1,delimiter=','))
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

exp=loadtxt('../exp/droplet_flux.csv',delimiter=',')
locs=array([1,5,10,15,20,25,30])
G=3.33e-5
mflux=[]
for n in locs:
	#file='exp/x'+str(n)+'/d5/UUx.csv'
	#ye=loadtxt(file,delimiter=',',skiprows=1)[:,0]
	#ue=loadtxt(file,delimiter=',',skiprows=1)[:,1]
	#file='exp/x'+str(n)+'/d5/UUy.csv'
	#ve=loadtxt(file,delimiter=',',skiprows=1)[:,1]
	
	#file='kEpsilon/sets/0.252/Profiles'+str(n)+'_k.xy'
	#yg=loadtxt(file,delimiter=' ')[:,0]/D
	#ug=sqrt(0.667*loadtxt(file,delimiter=' ')[:,1])
	
	#ke=0.5*(ue**2+2*ve**2)
	#ye=array([0.0431603,0.11217,0.175433,0.242521,0.305758,0.370879,0.43393,0.498768,0.563718])
	#ue=array([30.9176,30.9176,31.0588,30.9176,30.2118,28.5176,22.0235,11.4353,4.37647])
	
	ls=[]
	dx=0.2*D
	
	for j in range(times):
		for i in range(len(data[j])):
			drop=data[j][i,:]
			if drop[22] > n*D and drop[22] < n*D+dx:
				ls.append([drop[3],drop[13],drop[14],drop[19],drop[20],drop[22]/D,drop[23]/D,drop[9],drop[2]])
	
	ls=array(ls)
	ls=parcel(ls[:,0],ls[:,1],ls[:,2],ls[:,3],ls[:,4],ls[:,5],ls[:,6],ls[:,7],ls[:,8])
	print 'len(ls) = ',len(ls.x)
	
	r=0.
	for i in range(len(ls.m)):
		r=r+ls.m[i]*ls.u[i]
	
	mflux.append(r/2./dx*60e3*360./5./times)
	print 'mflux [g/min] = ', mflux[-1]
	
fig = figure()
# plot data
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )

# plot data
mflux=array(mflux)
locs[0]=0.
ax.plot(locs,mflux/mflux[0],'-k',label='numerical')
ax.errorbar(exp[:,0], exp[:,1]/exp[0,1],fmt='sr',markeredgewidth=1.5,markerfacecolor='None',markeredgecolor='r',yerr=0.1*exp[:,1]/exp[0,1],label='experimental')
xlim([0,30])
#print 'error liquid flow rate = ',1.-(mflux/mflux[0])/(exp[:,1]/exp[0,1])

ylim([0,1.2])
xlabel('Axial position - x/D')
ylabel('Liquid mass flux - $\dot{m} / \dot{m}0$')

legend().draw_frame(False)

savetxt('drop_mflux_num.dat',vstack([locs,mflux/mflux[0]]).T,delimiter='\t')
savetxt('drop_mflux_exp.dat',vstack([exp[:,0], exp[:,1]/exp[0,1]]).T,delimiter='\t')
savefig('drop_mflux.png')

show()
