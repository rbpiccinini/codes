#!/usr/bin/env python

from pylab import *

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
	
data=loadtxt('drops_RNG.csv',skiprows=1,delimiter=',')
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

for n in [5,10,15,20,25]:
	file='exp/x'+str(n)+'/Ux.csv'
	ye=loadtxt(file,delimiter=',')[:,0]
	ue=loadtxt(file,delimiter=',')[:,1]
	
	#ye=array([0.0431603,0.11217,0.175433,0.242521,0.305758,0.370879,0.43393,0.498768,0.563718])
	#ue=array([30.9176,30.9176,31.0588,30.9176,30.2118,28.5176,22.0235,11.4353,4.37647])
	
	ls=[]
	for i in range(len(data)):
		drop=data[i,:]
		if drop[22]/D > n-.35 and drop[22]/D < n+.35:
			ls.append([drop[3],drop[13],drop[14],drop[19],drop[20],drop[22]/D,drop[23]/D,drop[9],drop[2]])
	
	ls=array(ls)
	ls=parcel(ls[:,0],ls[:,1],ls[:,2],ls[:,3],ls[:,4],ls[:,5],ls[:,6],ls[:,7],ls[:,8])
	print 'len(ls) = ',len(ls.x)
	
	print max(ye)
	yc=linspace(0.,max(ls.y),30)
	ls.reynolds(yc,0.1)
	
	#print min(ls.x),max(ls.x)
	#figure(1)
	#plot(n,sum(ls.m[:]),'bo')
	#title('Liquid Mass [kg]')
	#xlabel('x/D')
	
	figure()
	plot(ls.yc,ls.uuc,'-b',linewidth=2)
	plot(ls.y,ls.uu,'.k')
	plot(ye,ue,'o-r')
	xlim([0,max(ye)])
	ylim([0,max(ue)])
	xlabel('y/D')
	ylabel('Axial Velocity - U (m/s)')
	title('Axial Position - x/D = '+str(n))

show()
