#!/usr/bin/env python

from pylab import *

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
			for i in range(len(self.y)):
				if self.y[i] > y[j]-dy and self.y[i] < y[j]+dy:
					ut.append(self.u[i])
					uut.append(self.uu[i])
					vt.append(self.v[i])
					vvt.append(self.vv[i])
			ut=array(ut)
			uut=array(uut)
			vt=array(vt)
			vvt=array(vvt)
			self.uc[j]=average(ut)
			self.uuc[j]=sqrt(average(uut**2))
			self.vc[j]=average(vt)
			self.vvc[j]=sqrt(average(vvt**2)) 
	
data=loadtxt('drops.csv',skiprows=1,delimiter=',')
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

for n in range(5,35,5):
	file='exp/diameter/z'+str(n)
	exp=loadtxt(file,delimiter=',')
	for i in range(1,exp.shape[1]):
		exp[1:,i]=exp[1:,i]/sum(exp[1:,i])
	
	ls=[]
	for i in range(len(data)):
		drop=data[i,:]
		if drop[22]/D > n-.25 and drop[22]/D < n+.25:
			ls.append([drop[3],drop[13],drop[14],drop[19],drop[20],drop[22]/D,drop[23]/D,drop[9],drop[2]])
	
	ls=array(ls)
	ls=parcel(ls[:,0],ls[:,1],ls[:,2],ls[:,3],ls[:,4],ls[:,5],ls[:,6],ls[:,7],ls[:,8])
	print 'len(ls) = ',len(ls.x)
	
	dc=exp[0,:]
	
	print min(ls.x),max(ls.x)
	#plot(yc,uc,'ob')
	figure()
#	plot(ls.yc,ls.vc,'.b')
	plot(exp[1:,0],exp[1:,1],'or')
	xlabel('Diameter (1E-6 m)')
	ylabel('Frquency')
	title('Axial Position - x/D = '+str(n))

show()
