#!/usr/bin/env python
from pylab import *
import matplotlib.delaunay.triangulate as delaunay
import os,sys,string,progress
import visit_writer

class grid:
	def __init__(self, f,B):
		xr=f[:,0].copy()*1e-3
		yr=f[:,1].copy()*1e-3
		r=(xr**2+yr**2)**0.5
		x=[]
		y=[]
		data=[]
		for i in range(len(xr)):
			if (r[i] < 0.5*B):
				x.append(xr[i])
				y.append(yr[i])
				data.append(f[i,:])
		x=array(x)
		y=array(y)
		self.circumcenters, self.edges, self.tri_points, self.tri_neighbors = delaunay.delaunay(x, y)
		self.triangles=range(len(self.tri_points))
		self.mesh = delaunay.Triangulation(x,y)
		self.xp=x
		self.yp=y
		self.data=array(data)
#-----------------------------------------------------------------------
#		Centroids, Velocity Field and Cell Area
#-----------------------------------------------------------------------
		N=len(self.triangles)
		self.uvw=zeros([N,3])
		self.centroids=zeros([N,2])
		self.area=zeros(N)
		
		for i in range(len(self.triangles)):
			points=self.tri_points[i]
			x=(self.mesh.x[points[0]]+self.mesh.x[points[1]]+self.mesh.x[points[2]])/3.
			y=(self.mesh.y[points[0]]+self.mesh.y[points[1]]+self.mesh.y[points[2]])/3.
			self.centroids[i,0]=x
			self.centroids[i,1]=y
			# Area
			p1=[self.mesh.x[points[0]],self.mesh.x[points[1]],self.mesh.x[points[2]]]
			p2=[self.mesh.y[points[0]],self.mesh.y[points[1]],self.mesh.y[points[2]]]
			M=array([p1,p2,[1,1,1]])
			self.area[i]=0.5*abs(det(M))
		
	def volflow(self):
		s=0.
		for i in self.triangles:
			x=self.centroids[i,0]
			y=self.centroids[i,1]
			w=self.uvw[i,2]
			s=s+self.area[i]*w
		return s
#-----------------------------------------------------------------------
# Iz
#-----------------------------------------------------------------------
	def Iz(self):
		s=0.
		for i in self.triangles:
			x=self.centroids[i,0]
			y=self.centroids[i,1]
			s=s+self.area[i]*(x**2+y**2)
		return s
#-----------------------------------------------------------------------
# Ix
#-----------------------------------------------------------------------
	def Ix(self):
		s=0.
		for i in self.triangles:
			y=self.centroids[i,1]
			s=s+self.area[i]*(y**2)
		return s
#-----------------------------------------------------------------------
# Iy
#-----------------------------------------------------------------------
	def Iy(self):
		s=0.
		for i in self.triangles:
			x=self.centroids[i,0]
			s=s+self.area[i]*(x**2)
		return s
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#  Lx
#-----------------------------------------------------------------------
def LxF(malha,field):
	s=0.
	for i in malha.triangles:
		x=malha.centroids[i,0]
		y=malha.centroids[i,1]
		u=field.uvw[i,0]
		v=field.uvw[i,1]
		w=field.uvw[i,2]
		s=s+malha.area[i]*(w*y)
	return s
#-----------------------------------------------------------------------
#  Ly
#-----------------------------------------------------------------------
def LyF(malha,field):
	s=0.
	for i in malha.triangles:
		x=malha.centroids[i,0]
		y=malha.centroids[i,1]
		u=field.uvw[i,0]
		v=field.uvw[i,1]
		w=field.uvw[i,2]
		s=s+malha.area[i]*(-w*x)
	return s
#-----------------------------------------------------------------------
#  Lz
#-----------------------------------------------------------------------
def LzF(malha,field):
	s=0.
	for i in malha.triangles:
		x=malha.centroids[i,0]
		y=malha.centroids[i,1]
		u=field.uvw[i,0]
		v=field.uvw[i,1]
		w=field.uvw[i,2]
		s=s+malha.area[i]*(v*x-u*y)
	return s
#-----------------------------------------------------------------------
#  VolFlow
#-----------------------------------------------------------------------
def volflow(malha,field):
	s=0.
	for i in malha.triangles:
		x=malha.centroids[i,0]
		y=malha.centroids[i,1]
		u=field.uvw[i,0]
		v=field.uvw[i,1]
		w=field.uvw[i,2]
		s=s+malha.area[i]*w
	return s
#-----------------------------------------------------------------------

class fields:
	def __init__(self,N):
		self.uvw=zeros([N,3])
#-----------------------------------------------------------------------

def normalize(malha,x,k):
	s=0.
	for i in range(len(malha.triangles)):
		s=s+(x.uvw[i,0]**2+x.uvw[i,1]**2+x.uvw[i,2]**2)*malha.area[i]
		
	s=sqrt(k/s)
	for i in range(len(malha.triangles)):
		x.uvw[i,0]=x.uvw[i,0]*s
		x.uvw[i,1]=x.uvw[i,1]*s
		x.uvw[i,2]=x.uvw[i,2]*s
	return x

def innerproduct(malha,x,y):
	s=0.
	for i in range(len(malha.triangles)):
		s=s+(x.uvw[i,0]*y.uvw[i,0]+x.uvw[i,1]*y.uvw[i,1]+x.uvw[i,2]*y.uvw[i,2])*malha.area[i]
	return s
#-----------------------------------------------------------------------

def Q(malha,field):
	N=len(field)
	Qr=zeros([N,N])
	
	prog = progress.progressBar(0, N*N, 50)
	print prog, "\r",
	for i in range(N):
		for k in range(N):
			if i<=k:
				Qr[i,k]=innerproduct(malha,field[i],field[k])
			else:
				Qr[i,k]=Qr[k,i]
		prog.updateAmount(i*N+N)
		print prog, "\r",
		sys.stdout.flush()
	return Qr
#-----------------------------------------------------------------------
def RMS(malha,field):
	N=len(field)
	RMS=fields(len(malha.triangles))
	mean=average(malha,field)
	kk=1./sqrt(N)
	prog = progress.progressBar(0, N, 50)
	print prog, "\r",
	for i in range(N):
		for j in range(len(malha.triangles)):
			RMS.uvw[j,0]=RMS.uvw[j,0]+(field[i].uvw[j,0]-mean.uvw[j,0])**2
			RMS.uvw[j,1]=RMS.uvw[j,1]+(field[i].uvw[j,1]-mean.uvw[j,1])**2
			RMS.uvw[j,2]=RMS.uvw[j,2]+(field[i].uvw[j,2]-mean.uvw[j,2])**2
		prog.updateAmount(i*N+N)
		print prog, "\r",
		sys.stdout.flush()
	
	for j in range(len(malha.triangles)):
		RMS.uvw[j,0]=kk*sqrt(RMS.uvw[j,0])
		RMS.uvw[j,1]=kk*sqrt(RMS.uvw[j,1])
		RMS.uvw[j,2]=kk*sqrt(RMS.uvw[j,2])
		
	return RMS, innerproduct(malha,RMS,RMS)
	
	#-----------------------------------------------------------------------
def average(malha,field):
	N=len(field)
	average=fields(len(malha.triangles))
	
	for j in range(len(malha.triangles)):
		su,sv,sw=0.,0.,0.
		for i in range(N):
			su=su+field[i].uvw[j,0]
			sv=sv+field[i].uvw[j,1]
			sw=sw+field[i].uvw[j,2]
		average.uvw[j,0]=su/N
		average.uvw[j,1]=sv/N
		average.uvw[j,2]=sw/N
		
	return average

#-----------------------------------------------------------------------
def PODmodes(outputdir,malha,B,field,v,w,fr):
	ind=w.argsort()
	ind=ind[::-1]
	# Sorting arrays (greatest eigenvalues first)
	e=cumsum(w[ind])/sum(w)
	ep=w[ind]/sum(w)
	# Getting number of mesh cells and modes
	N=len(field)
	Ntri=len(field[0].uvw[:,0])
	# Creating Tecplot mesh
	xi=linspace(-B/2.,B/2.,50)
	xi,yi=meshgrid(xi,xi)
	
	x=malha.centroids[:,0]
	y=malha.centroids[:,1]
	
	i=0
	POD=[]
	b=[]
	while e[i]<fr:
		POD.append(fields(Ntri))
		for j in range(N):
			POD[i].uvw[:,0]=POD[i].uvw[:,0]+field[j].uvw[:,0]*v[j,ind[i]]
			POD[i].uvw[:,1]=POD[i].uvw[:,1]+field[j].uvw[:,1]*v[j,ind[i]]
			POD[i].uvw[:,2]=POD[i].uvw[:,2]+field[j].uvw[:,2]*v[j,ind[i]]
		
		bj=zeros(N)
		for j in range(N):
			bj[j]=sqrt(ep[i])*innerproduct(malha,POD[i],field[j])/innerproduct(malha,POD[i],POD[i])
		
		POD[i]=normalize(malha,POD[i],1.)
		uu=POD[i].uvw[:,0]
		vv=POD[i].uvw[:,1]
		ww=POD[i].uvw[:,2]
		
		print 'wrinting field for mode E['+str(i)+'] = ',e[i]
		pts=[]
		vec=[]
		for m in range(len(x)):
			pts=pts+[x[m],y[m],0]
			vec=vec+[float(uu[m]),float(vv[m]),float(ww[m])]
			
		vars = (("data", 3, 1, vec), ("ptsvec", 3, 1, pts))
		visit_writer.WritePointMesh(outputdir+'POD'+str(i)+'.vtk', 1, pts, vars)
		savetxt(outputdir+'a'+str(i)+'.dat',v[:,ind[i]])
		b.append(bj)
		i=i+1
	b=array(b).T
	savetxt(outputdir+'b.dat',b)
#-----------------------------------------------------------------------


# Open field files
# --------------------
def readbase(file,bore):
	f=loadtxt(file,skiprows=4)
	malha=grid(f,bore)
	return malha

def main(B,N,inputdir,outputdir):
# -------------------------------
# SYGMA POD
# -------------------------------
	print('# SYGMA POD RC1')
# -------------------------------
# Reading mesh
# -------------------------------
	print '\nCleaning output dir: \t',outputdir
	print '\nOpening input dir: \t',inputdir
	os.system('rm '+outputdir+'*')
	f=loadtxt(inputdir+'/Export000.dat',skiprows=4)
	malha=grid(f,B)
	field=[]
	KE=[]
# -------------------------------
# Interpolating each field to the created mesh
# -------------------------------
	print '\nReading and interpolating each field:'
	prog = progress.progressBar(0, N, 50)
	for i in range(N):
		if i<10:
			file=inputdir+'Export00'+str(i)+'.dat'
		elif i>=10 and i<100:
			file=inputdir+'Export0'+str(i)+'.dat'
		elif i>=100:
			file=inputdir+'Export'+str(i)+'.dat'
		
		f=loadtxt(file,skiprows=4)
		xr=f[:,0].copy()*1e-3
		yr=f[:,1].copy()*1e-3
		r=(xr**2+yr**2)**0.5
		data=[]
		for k in range(len(xr)):
			if (r[k] < 0.5*B):
				data.append(f[k,:])
		data=array(data)
		field.append(fields(len(malha.triangles)))
		itpu=malha.mesh.linear_interpolator(data[:,3])
		itpv=malha.mesh.linear_interpolator(data[:,4])
		itpw=malha.mesh.linear_interpolator(data[:,5])
		
#		print 'Field:\t', i
		for j in range(len(malha.triangles)):
			x=malha.centroids[j,0]
			y=malha.centroids[j,1]
			field[i].uvw[j,0]=itpu.planes[j,0]*x + itpu.planes[j,1]*y + itpu.planes[j,2]
			field[i].uvw[j,1]=itpv.planes[j,0]*x + itpv.planes[j,1]*y + itpv.planes[j,2]
			field[i].uvw[j,2]=itpw.planes[j,0]*x + itpw.planes[j,1]*y + itpw.planes[j,2]
		
		prog.updateAmount(i+1)
		print prog, "\r",
		sys.stdout.flush()
		KE.append(innerproduct(malha,field[i],field[i]))
		
	print '\n Computing covariance matrix:'
	COV=Q(malha,field)
	print '\n Computing eigenvalues. Printing energy fraction of each mode: \n'
	w,v=eig(COV)
	ind=w.argsort()
	frac=w[ind]/sum(w)
	savetxt(outputdir+'eigenvalues.dat',w)
	savetxt(outputdir+'energy_fraction.dat',frac)
	KE=array(KE)
	TKE=sum(KE)
	
#	spectrum=fft(KE)
	print frac[:]
	print '\n kinetic energy (eigenvalue) [J/m/kg] =\t\t',sum(w)
	print ' kinetic energy (flow integral) [J/m/kg] =\t', TKE
	print '\n Writing POD modes to files.'
	PODmodes(outputdir,malha,B,field,v,w,0.95)
	
#	plot(arange(N),cumsum(frac[::-1]),'ob')
#	xlabel('POD Mode')
#	ylabel('Cumulative Energy Fraction')
#	title('POD Energy Spectrum')
#	title('Scania - Lift 12mm - 1D Plane')
#	ylim([0.8,1.])
#	savefig('energy_fraction.png')
	
#	show()

#inputdir='/home/piccinini/work/POD/python/action/snaps/'
#outputdir='/home/piccinini/work/POD/python/action/pods/'
#main(0.105,100,inputdir,outputdir)

#inputdir='/home/piccinini/work/POD/python/action_D_lift12/snaps/'
#outputdir='/home/piccinini/work/POD/python/action_D_lift12/pods/'
#main(0.105,100,inputdir,outputdir)

#inputdir='/home/piccinini/work/POD/python/ducato_D_lift2/snaps/'
#outputdir='/home/piccinini/work/POD/python/ducato_D_lift2/pods/'
#main(0.0944,100,inputdir,outputdir)

#inputdir='/home/piccinini/work/POD/python/ducato_D_lift12/snaps/'
#outputdir='/home/piccinini/work/POD/python/ducato_D_lift12/pods/'
#main(0.0944,100,inputdir,outputdir)

#inputdir='/home/piccinini/pod_encit/trunk/python/sygma/snaps/'
#outputdir='/home/piccinini/pod_encit/trunk/python/sygma/pods/'
#main(0.128,100,inputdir,outputdir)

inputdir='/home/piccinini/works/pod_encit/trunk/python/scania/snaps/'
outputdir='/home/piccinini/works/pod_encit/trunk/python/scania/pods/'
main(0.128,80,inputdir,outputdir)

#inputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/dir_close_esq_open/snaps/'
#outputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/dir_close_esq_open/pods/'
#main(0.0944,110,inputdir,outputdir)

#inputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/duas_abertas/snaps/'
#outputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/duas_abertas/pods/'
#main(0.0944,110,inputdir,outputdir)

#inputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/duas_meio_abertas/snaps/'
#outputdir='/home/piccinini/pod_encit/trunk/python/andreotti/dados/duas_meio_abertas/pods/'
#main(0.0944,110,inputdir,outputdir)
