from pylab import *
from numpy.linalg import eig
from scipy.sparse.linalg import eigsh
import numpy as np
from scipy import stats
from scipy import integrate

params = {
'backend': 'ps',
#'text.latex.preamble': ['\usepackage[latin1]{inputenc}'],
'axes.labelsize': 8, # fontsize for x and y labels (was 10)
'axes.titlesize': 8,
'axes.grid': False,
'font.size': 8, # was 10
'legend.fontsize': 8, # was 10
'xtick.labelsize': 8,
'ytick.labelsize': 8,
'text.usetex': True,
#'figure.figsize': [6.9,8.0],
'figure.figsize': [3.3, 2.04],
'font.family': 'sans-serif',
'font.sans-serif': 'Helvetica',
'lines.linewidth': 1,
'lines.markersize': 3
}

rcParams.update(params)

def householder(A):
	n = len(A)
	for k in range(n-2):
		alpha = np.sign(A[k+1,k])*np.sqrt(sum(A[k+1:,k]**2))
		r = np.sqrt(0.5*(alpha**2-A[k+1,k]*alpha))
		v = zeros(n)
		v[k+1] = (A[k+1,k]-alpha)/(2.0*r)
		v[k+2:] = A[k+2:,k]/(2.0*r)
		P = np.identity(n)-2.0*np.outer(v,v)
		A = dot(P,dot(A,P))

	a = zeros(n-1)
	b = zeros(n)
	c = zeros(n-1)
	
	# a \ b \ c
	for i in range(1,n-1):
		a[i-1] = A[i,i-1]
		b[i] = A[i,i]
		c[i] = A[i, i+1]
	
	b[0] = A[0,0]
	b[n-1] = A[n-1,n-1]

	c[0] = A[0,1]
	a[n-2] = A[n-1,n-2]
	
	A = np.diag(a,-1)+np.diag(b, 0)+np.diag(c,1)
	
	return A
	 

def eigCalc(n, L, k):
	x = linspace(0.,L,n)
	dx = L/n
	
	# wn=1
	# k=(1.0 + 0.7*sin(wn*pi*x))*1e-2

	kph = zeros(n-1)
	for i in range(n-1):
		kph[i]=2.0*k[i]*k[i+1]/(k[i]+k[i+1])
	
	kmh = zeros(n)
	for i in range(1,n):
		kmh[i]=2.0*k[i]*k[i-1]/(k[i]+k[i-1])
		
	A=zeros([n,n])
	a = zeros(n-1)
	b = zeros(n)
	c = zeros(n-1)
	
	# a \ b \ c
	for i in range(1,n-1):
		a[i-1] = kmh[i] 
		b[i] = -kph[i]-kmh[i]
		c[i] = kph[i]

#		A[i,i-1] = kmh[i]
#		A[i,i] = -kph[i]-kmh[i]
#		A[i,i+1] = kph[i] 

		b[0] = -kph[0]
		c[0] = kph[0]

		a[n-2] = kmh[n-1]
		b[n-1] =-kmh[n-1] #-A[n-1,n-2]


	a = -a/dx**2
	b = -b/dx**2
	c = -c/dx**2
	
	
	A = np.diag(a,-1)+np.diag(b, 0)+np.diag(c,1)
	
	nest=100
	w, v = eigsh(A, k=nest, which='SM')
	k = argsort(w)
	
	
#     reconstructing k(x)
#    s = zeros(n)
#    for i in range(nest):
#        s = s + w[i]*v[:,i]
#    s = integrate.cumtrapz(s, dx=dx, initial=0.)
#    kest = -s*dx/gradient(sum(v[:,1:],axis=1))
	kest = recK(A, dx)
	
	return w, v, kest, A

def recK(A, dx):
	#     reconstructing k(x)
	n = len(A)
	
	k=zeros(n)
	kph = zeros(n)
	kmh = zeros(n)

	for i in range(1,n-1):
		k[i] = 0.5*A[i,i]
	k[0] = A[0,0]
	k[n-1] = A[n-1,n-1] #A[n-1,n-1]
#	for i in range(n-1):
#		kph[i]= A[i, i+1]
#		kmh[i+1]= A[i+1, i]

#	kph[n-1] = kmh[n-1]
#	kmh[0] = kph[0]
#
#	k = -2.0*kph*kmh/(kph+kmh)*dx**2

	return k*dx**2
	

def calcvho(x,n):
	L=max(x)
	vho=zeros([len(x),n])
	for i in range(n):
		vho[:,i] = cos(i*pi/L*x)
		vho[:,i] = vho[:,i]/sqrt(dot(vho[:,i],vho[:,i]))
	return -vho
				
#print(eigs[1:10:1]/(pi*arange(1,10,1)/L)**2)

#for n in [5,10,50,100,250,500,1000,2000]:
L = 1.
wn = 3
f=1.0

n = 500
m = 200
x = linspace(0.,L,n) 
khe1=(1.0 + 0.7*sin(wn*pi*x))*1e-2 #np.random.random(n)
kho = L/trapz(1.0/khe1,x)
khe2=array(349*[1.]+[35.]+[70.]+149*[100.])
khe2=0.01+(100-0.01)/(1.0+np.exp(-100.0*(x-0.7)))
khe2=khe2/(1.0/trapz(1.0/khe2,x))*kho

print('kho = ', kho, 1.0/trapz(1.0/khe2,x))

w1, v1, kest1, A1 = eigCalc(n,L=L, k=khe1)
w2, v2, kest2, A2 = eigCalc(n,L=L, k=khe2)
who, vho, kestho, Aho = eigCalc(n,L=L, k=kho*ones(n))

#Arec = Aho #zeros([n,n])
#nest = 100
#for i in range(nest):
#	Arec = Arec + (w[i]-who[i])*outer(vho[:,i],vho[:,i])
#Arec = householder(Arec)
#wrec, vrec = eigsh(Arec, k=10, which='SM')
#krec = recK(Arec, L/n)

figure(figsize=[3.3, 3.3])
plot(x,kho*gradient(vho[:,5])**2/gradient(x)**2,'k-')
fill(x,khe1*gradient(v1[:,5])**2/gradient(x)**2, color='black', alpha=0.25)
fill(x,khe2*gradient(v2[:,5])**2/gradient(x)**2, color='black', alpha=0.50)
xlim([0.,1.])
xlabel('Position ($m$)')
ylabel('$\eta(x)(dPsi/dx)^2$ ($kPa^2/s$)')
subplots_adjust(top=0.8)
legend(('homogeneous', 'weakly heterogeneous', 'strongly heterogeneous'), loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=1, borderaxespad=0, frameon=False)
#tight_layout()
savefig('localization.jpg', dpi=300, additional_artists=[], bbox_inches="tight")
show()



