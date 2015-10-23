from pylab import *
from scipy import stats
from scipy.interpolate import splrep, splev
from scipy import optimize as opt

class classBuildup:
	def __init__(self,t,p,sampling=100,nEV=2,r2=0.95):
		self.p0=p
		self.t0=t
		self.t=logspace(-2.,log10(t[-1]),sampling,base=10.)
		self.p=interp(self.t,t,p)
		self.lam=zeros(nEV+1)
		self.c=zeros(nEV+1)
		self.r2=r2
		self.k=[len(self.t)]
		self.dt=gradient(self.t)
		self.dp=gradient(self.p,self.dt)
		self.dpm=[]
		self.dpm.append(self.dp)
		for i in range(nEV):
			lam,c=self.getEV(self.t,self.dpm[-1],r2=self.r2)
			self.dpm.append(self.dpm[-1]-lam*c*exp(-lam*self.t))
			self.lam[i+1]=lam
			self.c[i+1]=c


	def residual(self,c,lam,lndp,t):
		return min(lndp+lam*t-log(c*lam))**2

	def getEV(self,t, dp, r2=0.95):
		lndp=log(dp)
		n=len(lndp)
		r2_=0.
		plot(t[:self.k[-1]],lndp[:self.k[-1]], 'o')
#		show()
		self.k.append(argmin(abs(t-input('Qual o tempo inicial para ajustar exponencial?'))))
		r=stats.linregress(t[self.k[-1]:self.k[-2]],lndp[self.k[-1]:self.k[-2]])
		r2_=r[2]**2
		lam=-r[0]
		c=exp(r[1])/lam
#20		c=opt.fmin(self.residual,c,args=(lam,t[self.k[-1]:self.k[-2]],lndp[self.k[-1]:self.k[-2]]),ftol=0.001)
		plot(t[self.k[-1]:self.k[-2]],lndp[self.k[-1]:self.k[-2]], 'o')
		plot(t[self.k[-1]:self.k[-2]],-lam*t[self.k[-1]:self.k[-2]]+log(c*lam), '-')
#		show()
		return lam,c

	def getPress(self,t):
		self.textrap=t
		self.pi=self.p[-1]
		for i in range(len(self.lam)):
			self.pi=self.pi+self.c[i]*exp(-self.lam[i]*self.t[-1])

		self.pextrap=self.pi
		for i in range(len(self.lam)):
			self.pextrap=self.pextrap-self.c[i]*exp(-self.lam[i]*t)



data=loadtxt('dados_pressao_TLD.dat', skiprows=1)
tx=data[:,0]-data[0,0]
px=data[:,1]/14.23

siri=classBuildup(tx,px,sampling=1000,nEV=4,r2=0.95)
siri.getPress(logspace(-2.,log10(tx[-1])+1,200,base=10.))

print 'lambda \t = \t ',siri.lam
print 'c \t = \t',siri.c


## pi = siri.p[-1]+c*exp(-lam*siri.t[-1])
## print 'pf [kgf/cm2] = ', siri.p[-1]

print 'pfinal  [kgf/cm2] = ', siri.p[-1]
print 'pextrap [kgf/cm2] = ', siri.pi
print 'tempo final [h] = ', tx[-1]

figure()
plot(siri.textrap, siri.pextrap,'or')
plot(tx,px,'-')
grid(1)
ylim([px[0],100])
xlim([0.,siri.textrap[-1]])
show()
