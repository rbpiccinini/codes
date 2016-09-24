#!/usr/bin/env python
from pylab import *
import scipy.optimize as opt
import numpy as np
import glob

#-----------------------------------------------------------------------

class case:
    def __init__(self, f,name,style,outname,color):
        self.file=f
        self.name=name
        self.style=style
        self.outname=outname
        self.color=color


def func(p,x):
   w=p[1]*np.exp(p[0]*x)
   print w.size
   return w

def residuals(p,x,y):
   w=p[1]*np.exp(p[0]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B

# MWM Action = Head 1
# Scania DC-9/12 = Head 2
# Sygma = Head 3

file=[]
file.append(case('/home/piccinini/work/POD/python/action','Head 1 L=2mm','o-','H1L2','#1518BE'))
file.append(case('/home/piccinini/work/POD/python/action_D_lift12','Head 1 L=12mm','x-','H1L12','#1518BE'))
file.append(case('/home/piccinini/work/POD/python/scania','Head 2 L=12mm','s-','H2L12','#F79B1C'))
file.append(case('/home/piccinini/work/POD/python/sygma','Head 3 L=12mm','+-','H3L12','#008F2D'))
file.append(case('/home/piccinini/work/POD/python/ducato_D_lift2','Head 4 L=2mm','v-','H4L2','#D3291A'))
file.append(case('/home/piccinini/work/POD/python/ducato_D_lift12','Head 4 L=12mm','^-','H4L12','#D3291A'))

#-----------------------------------------------------------------------

fig=figure()
ax1 = fig.add_axes([0.1, 0.1, 0.45, 0.8])
ax2 = fig.add_axes([0.65, 0.1, 0.3, 0.4])
plots1=[]
plots2=[]
leg=[]
for i in file:
    x=loadtxt(i.file+'/energy_fraction.dat')
    ind=x.argsort()
    ind=ind[::-1]
    # Sorting arrays (greatest eigenvalues first)
    e=x[ind]/sum(x)
    plots1.append(ax1.semilogy(e,i.style,color=i.color))
    plots2.append(ax2.semilogy(e,i.style,color=i.color))
    leg.append(i.name)
    
    print len(arange(30,len(x))),len(x[ind[30:]])
    p0=[0.3,2.]
    plsq = opt.fmin(residuals, p0, args=(arange(30,len(x)),x[ind[30:]]), maxiter=10000, maxfun=10000)
    #plsq = opt.leastsq(residuals, p0, args=(x,y))
    print 'tau= ',-1./plsq[0], plsq[1]
    

fig.legend(plots1,leg,'upper right')
setp(ax1, xlim=(0,100))
setp(ax2, xlim=(-1,10),ylim=(1e-3,1))
ax1.grid(1)
ax2.grid(1)

setp(ax1,ylabel='Fraction of Energy Content',xlabel='POD mode')
savefig('energy_spectra.pdf')



#-----------------------------------------------------------------------
fig=figure()

ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.6])

leg=[]
ind=arange(len(file))
phi0=[]
phi1=[]
ticname=[]
for i in file:
    x=loadtxt(i.file+'/energy_fraction.dat')
    phi0.append(x[-1])
    phi1.append(x[-2])
    leg.append(i.name)
    ticname.append(i.outname)

phi0=array(phi0)
phi1=array(phi1)
width=0.6
p1=ax1.bar(ind, phi0,   width, color='#27408B')
p2=ax1.bar(ind, phi1,   width, color='#436EEE',bottom=phi0)
p3=ax1.bar(ind, 1.-phi0-phi1,   width, color='#5CACEE',bottom=phi1+phi0)

xticks(ind+width/2., ticname)

fig.legend((p1[0],p2[0],p3[0]),('$\Phi_1$','$\Phi_2$','$\Phi_3+...+\Phi_N$'),'upper center')

#plot(x[-1::-1],i.style,color=i.color)
#leg.append(i.name)
#legend(leg)
#xlim([-0.5,100])
#ylim([1e-3,1])
#grid(1)

ylabel('Fraction of Average Energy Content')
#xlabel('POD mode')
savefig('energy_fraction.pdf')

#-----------------------------------------------------------------------
figure()
leg=[]
for i in file:
    x=loadtxt(i.file+'/energy_fraction.dat')
    ind=x.argsort()
    ind=ind[::-1]
    # Sorting arrays (greatest eigenvalues first)
    e=cumsum(x[ind])/sum(x)
    plot(range(1,len(e)+1),e,i.style,color=i.color)
    leg.append(i.name)

legend(leg,'lower right')
xlim([0,100])
ylim([0.4,1])
grid(1)

ylabel('Cumulative Sum of Energy Content')
xlabel('POD modes')
savefig('cumulative_energy.pdf')

#-----------------------------------------------------------------------
leg=[]
k=0
for i in file:
    k=k+1
    b=loadtxt(i.file+'/b.dat')
    e=sum(b**2,axis=1)*len(b)*0.95
    e0=b[:,0]**2*len(b)*0.95
    e5=sum(b[:,:5]**2,axis=1)*len(b)*0.95
    e10=sum(b[:,:10]**2,axis=1)*len(b)*0.95
    
    figure()
#    Sorting arrays (greatest eigenvalues first)
#    plot(e,'-k')
#    plot(b0**2*len(b0)+b1**2*len(b1),'--k')
    plot(range(1,len(e0)+1),e0/e,'0.55')
    plot(range(1,len(e10)+1),e10/e,'k')
#    plot(e5/e,'k')
#    plot(e,'k')
    title(i.name)
#    ylim([.6,1.])
    grid(1)
    ylim([0.3,1.0])
#   xlim([1,len(e0)+1])
    ylabel('Energy Fraction')
    xlabel('Time [ms]')
    if k!=1:
        legend(('$\Phi_1$', '$\Phi_1+...+\Phi_{10}$'),'lower right')
    else:
        legend(('$\Phi_1$', '$\Phi_1+...+\Phi_{10}$'),'upper right')
    savefig('energy'+str(k)+'.pdf')


#-----------------------------------------------------------------------

#leg=[]
#k=0
#for i in file:
    #k=k+1
    #b=loadtxt(i.file+'/b.dat')
    #e=sum(b**2,axis=1)*len(b)*0.95
    #e0=b[:,0]**2*len(b)*0.95
    #e5=sum(b[:,:5]**2,axis=1)*len(b)*0.95
    #e10=sum(b[:,:10]**2,axis=1)*len(b)*0.95
    
    #figure()
    ## Sorting arrays (greatest eigenvalues first)
##    plot(e,'-k')
##    plot(b0**2*len(b0)+b1**2*len(b1),'--k')
    #plot(sqrt(e0/e),'0.8')
    #plot(sqrt(e10/e),'k')
##    plot(e,'k')
    #title(i.name)
##    ylim([.6,1.])
    #grid(1)
    #ylim([0.3,1.0])
    #ylabel('Energy Fraction')
    #xlabel('Time (ms)')
    #legend(('First Mode', '10 First Modes'),'lower right')
    #savefig('ulinha'+str(k)+'.pdf')


#-----------------------------------------------------------------------
# Lz components
#-----------------------------------------------------------------------

mode_number=3
color=['#3A5FCD','0.2','0.4','0.6','0.8']
for mode_number in range(len(file)):
    fig=figure()
    ax1 = fig.add_axes([0.1, 0.1, 0.6, 0.8])
    plots=[]
    
    b=loadtxt(file[mode_number].file+'/b.dat')
    Lz=b[:,:5]*len(b)
    
    for i in range(5):
        Lz[:,i]=Lz[:,i]*float(loadtxt(file[mode_number].file+'/Lz/Lz.'+str(i)+'.csv',skiprows=1,delimiter=',')[0])
        plots.append(ax1.plot(range(1,len(Lz[:,i])+1),1.e3*Lz[:,i],color=color[i]))
    
    xlabel('Times [ms]')
    ylabel('$10^{-3}$ Angular Momentum [$m^2/s$]')
    title(file[mode_number].name)
    
    fig.legend(plots,('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5'),'center right')
    
    grid(1)
    savefig('Lz_'+file[mode_number].outname+'.pdf')


#-----------------------------------------------------------------------

#show()


