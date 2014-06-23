#!/usr/bin/python

from pylab import *

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
## rc('font',family='serif')

rc('font',family='serif')
rc('text', usetex=True)
rc('font', size=16)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
rc('legend', fontsize=16)
rc('figure',figsize=[8,6])

fbottom=0.15
fleft=0.15

def errol(rho,c,z):
    g=9.81
    
    dp0=rho*g*z
    dp1=1/c*(exp(c*dp0)-1.)
    return (dp0-dp1)/dp1

def errog(rho,p0,z):
    g=9.81
    
    dp0=rho*g*z
    dp1=p0*(exp(dp0/p0)-1.)
    return (dp0-dp1)/dp1


g=9.81
z=linspace(1e-3,500.,10)
c1=3.6e-6/(g*1e4)*14.23
c2=50e-6/(g*1e4)*14.23

rhoo=1e3
l1=errol(rhoo,c1,z)
l2=errol(rhoo,c2,z)

p0=100e5
rhog=p0/1e5
g1=errog(rhog,p0,z)
g2=errog(rhog,p0*5,z)

figure()
#plot(z, y2-1,'r',label='approx 0th')
plot(z, 100*l1,'k-',label=r'$\rho_o =1000\ kg/m^3$ and $c_o=10^{-3}\ 1/MPa$')
plot(z, 100*l2,'sk-',label=r'$\rho_o =1000\ kg/m^3$ and $c_o=2 \times 10^{-3}\ 1/MPa$')

plot(z, 100*g1,'sk--',label=r'$\rho_g =100\ kg/m^3$ and $p^{ref}_g = 10 MPa$')
plot(z, 100*g2,'k--',label=r'$\rho_g =100\ kg/m^3$ and $p^{ref}_g = 50 MPa$')


legend(loc='lower left')
xlabel('Depth, m')
ylabel('Error, $\%$')
ylim([-4.,0])
grid(1)
subplots_adjust(bottom=fbottom,left=fleft,right=0.85)
print max(l1), max(l2), c1, c2
grid(1)

savefig('density.pdf',dpi=300)
show()
