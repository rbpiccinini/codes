from pylab import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
## rc('font',family='serif')

rc('text', usetex=True)
rc('font', size=16)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
rc('legend', fontsize=16)
rc('figure',figsize=[8,6])

fbottom=0.15
fleft=0.15

rho=1e3
g=9.81
z=linspace(0,500.,100)
c1=100e-11
c2=50e-6/6894.75

x1=rho*g*z*c1
x2=rho*g*z*c2
## x=0.1*100e-6*z

y1=1./(1.-x1)
y2=1./(1.-x2)

figure()
#plot(z, y2-1,'r',label='approx 0th')
plot(z, 100*(y1-1),'k-',label=r'$\rho =1000\ kg/m^3$ and $c=10^{-9}\ Pa$')
plot(z, 100*(y2-1),'k--',label=r'$\rho =1000\ kg/m^3$ and $c=2 \times 10^{-9}\ Pa$')
legend(loc='upper left')
xlabel('Depth, m')
ylabel('Error, $\%$')
grid(1)
subplots_adjust(bottom=fbottom,left=fleft,right=0.85)
print max(y1), max(y2), c2

grid(1)
show()
