#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      BFV8
#
# Created:     08/11/2013
# Copyright:   (c) BFV8 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Sidafa Conde
# Step 1: 1D convection equation

import numpy as np
import matplotlib.pyplot as pl

Lx=2.
c = 1
nx = 50
nt = 40
dt = 0.001
dx = 2./(nx-1)

print 't = ',Lx/c
print 's = ',dt*nt*c
print 'CFL = ',c*dt/dx

x=np.linspace(0.,2.,nx)
x = 1.0*x
u0 = 2.0*(x>0.5)*(x<1.0)
pl.figure()
PLT = pl.plot(x,u0,'-sk',linewidth=1)
pl.hold(True)
#pl.clf()

for it in range(1,nt):
    un = u0
    for i in range(2,nx):
        u0[i] = un[i] - c*dt/dx*(un[i] - un[i-1])

pl.plot(x,u0,'-k',linewidth=1)
pl.hold(True)

pl.show()