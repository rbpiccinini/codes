
import firedrake as fd
import numpy as np
from slepc4py import SLEPc
from firedrake.petsc import PETSc
import scipy.io as sp
import pandas as pd
import matplotlib.pyplot as plt


def coords2ijk(x, y, z, Delta, data_array):
    i = np.floor(x / Delta[0]).astype(int)
    j = np.floor(y / Delta[1]).astype(int)
    k = np.floor(z / Delta[2]).astype(int)
    if i == data_array.shape[0]:
        i = i-1
    if j == data_array.shape[1]:
        j = i-1
    if k == data_array.shape[2]:
        k = k-1

    return data_array[i, j, k]

# Define time steps
tp = 48. # prod time in hours
t = np.linspace(0.,10*24*3600, 4801)
#t = list(np.logspace(np.log10(0.01*3600), np.log10(tp*3600), 80))
#t = t+ list(tp*3600+np.logspace(np.log10(0.01*3600), np.log10(20.*24*3600), 150))
#t = np.array(t)


dts = np.diff(t)
dt = fd.Constant(dts[0]) # time step size

# Define mesh
Delta_x = 20*0.3048 # ft -> m
Delta_y = 10*0.3048
Delta_z = 2*0.3048

Nx = int(60)
Ny = int(220)

layers = range(5) # range(35,40)
Nz = len(layers)

Lx = Nx*Delta_x
Ly = Ny*Delta_y
Lz = Nz*Delta_z

mesh = fd.utility_meshes.RectangleMesh(nx=Nx, ny=Ny, Lx=Lx, Ly=Ly,
                                       quadrilateral=True)
mesh = fd.ExtrudedMesh(mesh, layers=Nz, layer_height=Delta_z)

# We need to decide on the function space in which we'd like to solve the
# problem. Let's use piecewise linear functions continuous between
# elements::

V = fd.FunctionSpace(mesh, "CG", 1)
Vvec = fd.VectorFunctionSpace(mesh, "CG", 1)

# We'll also need the test and trial functions corresponding to this
# function space::

u = fd.TrialFunction(V)
v = fd.TestFunction(V)

# We declare a function over our function space and give it the
# value of our right hand side function::

# Permeability tensor
print("START: Read in reservoir fields")
Delta = np.array([Delta_x, Delta_y, Delta_z])
ct = 79.08e-11 # (1.0+0.2*3+0.8*4.947)*14.2 * 10**-6 kgf/cm2
mu = 0.003 # Pa-s

coords = fd.project(mesh.coordinates, Vvec).dat.data
spe10_2 = sp.loadmat('./spe10/spe10_2.mat')

kx_array = spe10_2['Kx'][:,:, layers]
ky_array = spe10_2['Ky'][:,:, layers]
kz_array = spe10_2['Kz'][:,:, layers]
phi_array = spe10_2['p'][:,:, layers]

Kx = fd.Function(V)
Ky = fd.Function(V)
Kz = fd.Function(V)
phi = fd.Function(V)

Tx = fd.Function(V)
Ty = fd.Function(V)
Tz = fd.Function(V)
w = fd.Function(V)

Kx.rename('Kx')
Ky.rename('Ky')
Kz.rename('Kz')
phi.rename('phi')

Tx.rename('Kx')
Ty.rename('Ky')
Tz.rename('Kz')
w.rename('w')

coords2ijk = np.vectorize(coords2ijk, excluded=['data_array', 'Delta'])

Kx.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kx_array)
Ky.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=ky_array)
Kz.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kz_array)
phi.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                       coords[:, 2], Delta=Delta, data_array=phi_array)
Tx.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kx_array/mu)
Ty.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=ky_array/mu)
Tz.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kz_array/mu)

w.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                       coords[:, 2], Delta=Delta, data_array=phi_array*ct)
print("END: Read in reservoir fields")

# Define initial value
x, y, z = fd.SpatialCoordinate(mesh)
eps = 1.
xc = 0.5*Delta_x*Lx
yc = 0.5*Delta_y*Ly
r = ((x-xc)**2+(y-yc)**2)**0.5 

# u_0 = fd.Function(V).interpolate(10*y/(Delta_y*Ny))
pi = 400e5
u_0 = fd.Function(V).interpolate(fd.Constant(pi)) # 400 bar
u_0.rename('Pressure')

u_n = fd.Function(V).assign(u_0)
u_n.rename('Pressure')

# Compute volume
vol = fd.Function(V).assign(fd.Constant(1.0))
vol = fd.assemble(vol*fd.dx)

# Write properties
fd.File('spe10_u0.pvd').write(u_n, Kx, Ky, Kz, phi)
outfile = fd.File("spe10_trans.pvd")
outfile.write(u_n)

# Compute transmissibility with harmonic interpolation of K to facets
Tx_facet = fd.conditional(fd.gt(fd.avg(Tx), 0.0), Tx('+')*Tx('-') / fd.avg(Tx), 0.0)
Ty_facet = fd.conditional(fd.gt(fd.avg(Ty), 0.0), Ty('+')*Ty('-') / fd.avg(Ty), 0.0)
Tz_facet = fd.conditional(fd.gt(fd.avg(Tz), 0.0), Tz('+')*Tz('-') / fd.avg(Tz), 0.0)

# Define the well source term
#f = fd.Function(V).interpolate(10*fd.exp(-r**2/eps)/eps**0.5)
xw = Delta_x*Nx*0.5
yw = Delta_y*Ny*0.5

# line = np.array([xw, yw, 1])*np.ones([10,3])
# line[:, 2] = np.linspace(0., Lz, len(line)) 

radius = 0.1 #0.1875*0.3048 # 0.1 radius of well # 0.1875*0.3048 from ChenZhang2009
f = fd.Function(V)
f.assign(fd.interpolate(fd.conditional(pow(x-xw,2)+pow(y-yw,2)<pow(radius,2), fd.exp(-(1.0/(-pow(x-xw,2)-pow(y-yw,2)+pow(radius,2))))*Kx, 0.0), V))
norm = fd.assemble(f*fd.dx)

pv = fd.assemble(phi*fd.dx) # pore volume
q = -0.1*pv/(30*3600*24) # prod rate of 0.01*pv / month
f.assign(q*f/norm)

# Plot source
# rr = np.logspace(-5, np.log(Delta_x*Nx),100)
# ff = np.exp(-rr**2/eps)/eps**0.5
# plt.loglog(rr, ff)
# plt.show()

# We can now define the bilinear and linear forms for the left and right
dx = fd.dx
TdivU = fd.as_vector((Tx_facet*u.dx(0), Ty_facet*u.dx(1), Tz_facet*u.dx(2)))
# F = (u - u_n)/dt*w*v*dx + (fd.dot(TdivU, fd.grad(v)))*dx + f*v*dx
F = w*u*v*dx + dt*(fd.dot(TdivU, fd.grad(v)))*dx -(w*u_n + dt*f)*v*dx
a = fd.lhs(F)
L = fd.rhs(F)

# a = inner(u,v)*dx+ dt*(fd.dot(TdivU, fd.grad(v)))*dx
# L = (u_n + dt*f)*v*dx
# m = u * v * w * dx

u = fd.Function(V)
print('Start solver')

t = []
pwf = []
pavg = []
q = []

# t =0
t.append(0.)
pwf.append(u_n.at([xw, yw, 0]))
pavg.append(fd.assemble(phi*u_n*fd.dx)/pv)
q.append(0)

print("t [h] = {:6.3f}  \t pavg [bar] = {:6.3f}".format(t[-1]/3600, pavg[-1]/1e5))
shutin = False
for n in range(len(dts)):
    fd.solve(a == L, u, solver_parameters={'ksp_type': 'cg'})
    
    # Update previous solution
    u_n.assign(u)
    dt.assign(fd.Constant(dts[n]))
    
    # Save t, pwf and pavg to list
    t.append(t[-1] + dts[n])
    pwf.append(u_n.at([xw, yw, 0]))
    pavg.append(fd.assemble(u_n*phi*dx)/pv)
    q.append((pavg[-1]-pavg[-2])*ct*pv/dts[n])
    
    # Update source term 
    if t[-1]>=tp*3600 and not shutin:
        f.assign(fd.Constant(0))
        shutin=True
        print('well shutin')

    # Save to VTK file
    if n % 10 == 0:
        outfile.write(u_n)

    # print info 
    print("t [h] = {:6.3f}  \t pavg [bar] = {:6.3f} \t pwf = {:6.3f} \t q [m3/d] = {:4.1f}".format(t[-1]/3600, pavg[-1]/1e5, pwf[-1]/1e5, q[-1]*24*3600))

    # defining dataframe for results
    df = pd.DataFrame(columns=['t (h)', 'pwf (bar)','pavg (bar)', 'q (m3/d)'])
    df['t (h)'] = np.array(t)/3600
    df['pwf (bar)'] = np.array(pwf)/1e5
    df['pavg (bar)'] = np.array(pavg)/1e5
    df['q (m3/d)'] = np.array(q)*24*3600
    df.to_csv('solution.dat')

    if abs(pwf[-1]-pavg[-1]) < 1.:
        break
