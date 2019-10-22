
import firedrake as fd
import numpy as np
from slepc4py import SLEPc
from firedrake.petsc import PETSc
import scipy.io as sp
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
T = 20.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size

# Define mesh
Delta_x = 20*0.3048 # ft -> m
Delta_y = 10*0.3048
Delta_z = 2*0.3048

Nx = int(60)
Ny = int(220)

layers = range(35,40)
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
                                    coords[:, 2], Delta=Delta, data_array=kx_array/mu*24*3600)
Ty.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=ky_array/mu*24*3600)
Tz.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kz_array/mu*24*3600)

w.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                       coords[:, 2], Delta=Delta, data_array=(phi_array+1e-4)*ct)
print("END: Read in reservoir fields")

# Define initial value
x, y, z = fd.SpatialCoordinate(mesh)
eps = 1.
xc = 0.5*Delta_x*Lx
yc = 0.5*Delta_y*Ly
r = ((x-xc)**2+(y-yc)**2)**0.5 

# u_0 = fd.Function(V).interpolate(10*y/(Delta_y*Ny))
u_0 = fd.Function(V).interpolate(fd.Constant(0.))
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
radius = 0.1 #0.1875*0.3048 # 0.1 radius of well # 0.1875*0.3048 from ChenZhang2009
f = fd.Function(V)
f.assign(fd.interpolate(fd.conditional(pow(x-xw,2)+pow(y-yw,2)<pow(radius,2), fd.exp(-(1.0/(-pow(x-xw,2)-pow(y-yw,2)+pow(radius,2)))), 0.0), V))
norm = fd.assemble(f*fd.dx)
f.assign(1e-4*f/norm)

# Plot source
rr = np.logspace(-5, np.log(Delta_x*Nx),100)
ff = np.exp(-rr**2/eps)/eps**0.5
plt.loglog(rr, ff)
# plt.show()

# We can now define the bilinear and linear forms for the left and right
dx = fd.dx
TdivU = fd.as_vector((Tx_facet*u.dx(0), Ty_facet*u.dx(1), Tz_facet*u.dx(2)))
# F = (u - u_n)/dt*w*v*dx + (fd.dot(TdivU, fd.grad(v)))*dx + f*v*dx
F = u*v*dx + dt/w*(fd.dot(TdivU, fd.grad(v)))*dx -(u_n + dt/w*f)*v*dx
a = fd.lhs(F)
L = fd.rhs(F)

# a = inner(u,v)*dx+ dt*(fd.dot(TdivU, fd.grad(v)))*dx
# L = (u_n + dt*f)*v*dx
# m = u * v * w * dx

u = fd.Function(V)
t = 0
print('Start solver')
pavg0 = fd.assemble(u_0*dx)/vol
print("t= {:6.3f}  \t pavg = {:1.4e}".format(t, pavg0))

pwf = []
pwf.append([t, u_0.at([xw, yw, 0])])
for n in range(num_steps):
    fd.solve(a == L, u, solver_parameters={'ksp_type': 'cg'})
    
    # Update previous solution
    u_n.assign(u)
    t += dt
    
    # Save to file and plot solution
    # if n % 20 == 0:
    outfile.write(u_n)
    pwf.append([t, u_n.at([xw, yw, 0])])
    pavg = fd.assemble(u_n*dx)/vol
    print("t= {:6.3f}  \t pavg = {:1.4e}".format(t, pavg))


pwf = np.array(pwf)
