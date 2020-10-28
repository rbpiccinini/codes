
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

# Define mesh
Delta_x = 20*0.3048 # ft -> m
Delta_y = 10*0.3048
Delta_z = -2*0.3048

Nx = int(60)
Ny = int(220)

layers = range(1) # range(35,46,1)
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

# Create extruded mesh element
horiz_elt = fd.FiniteElement("DG", fd.quadrilateral, 0)
vert_elt =  fd.FiniteElement("DG", fd.interval, 0)
elt = fd.TensorProductElement(horiz_elt, vert_elt)

# Create function space
V = fd.FunctionSpace(mesh, elt)
Vvec = fd.VectorFunctionSpace(mesh, "DG", 0)

# We'll also need the test and trial functions corresponding to this
# function space::

u = fd.TrialFunction(V)
v = fd.TestFunction(V)

# We declare a function over our function space and give it the
# value of our right hand side function::

# Permeability tensor
print("START: Read in reservoir fields")
Delta = np.array([Delta_x, Delta_y, Delta_z])
ct = 79.08e-6/9.80665e4 # (1.0+0.2*3+0.8*4.947)*14.2 * 10**-6 kgf/cm2
mu = 0.003 # Pa-s

# coords = fd.project(mesh.coordinates, C).dat.data
coords = fd.project(mesh.coordinates, Vvec).dat.data
spe10_2 = sp.loadmat('../../input/spe10/spe10_2.mat')

# well location
xw = np.array([118.872, 181.356, Lz/2])

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

I = fd.Function(V)
J = fd.Function(V)
K = fd.Function(V)

I.rename('I')
J.rename('J')
K.rename('K')

Kx.rename('Kx')
Ky.rename('Ky')
Kz.rename('Kz')
phi.rename('phi')

Tx.rename('Kx')
Ty.rename('Ky')
Tz.rename('Kz')
w.rename('w')

coords2ijk = np.vectorize(coords2ijk, excluded=['data_array', 'Delta'])

to_days = 24*3600

I.dat.data[...] = np.floor(coords[:, 0] / Delta[0]).astype(int)
J.dat.data[...] = np.floor(coords[:, 1] / Delta[1]).astype(int)
K.dat.data[...] = np.floor(coords[:, 2] / Delta[2]).astype(int)          

Kx.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                   coords[:, 2], Delta=Delta, data_array=kx_array)
Ky.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                   coords[:, 2], Delta=Delta, data_array=ky_array)
Kz.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                   coords[:, 2], Delta=Delta, data_array=kz_array)
phi.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                      coords[:, 2], Delta=Delta, data_array=phi_array)
Tx.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                   coords[:, 2], Delta=Delta, data_array=kx_array/mu*to_days)
Ty.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                   coords[:, 2], Delta=Delta, data_array=ky_array/mu*to_days)
Tz.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                   coords[:, 2], Delta=Delta, data_array=kz_array/mu*to_days)
w.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                      coords[:, 2], Delta=Delta, data_array=phi_array*ct)

print("END: Read in reservoir fields")
fd.File("../../output/spatial_heterog.pvd").write(Kx, Ky, Kz, phi)

# Permeability field harmonic interpolation to facets
n = fd.FacetNormal(mesh)

Tx_facet = fd.conditional(fd.gt(fd.avg(Tx), 0.0), Tx('+')*Tx('-') / fd.avg(Tx), 0.0)
Ty_facet = fd.conditional(fd.gt(fd.avg(Ty), 0.0), Ty('+')*Ty('-') / fd.avg(Ty), 0.0)
Tz_facet = fd.conditional(fd.gt(fd.avg(Tz), 0.0), Tz('+')*Tz('-') / fd.avg(Tz), 0.0)

T_facet = (Tx_facet*(abs(n[0]('+'))+abs(n[0]('-')))/2 +
             Ty_facet*(abs(n[1]('+'))+abs(n[1]('-')))/2 +
             Tz_facet*(abs(n[2]('+'))+abs(n[2]('-')))/2)


# We can now define the bilinear and linear forms for the left and right
x,y,z = mesh.coordinates

x_func = fd.interpolate(x, V)
y_func = fd.interpolate(y, V)
z_func = fd.interpolate(z, V)

Delta_x = fd.jump(x_func)
Delta_y = fd.jump(y_func)
Delta_z = fd.jump(z_func)

Delta = fd.sqrt(Delta_x**2+Delta_y**2+Delta_z**2)

dx = fd.dx
# a = Txy_facet*fd.jump(v,n)*fd.jump(u,n)/Delta_h*fd.dS_h + Tz_facet*fd.jump(v)*fd.jump(u)/Delta_z*fd.dS_v
# m = u * v * w * dx

alpha = 1
gamma = 1

m = u * v * w * dx
a = T_facet/Delta*fd.dot(fd.jump(v,n), fd.jump(u,n))*fd.dS_h \
  + T_facet/Delta*fd.dot(fd.jump(v,n), fd.jump(u,n))*fd.dS_v 

# Defining the eigenvalue problem
petsc_a = fd.assemble(a).M.handle
petsc_m = fd.assemble(m).M.handle

print('## Matrix assembled.')

## Save *.h5 file
# ViewHDF5 = PETSc.Viewer()     # Init. Viewer
# ViewHDF5.createHDF5('A.h5', mode=PETSc.Viewer.Mode.WRITE,comm= PETSc.COMM_WORLD)
# ViewHDF5(petsc_a)   # Put PETSc object into the viewer
# ViewHDF5.destroy()            # Destroy Viewer
# 
# ViewHDF5 = PETSc.Viewer()     # Init. Viewer
# ViewHDF5.createHDF5('M.h5', mode=PETSc.Viewer.Mode.WRITE,comm= PETSc.COMM_WORLD)
# ViewHDF5(petsc_m)   # Put PETSc object into the viewer
# ViewHDF5.destroy()            # Destroy Viewer

## Load *.h5 file
#ViewHDF5 = PETSc.Viewer()     # Init. Viewer
#ViewHDF5.createHDF5('grid.h5', mode=PETSc.Viewer.Mode.READ,comm= PETSc.COMM_WORLD)
#ViewHDF5.view(obj=petsc_a)   # Put PETSc object into the viewer
#ViewHDF5.destroy()            # Destroy Viewer

# Set solver options
num_eigenvalues = 1500 

opts = PETSc.Options()

opts.setValue("eps_gen_hermitian", None)
opts.setValue("eps_monitor", None)

opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue('st_ksp_type', 'preonly')
opts.setValue('st_pc_type', 'lu')
opts.setValue('st_type', 'sinvert')
opts.setValue("eps_type", "krylovschur")
opts.setValue('pc_factor_mat_solver_type', 'mumps')

# opts.setValue('st_ksp_type', 'gmres')
# opts.setValue('st_pc_type', 'bjacobi')

opts.setValue("eps_target_magnitude", None)
opts.setValue("eps_target", 0.)
opts.setValue("eps_tol", 1e-13)
opts.setValue("st_ksp_max_it", 10000)
opts.setValue("st_ksp_rtol", 1e-8)

# opts.setValue("eps_gen_hermitian", None)
# opts.setValue("st_pc_factor_shift_type", "NONZERO")
# opts.setValue("eps_type", "krylovschur")
# opts.setValue("eps_smallest_real", None)
# opts.setValue("eps_tol", 1e-10)
# opts.setValue("eps_ncv", 40)

# Solve for eigenvalues
print('Computing eigenvalues...')
es = SLEPc.EPS().create().create(comm=SLEPc.COMM_WORLD)
es.setDimensions(num_eigenvalues)
es.setOperators(petsc_a, petsc_m)
es.setFromOptions()
print(es.getDimensions())

es.solve()

# Number of converged eigenvalues
nconv = es.getConverged()
eigvecs = []
eigvalues = []

for i in range(nconv):
    print(es.getEigenvalue(i).real)
    vr, vi = petsc_a.getVecs()
    lam = es.getEigenpair(i, vr, vi)
    
    eigvecs.append(fd.Function(V))
    eigvecs[-1].vector()[:] = vr*vr
    eigvalues.append(lam.real)
    eigvecs[-1].rename('eigvec^2'+str('{:2d}').format(i))

eigvalues = np.array(eigvalues)
eigvecs = np.array([x.vector()[:] for x in eigvecs])

# Elimina autovalor nulo
eigvalues = eigvalues[1:]
eigvecs = eigvecs[1:]

psi2_lam = np.einsum('ij,i->j', eigvecs, 1.0/eigvalues)
IP_array = 1.0/psi2_lam
IP_array = IP_array/IP_array.max()

IP = fd.Function(V)
IP.dat.data[...] = IP_array
IP.rename('IP')

Print = PETSc.Sys.Print

Print()
Print("******************************")
Print("*** SLEPc Solution Results ***")
Print("******************************")
Print()

its = es.getIterationNumber()
Print("Number of iterations of the method: %d" % its)

eps_type = es.getType()
Print("Solution method: %s" % eps_type)

nev, ncv, mpd = es.getDimensions()
Print("Number of requested eigenvalues: %d" % nev)

tol, maxit = es.getTolerances()
Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))

# For more details on how to specify solver parameters, see the section
# of the manual on :doc:`solving PDEs <../solving-interface>`.
#
# Next, we might want to look at the result, so we output our solution
# to a file::
fd.File("../../output/spatial_IP.pvd").write(I, J, K, phi, Kx, Ky, Kz, IP)
