
import firedrake as fd
import numpy as np
from slepc4py import SLEPc
from firedrake.petsc import PETSc
import scipy.io as sp


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
Delta_z = 2*0.3048

Nx = int(60)
Ny = int(220)

layers = range(35,46,1)
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

to_days = 24*3600

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

fd.File("spe10_small.pvd").write(Kx, Ky, Kz, phi)

# Permeability field harmonic interpolation to facets
Tx_facet = fd.conditional(fd.gt(fd.avg(Tx), 0.0), Tx('+')*Tx('-') / fd.avg(Tx), 0.0)
Ty_facet = fd.conditional(fd.gt(fd.avg(Ty), 0.0), Ty('+')*Ty('-') / fd.avg(Ty), 0.0)
Tz_facet = fd.conditional(fd.gt(fd.avg(Tz), 0.0), Tz('+')*Tz('-') / fd.avg(Tz), 0.0)

# We can now define the bilinear and linear forms for the left and right
dx = fd.dx
TdivU = fd.as_vector((Tx_facet*u.dx(0), Ty_facet*u.dx(1), Tz_facet*u.dx(2)))
a = (fd.dot(TdivU, fd.grad(v))) * dx
m = u * v * w * dx

# Defining the eigenvalue problem

petsc_a = fd.assemble(a).M.handle
petsc_m = fd.assemble(m).M.handle

# Save *.h5 file
ViewHDF5 = PETSc.Viewer()     # Init. Viewer
ViewHDF5.createHDF5('A.h5', mode=PETSc.Viewer.Mode.WRITE,comm= PETSc.COMM_WORLD)
ViewHDF5(petsc_a)   # Put PETSc object into the viewer
ViewHDF5.destroy()            # Destroy Viewer

ViewHDF5 = PETSc.Viewer()     # Init. Viewer
ViewHDF5.createHDF5('M.h5', mode=PETSc.Viewer.Mode.WRITE,comm= PETSc.COMM_WORLD)
ViewHDF5(petsc_m)   # Put PETSc object into the viewer
ViewHDF5.destroy()            # Destroy Viewer

# Load *.h5 file
#ViewHDF5 = PETSc.Viewer()     # Init. Viewer
#ViewHDF5.createHDF5('grid.h5', mode=PETSc.Viewer.Mode.READ,comm= PETSc.COMM_WORLD)
#ViewHDF5.view(obj=petsc_a)   # Put PETSc object into the viewer
#ViewHDF5.destroy()            # Destroy Viewer

# Set solver options
num_eigenvalues = 5

opts = PETSc.Options()
opts.setValue("eps_gen_hermitian", None)
opts.setValue("eps_monitor", None)

opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue('st_ksp_type', 'preonly')
opts.setValue('st_pc_type', 'lu')
opts.setValue('st_type', 'sinvert')
opts.setValue("eps_type", "krylovschur")
opts.setValue('pc_factor_mat_solver_type', 'mumps')

#opts.setValue('st_ksp_type', 'gmres')
#opts.setValue('st_pc_type', 'bjacobi')

opts.setValue("eps_target_magnitude", None)
opts.setValue("eps_target", 1e-6)
opts.setValue("eps_tol", 1e-6)
opts.setValue("st_ksp_max_it", 200)
opts.setValue("st_ksp_rtol", 1e-6)

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
    eigvalues.append(lam.real)
    
    eigvecs.append(fd.Function(V))
    eigvecs[-1].vector()[:] = vr
    eigvecs[-1].rename('eigvec'+str('{:2d}').format(i))


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
np.savetxt('eigvalues.txt',np.array(eigvalues))
fd.File("spe10_small.pvd").write(phi, Kx, Ky, Kz, *eigvecs)