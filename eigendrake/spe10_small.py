
import firedrake as fd
import numpy as np
from slepc4py import SLEPc
from firedrake.petsc import PETSc


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
Delta_x = 20
Delta_y = 10
Delta_z = 2

Nx = int(60/10)
Ny = int(220/10)
Nz = int(85/17)

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

coords = fd.project(mesh.coordinates, Vvec).dat.data
kx_array = np.load('./spe10/spe10_kx.npy')[:Nx, :Ny, :Nz]
ky_array = np.load('./spe10/spe10_ky.npy')[:Nx, :Ny, :Nz]
kz_array = np.load('./spe10/spe10_kz.npy')[:Nx, :Ny, :Nz]
phi_array = np.load('./spe10/spe10_po.npy')[:Nx, :Ny, :Nz]

Kx = fd.Function(V)
Ky = fd.Function(V)
Kz = fd.Function(V)
phi = fd.Function(V)

Kx.rename('Kx')
Ky.rename('Ky')
Kz.rename('Kz')
phi.rename('phi')

coords2ijk = np.vectorize(coords2ijk, excluded=['data_array', 'Delta'])

Kx.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kx_array)
Ky.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=ky_array)
Kz.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    coords[:, 2], Delta=Delta, data_array=kz_array)
phi.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                       coords[:, 2], Delta=Delta, data_array=phi_array)
print("END: Read in reservoir fields")

fd.File("spe10_small.pvd").write(Kx, Ky, Kz, phi)

# Permeability field harmonic interpolation to facets
Kx_facet = fd.conditional(fd.gt(fd.avg(Kx), 0.0), Kx('+')*Kx('-') / fd.avg(Kx), 0.0)
Ky_facet = fd.conditional(fd.gt(fd.avg(Ky), 0.0), Ky('+')*Ky('-') / fd.avg(Ky), 0.0)
Kz_facet = fd.conditional(fd.gt(fd.avg(Kz), 0.0), Kz('+')*Kz('-') / fd.avg(Kz), 0.0)

# We can now define the bilinear and linear forms for the left and right
dx = fd.dx
KdivU = fd.as_vector((Kx_facet*u.dx(0), Ky_facet*u.dx(1), Kz_facet*u.dx(2)))
a = (fd.dot(KdivU, fd.grad(v))) * dx
m = u * v * phi * dx

# Defining the eigenvalue problem

petsc_a = fd.assemble(a).M.handle
petsc_m = fd.assemble(m).M.handle

num_eigenvalues = 10

# Set solver options
opts = PETSc.Options()
opts.setValue("eps_gen_hermitian", None)
opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue('st_type', 'sinvert')
opts.setValue("eps_type", "krylovschur")
opts.setValue("eps_target_magnitude", None)
opts.setValue("eps_target", 0)
opts.setValue("eps_tol", 1e-10)


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
for i in range(nconv):
    print(es.getEigenvalue(i).real)
    
    vr, vi = petsc_a.getVecs()
    lam = es.getEigenpair(i, vr, vi)
    
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

fd.File("helmholtz.pvd").write(Kx, Ky, Kz, *eigvecs)
