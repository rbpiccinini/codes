import firedrake as fd
import numpy as np
from slepc4py import SLEPc
from firedrake.petsc import PETSc

Lx = 1.
Ly = Lx/2.
Lz = Lx/3.
mesh = fd.utility_meshes.RectangleMesh(nx=30, ny=15, Lx=Lx, Ly=Ly,
                                       quadrilateral=True)
nz = 10
mesh = fd.ExtrudedMesh(mesh, layers=nz, layer_height=Lz/nz)

# We need to decide on the function space in which we'd like to solve the
# problem. Let's use piecewise linear functions continuous between
# elements::

V = fd.FunctionSpace(mesh, "CG", 1)
W = fd.TensorFunctionSpace(mesh, "CG", 1)

# We'll also need the test and trial functions corresponding to this
# function space::

u = fd.TrialFunction(V)
v = fd.TestFunction(V)

# We declare a function over our function space and give it the
# value of our right hand side function::


# Permeability tensor
# Homogeneous
one = fd.Constant(1.0)
zero = fd.Constant(0.0)
Kx = fd.interpolate(one, V)
Ky = fd.interpolate(one, V)
Kz = fd.interpolate(one, V)

# Heterogeneous
x, y, z = fd.SpatialCoordinate(W.mesh()) 

#Kx = fd.Function(V).interpolate(1.0+y)
#Ky = fd.Function(V).interpolate(1.0+0.25*fd.sin(fd.pi/Ly*y/2))
#Kz = fd.Function(V).interpolate(z)

Kx.rename('Kx', 'Permeability in x direction.')
Ky.rename('Ky', 'Permeability in y direction.')
Kz.rename('Kz', 'Permeability in z direction.')

#Kx = fd.Constant(1.0)
#Ky = fd.Constant(1.0)
#Kz = fd.Constant(1.0)

# Permeability field harmonic interpolation to facets
Kx_facet = fd.conditional(fd.gt(fd.avg(Kx), 0.0), Kx('+')*Kx('-') / fd.avg(Kx), 0.0)
Ky_facet = fd.conditional(fd.gt(fd.avg(Ky), 0.0), Ky('+')*Ky('-') / fd.avg(Ky), 0.0)
Kz_facet = fd.conditional(fd.gt(fd.avg(Kz), 0.0), Kz('+')*Kz('-') / fd.avg(Kz), 0.0)

# Porosity
por = fd.Constant(1.0)

# We can now define the bilinear and linear forms for the left and right
dx = fd.dx
KdivU = fd.as_vector((Kx_facet*u.dx(0), Ky_facet*u.dx(1), Kz_facet*u.dx(2)))
a = (fd.dot(KdivU, fd.grad(v))) * dx
m = u * v * por * dx

# Defining the eigenvalue problem

petsc_a = fd.assemble(a).M.handle
petsc_m = fd.assemble(m).M.handle
num_eigenvalues = 20

# Set solver options
opts = PETSc.Options()
opts.setValue("eps_gen_hermitian", None)
opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue("eps_type", "krylovschur")
opts.setValue("eps_smallest_real", None)
opts.setValue("eps_tol", 1e-10)
opts.setValue("eps_ncv", 40)


# Solve for eigenvalues
es = SLEPc.EPS().create()
es.setDimensions(num_eigenvalues)
es.setOperators(petsc_a, petsc_m)
es.setFromOptions()
es.solve()

# Number of converged eigenvalues
nconv = es.getConverged()
eigvecs = []
eigvalues = []
for i in range(nconv):
    print(es.getEigenvalue(i).real)
    eigvalues.append(es.getEigenvalue(i).real)
    
    vr, vi = petsc_a.getVecs()
    lam = es.getEigenpair(i, vr, vi)
    
    eigvecs.append(fd.Function(V))
    eigvecs[-1].vector()[:] = vr
    eigvecs[-1].rename('eigvec'+str('{:2d}').format(i))

Print = PETSc.Sys.Print

print()
print("******************************")
print("*** SLEPc Solution Results ***")
print("******************************")
print()

its = es.getIterationNumber()
print("Number of iterations of the method: %d" % its)

eps_type = es.getType()
print("Solution method: %s" % eps_type)

nev, ncv, mpd = es.getDimensions()
print('nev = {:d}, ncv = {:d}, mpd = {:d}'.format(nev, ncv, mpd))

tol, maxit = es.getTolerances()
print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))

# For more details on how to specify solver parameters, see the section
# of the manual on :doc:`solving PDEs <../solving-interface>`.
#
# Next, we might want to look at the result, so we output our solution
# to a file::

fd.File("helmholtz.pvd").write(Kx, Ky, Kz, *eigvecs)
