import firedrake as fd
import numpy as np
from slepc4py import SLEPc
from firedrake.petsc import PETSc


def coords2ijk(x, y, Delta, data_array):
    i = np.floor(x / Delta[0]).astype(int)
    j = np.floor(y / Delta[1]).astype(int)
    if i == 60:
        i = i-1
    if j == 220:
        j = i-1

    return data_array[i, j]

# Define mesh
Delta_x = 20 * 0.3048 # ft -> m
Delta_y = 10 * 0.3048 # ft -> m
Delta_z =  2 * 0.3048 # ft -> m

Nx = 60
Ny = 220
Nz = 1 #85

Lx = Nx*Delta_x
Ly = Ny*Delta_y
Lz = Nz*Delta_z

mesh = fd.utility_meshes.RectangleMesh(nx=Nx, ny=Ny, Lx=Lx, Ly=Ly,
                                       quadrilateral=True)
# mesh = fd.ExtrudedMesh(mesh, layers=Nz, layer_height=Delta_z)

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
Delta = np.array([Delta_x, Delta_y])

# compressibility
ct = 14.54E-11 # 1/Pa

# viscosity
mu = 1e-3# Pa

coords = fd.project(mesh.coordinates, Vvec).dat.data
kx_array = np.load('./spe10/spe10z85_kx.npy') * 1e-15/ct/mu # mD -> m^2
ky_array = np.load('./spe10/spe10z85_ky.npy') * 1e-15/ct/mu # mD -> m^2
phi_array = np.load('./spe10/spe10z85_phi.npy')



#kx_array = np.load('./spe10/spe10_kx.npy')
#ky_array = np.load('./spe10/spe10_ky.npy')
#phi_array = np.load('./spe10/spe10_po.npy')
#
#kx_array = kx_array[:,:,-1]
#ky_array = ky_array[:,:,-1]
#phi_array = phi_array[:,:,-1]
#
#np.save('./spe10z85_kx.npy', kx_array)
#np.save('./spe10z85_ky.npy', ky_array)
#np.save('./spe10z85_phi.npy', phi_array)
kx_geomean = np.exp(np.log(kx_array).sum()/(Nx*Ny))
ky_geomean = np.exp(np.log(ky_array).sum()/(Nx*Ny))
Kx_hom = fd.Constant(kx_geomean)
Ky_hom = fd.Constant(ky_geomean)
phi_hom = fd.Constant(phi_array.mean())
print("END: Read in reservoir fields")
# -----------------------------------------------------------------------------
# Solve homogeneous problem
# -----------------------------------------------------------------------------
print("START: Solve homogeneous problem")
dx = fd.dx
KdivU = fd.as_vector((Kx_hom*u.dx(0), Ky_hom*u.dx(1)))
a = (fd.dot(KdivU, fd.grad(v))) * dx
m = u * v * phi_hom * dx

# Defining the eigenvalue problem

petsc_a = fd.assemble(a).M.handle
petsc_m = fd.assemble(m).M.handle

# save files
viewer = PETSc.Viewer().createBinary('A_small.dat', 'w')
viewer.pushFormat(viewer.Format.NATIVE)
viewer(petsc_a)
viewer = PETSc.Viewer().createBinary('M_small.dat', 'w')
viewer.pushFormat(viewer.Format.NATIVE)
viewer(petsc_m)

viewer_new = PETSc.Viewer().createBinary('A_2d.dat', 'r')
viewer_new.pushFormat(viewer.getFormat())
petsc_a_new = PETSc.Mat().load(viewer_new)
viewer_new = PETSc.Viewer().createBinary('M_2d.dat', 'r')
viewer_new.pushFormat(viewer.getFormat())
petsc_m_new = PETSc.Mat().load(viewer_new)
assert petsc_a_new.equal(petsc_a), "Reload unsuccessful"
assert petsc_m_new.equal(petsc_m), "Reload unsuccessful"
print("Reload successful")



num_eigenvalues = 30

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
es = SLEPc.EPS().create()
es.setDimensions(num_eigenvalues)
es.setOperators(petsc_a, petsc_m)
es.setFromOptions()
es.solve()

# Number of converged eigenvalues
nconv = es.getConverged()
eigvecs = []
eigvalues_hom = []
lamvec_hom = [fd.Function(V)]
lamvec_hom[0].vector()[:] = 0
for i in range(nconv):
    print(es.getEigenvalue(i).real)
    eigvalues_hom.append(es.getEigenvalue(i).real)
    
    vr, vi = petsc_a.getVecs()
    lam = es.getEigenpair(i, vr, vi)
    
    eigvecs.append(fd.Function(V))
    eigvecs[-1].vector()[:] = vr
    eigvecs[-1].rename('eigvec'+str('{:2d}').format(i))
    
    lamvec_hom.append(fd.Function(V))
    lamvec_hom[-1].vector()[:] = lamvec_hom[-2].vector()[:] + vr*lam.real*np.exp(-lam.real*3600*24)
    lamvec_hom[-1].rename('lamvec'+str('{:2d}').format(i))
    
eigvalues_hom = np.array(eigvalues_hom)*3600*24 # 1/s -> 1/day

Kx_hom = fd.Function(V).interpolate(Kx_hom)
Ky_hom = fd.Function(V).interpolate(Ky_hom)
phi_hom = fd.Function(V).interpolate(phi_hom)
Kx_hom.rename('Kx')
Ky_hom.rename('Ky')
phi_hom.rename('phi')

outfile = fd.File("./results/homogeneous.pvd")
eigvec = fd.Function(V)
lamvec = fd.Function(V)
for i in range(len(eigvecs)):
    eigvec = eigvecs[i]
    lamvec = lamvec_hom[i]
    eigvec.rename('eigvec')
    lamvec.rename('lamvec')
    outfile.write(phi_hom, Kx_hom, Ky_hom, lamvec, eigvec, time=i)

# -----------------------------------------------------------------------------
# Solve heterogenous problem
# -----------------------------------------------------------------------------
print("START: Solve heterogenous problem")
Kx = fd.Function(V)
Ky = fd.Function(V)
phi = fd.Function(V)

Kx.rename('Kx')
Ky.rename('Ky')
phi.rename('phi')

coords2ijk = np.vectorize(coords2ijk, excluded=['data_array', 'Delta'])

Kx.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    Delta=Delta, data_array=kx_array)
Ky.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    Delta=Delta, data_array=ky_array)
phi.dat.data[...] = coords2ijk(coords[:, 0], coords[:, 1],
                                    Delta=Delta, data_array=phi_array)

fd.File("grid.pvd").write(Kx, Ky, phi)

# Permeability field harmonic interpolation to facets
Kx_facet = fd.conditional(fd.gt(fd.avg(Kx), 0.0), Kx('+')*Kx('-') / fd.avg(Kx), 0.0)
Ky_facet = fd.conditional(fd.gt(fd.avg(Ky), 0.0), Ky('+')*Ky('-') / fd.avg(Ky), 0.0)

# We can now define the bilinear and linear forms for the left and right
dx = fd.dx
KdivU = fd.as_vector((Kx_facet*u.dx(0), Ky_facet*u.dx(1)))
a = (fd.dot(KdivU, fd.grad(v))) * dx
m = u * v * phi* dx

# Defining the eigenvalue problem

petsc_a = fd.assemble(a).M.handle
petsc_m = fd.assemble(m).M.handle

num_eigenvalues = 30

# Set solver options

#opts = PETSc.Options()
#opts.setValue("eps_gen_hermitian", None)
#opts.setValue("st_pc_factor_shift_type", "NONZERO")
#opts.setValue("eps_type", "krylovschur")
#opts.setValue("eps_smallest_magnitude", None)
#opts.setValue("eps_tol", 1e-10)
#opts.setValue("eps_ncv", 100)
#opts.setValue("eps_max_it", 120)

opts = PETSc.Options()
opts.setValue("eps_gen_hermitian", None)
opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue('sp_initial_guess_nonzero', True)
opts.setValue('st_type', 'sinvert')
opts.setValue("eps_type", "krylovschur")
opts.setValue("eps_target_magnitude", None)
opts.setValue("eps_target", 0)
opts.setValue("eps_tol", 1e-20)

# Solve for eigenvalues
print('Computing eigenvalues...')
es = SLEPc.EPS().create().create(comm=SLEPc.COMM_WORLD)
es.setDimensions(num_eigenvalues)
es.setOperators(petsc_a, petsc_m)
es.setFromOptions()
nev, ncv, mpd = es.getDimensions()

es.solve()

# Number of converged eigenvalues
nconv = es.getConverged()
eigvecs = []
eigvalues_het = []
lamvec_het = [fd.Function(V)]
lamvec_het[0].vector()[:] = 0
for i in range(min(num_eigenvalues, nconv)):
    print(es.getEigenvalue(i).real)
    eigvalues_het.append(es.getEigenvalue(i).real)
    
    vr, vi = petsc_a.getVecs()
    lam = es.getEigenpair(i, vr, vi)
    
    eigvecs.append(fd.Function(V))
    eigvecs[-1].vector()[:] = vr
    eigvecs[-1].rename('eigvec'+str('{:2d}').format(i))
    
    lamvec_het.append(fd.Function(V))
    lamvec_het[-1].vector()[:] = lamvec_het[-2].vector()[:] + vr*lam.real*np.exp(-lam.real*3600*24)
    lamvec_het[-1].rename('lamvec'+str('{:2d}').format(i))

eigvalues_het = np.array(eigvalues_het)*3600*24 # 1/s -> 1/day
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

outfile = fd.File("./results/heterogeneous.pvd")
eigvec = fd.Function(V)
lamvec = fd.Function(V)
for i in range(len(eigvecs)):
    eigvec = eigvecs[i]
    lamvec = lamvec_het[i]
    eigvec.rename('eigvec')
    lamvec.rename('lamvec')
    outfile.write(phi, Kx, Ky, lamvec, eigvec, time=i)

plt.plot(eigvalues_het[1:20]/eigvalues_hom[1:20],':s')
