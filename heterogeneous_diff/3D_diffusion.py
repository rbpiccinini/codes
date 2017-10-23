from dolfin import *
import numpy as np
import time

# Read mesh from file and create function space
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(20.0, 10.0, 1.0), 30, 20, 4)
V = FunctionSpace(mesh, "Lagrange", 1)

# Code for C++ evaluation of conductivity
conductivity_code = """

class Conductivity : public Expression
{
public:

  // Create expression with 3 components
  Conductivity() : Expression(3) {}

  // Function for evaluating expression on each cell
  void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.index;
    values[0] = (*kx)[cell_index];
    values[1] = (*ky)[cell_index];
    values[2] = (*kx)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<MeshFunction<double> > kx;
  std::shared_ptr<MeshFunction<double> > ky;
  std::shared_ptr<MeshFunction<double> > kz;

};
"""

# Create mesh functions for c00, c01, c11
kx = MeshFunction("double", mesh, 3)
ky = MeshFunction("double", mesh, 3)
kz = MeshFunction("double", mesh, 3)

# Iterate over mesh and set values
for cell in cells(mesh):
    kx[cell] = 1.0
    ky[cell] = 1.0
    kz[cell] = 0.1 

c = Expression(cppcode=conductivity_code, degree=0)
c.kx = kx
c.ky = ky
c.kz = kz
C = as_matrix(((c[0], 0.0, 0.0), (0.0, c[1], 0.0), (0.0, 0.0, c[2])))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
a = inner(C*grad(u), grad(v))*dx
L = Constant(0.0)*v*dx
m = u*v*dx

# Assemble stiffness form
A, _ = assemble_system(a, L)
B = assemble(m)

# Create eigensolver
eigensolver = SLEPcEigenSolver(as_backend_type(A), as_backend_type(B))
eigensolver.parameters['spectrum']='smallest magnitude'
eigensolver.parameters['tolerance']=1E-10

# Compute all eigenvalues of A x = \lambda x
N=20
print("Computing eigenvalues with SLEPc. This can take a minute.")
start_time = time.time()
eigensolver.solve(N)
print('CPU time: %s seconds' % (time.time()-start_time))

eig = Function(V)
eig_vec = eig.vector()
eigvals=[]
print('converged: %d', eigensolver.get_number_converged())
for j in range(eigensolver.get_number_converged()):
    r, c, rx, cx = eigensolver.get_eigenpair(j)
    eigvals.append(r)
    eig_vec[:] = rx
    print(r, c)

np.savetxt('eigs.txt',eigvals)
