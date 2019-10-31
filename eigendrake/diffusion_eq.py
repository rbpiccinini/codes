# Simple Diffusion equation
# =========================
#
# Let's start by considering the modified Helmholtz equation on a unit square,
# :math:`\Omega`, with boundary :math:`\Gamma`:
#
# .. math::
#
#    -\nabla^2 u &= f
#
#    \nabla u \cdot \vec{n} &= 0 \quad \textrm{on}\ \Gamma
#
# for some known function :math:`f`. The solution to this equation will
# be some function :math:`u\in V`, for some suitable function space
# :math:`V`, that satisfies these equations. 

from firedrake import *

T = 2.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size

# Create mesh and define function space
nx = ny = 30
mesh = RectangleMesh( nx, ny, 4, 4, quadrilateral=True)
V = FunctionSpace(mesh, "CG", 1)

# Define initial value
x, y = SpatialCoordinate(mesh)
u_0 = Constant(0) # Function(V).interpolate(exp(-5*pow(x-2, 2) - 5*pow(y-2, 2)))
u_n = Function(V).assign(u_0)

# We declare the output filename, and write out the initial condition.
outfile = File("diff_eq.pvd")
outfile.write(u_n)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
# f = Constant(1)

# f = Function(V)
x, y = SpatialCoordinate(mesh)
# f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))
eps=1e-4
r = sqrt((x-2)**2+(y-2)**2) 
f = Function(V).interpolate(eps/(r**2+eps**2))

# We can now define the bilinear and linear forms for the left and right
# hand sides of our equation respectively::
a =  inner(u, v)*dx + dt*dot(grad(u), grad(v))*dx
L = (u_n + dt*f)*v*dx

# Finally we solve the equation. We redefine `u` to be a function
# holding the solution:: 

u = Function(V)
t = 0
for n in range(num_steps):
    solve(a == L, u, solver_parameters={'ksp_type': 'cg'})
    
    # Update previous solution
    u_n.assign(u)
    
    # Save to file and plot solution
    # if n % 20 == 0:
    outfile.write(u_n)
    print("t= {:6.3f}".format(t))
    t += dt

# For more details on how to specify solver parameters, see the section
# of the manual on :doc:`solving PDEs <../solving-interface>`.
#
# Next, we might want to look at the result, so we output our solution
# to a file::

