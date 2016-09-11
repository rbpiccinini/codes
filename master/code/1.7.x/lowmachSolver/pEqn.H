//updates rho after solving for enthalpy
rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();

//Computes U without the pressure contribution
U = rUA*UEqn.H();
//Updates phi
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        );

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
	fvScalarMatrix pdEqn
	(
		fvc::ddt(psi,p)
      + fvc::div(phi)
	  - fvm::laplacian(rho*rUA, pd)
	 ==
		Sevap
	);

	pdEqn.solve();

	if (nonOrth == nNonOrthCorr)
	{
		phi += pdEqn.flux();
	}
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rUA*fvc::grad(pd);
U.correctBoundaryConditions();
