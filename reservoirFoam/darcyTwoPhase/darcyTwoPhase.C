/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "rockFluidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"        
        #include "CourantNo.H"

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                
                M =  Kf*rockFluid->M();
                L =  Kf*rockFluid->L();
                
//                phiGG = (L & g) & mesh.Sf();
                
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(-M, p)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();
                
                if (nonOrth == nNonOrthCorr)
                {
                    phi = pEqn.flux();
                    U = fvc::reconstruct(phi);
                    U.correctBoundaryConditions();

                }
                
                
                for (int scorr=0; scorr<sCorr; scorr++)
                {
//                  volVectorField Uw = - rockFluid->krw() / muw * K & (fvc::grad(p)-rhow*g);
//                    volVectorField Uw = - rockFluid->krw() / muw * K & fvc::grad(p);
//                    surfaceVectorField Uww =  - fvc::interpolate(Sw) / muw * linearInterpolate(K) & (fvc::snGrad(p)*mesh.Sf()/mesh.magSf());
//                  surfaceScalarField phiW = linearInterpolate(Uw) & mesh.Sf();
                    
//                    phi += linearInterpolate(rhow - rhon) * rockFluid->krn() / munf * (Kf & g) & mesh.Sf();
//                    phiS = Fw*phi;

/*                    
                    phiS = fvc::flux(phi,Sw,"div(phiS)");
                    
                    fvScalarMatrix SwEqn
                     (
                        Eps * fvm::ddt(Sw) 
                        + fvc::div(phiS)
                     );
                     SwEqn.relax();
                     SwEqn.solve();
                     
                     Sn=1.-Sw;
                     Swf = fvc::interpolate(Sw);
                     rockFluid->correct();
                     
                     Fw=rockFluid->Fw(); 
                     
                     
                    
*/
                    #include "SwEqn.H"
                    
                    

                
                }
                
            }

            #include "continuityErrs.H"
        }
      
                
        runTime.write();
        Info<< "Sw [min, max] = " << min(Sw).value() << ", " << max(Sw).value() << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
