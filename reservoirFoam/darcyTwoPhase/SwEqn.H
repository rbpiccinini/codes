        
        phiS = fvc::flux(Fw*phi,Sw,"div(phiS)");
        surfaceScalarField phiSB = fvc::flux(Fw*phi,Sw,"upwind");
        
        for
        (
            subCycle<volScalarField> SwSubCycle(Sw, nSubCycles);
            !(++SwSubCycle).end();
        )
        {
            phiS = fvc::flux(Fw*phi,Sw,"div(phiS)");
            
            MULES::explicitSolve
            (
                Sw,
                phiSB,
                phiS,
                1.0,
                0.0
            );
            
            
            Swf = fvc::interpolate(Sw);
            rockFluid->correct();
                   
            Fw=rockFluid->Fw();
        }
        
        Sn=1.-Sw;
