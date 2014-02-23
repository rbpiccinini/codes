    #include "readGravitationalAcceleration.H"
        
    
    Info<< "Reading wettingPhaseProperties\n" << endl;

    IOdictionary wettingPhaseProperties
    (
        IOobject
        (
            "wettingPhaseProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    dimensionedScalar const_muw(wettingPhaseProperties.lookup("muw"));    
    volScalarField muw 
    (
        IOobject
        (
            "muw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        const_muw
    );
    
    dimensionedScalar const_rhow(wettingPhaseProperties.lookup("rhow"));  
    volScalarField rhow 
    (
        IOobject
        (
            "rhow",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        const_rhow
    );

    
    Info<< "Reading nonWettingPhaseProperties\n" << endl;

    IOdictionary nonWettingPhaseProperties
    (
        IOobject
        (
            "nonWettingPhaseProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar const_mun(nonWettingPhaseProperties.lookup("mun"));
    volScalarField mun 
    (
        IOobject
        (
            "mun",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        const_mun
    );
    
    
    dimensionedScalar const_rhon(nonWettingPhaseProperties.lookup("rhon"));
    volScalarField rhon 
    (
        IOobject
        (
            "rhon",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        const_rhon
    );

    
    
    
    IOdictionary rockProperties
    (
        IOobject
        (
            "rockProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar cr
    (
        rockProperties.lookup("cr")
    );

    Info<< "Reading field p\n" << endl;
    
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field K\n" << endl;
    volTensorField K
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    Info<< "Reading field Eps\n" << endl;
    volScalarField Eps
    (
        IOobject
        (
            "Eps",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    Info<< "Reading field Sw\n" << endl;
    volScalarField Sw
    (
        IOobject
        (
            "Sw",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField Sn
    (
        IOobject
        (
            "Sn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.-Sw
    );
    
    IOdictionary rockFluidProperties
    (
        IOobject
        (
            "rockFluidProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    


    
    Info<< "Reading/calculating face flux field phi\n" << endl;
    
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );
    
  
    surfaceScalarField Swf
    (
        IOobject
        (
            "Swf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(Sw)
    );
    
//    surfaceScalarField Swf = fvc::interpolate(Sw);
    surfaceScalarField rhowf = fvc::interpolate(rhow);
    surfaceScalarField rhonf = fvc::interpolate(rhon);
    surfaceScalarField muwf = fvc::interpolate(muw);
    surfaceScalarField munf = fvc::interpolate(mun);
    
    autoPtr<rockFluidModel> rockFluid = rockFluidModel::New
    (
        "rock",
        rockFluidProperties,
        Swf,
        rhowf,
        rhonf,
        muwf,
        munf
    );
    
    rockFluid->correct();
    
/*
    volScalarField wKrw
    (
        IOobject
        (
            "wKrw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/
//    krw = rockFluid->krw(); ??? doesnt work
    
    Info<< "Computing fields M and L\n" << endl;
    
    surfaceTensorField Kf = fvc::interpolate(K,"harmonic");
    surfaceTensorField M =  Kf*rockFluid->M();
    surfaceTensorField L =  Kf*rockFluid->L();

    surfaceScalarField Fw
    (
        IOobject
        (
            "Fw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rockFluid->Fw()
    );

    surfaceScalarField phiGG
    (
        IOobject
        (
            "phiGG",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (L & g) & mesh.Sf()
    );
    
    surfaceScalarField phiS
    (
        IOobject
        (
            "phiS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*rockFluid->Fw()
    );

    label pRefCell = 0;
    scalar pRefValue = 0.;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
