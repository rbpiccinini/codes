Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<psiChemistryModel> pChemistry
(
    psiChemistryModel::New(mesh)
);
psiChemistryModel& chemistry = pChemistry();

hsCombustionThermo& thermo = chemistry.thermo();

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

if (!composition.contains(inertSpecie))
{
    FatalErrorIn(args.executable())
        << "Specified inert specie '" << inertSpecie << "' not found in "
        << "species list. Available species:" << composition.species()
        << exit(FatalError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
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


Info<< "Reading field pd\n" << endl;
volScalarField pd
(
    IOobject
    (
        "pd",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
volScalarField& hs = thermo.hs();

    
#include "compressibleCreatePhi.H"

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "chemistry::Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Interpolating droplet properties to grid cells
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volScalarField drop_smd
(
 IOobject
 (
  "drop_smd",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("drop_smd", dimensionSet(0, 1, 0, 0, 0), 0.0)
 );
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volVectorField class0_U
(
 IOobject
 (
  "class0_U",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class0_U", dimensionSet(0, 1, -1, 0, 0),vector::zero)
 );

volVectorField class0_mU
(
 IOobject
 (
  "class0_mU",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class0_mU", dimensionSet(1, 1, -1, 0, 0),vector::zero)
 );

volScalarField class0_T
(
 IOobject
 (
  "class0_T",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class0_T", dimensionSet(0, 0, 0, 1, 0), 0.0)
 );

volScalarField class0_m
(
 IOobject
 (
  "class0_m",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class0_m", dimensionSet(1, 0, 0, 0, 0), 0.0)
 );

volScalarField class0_N
(
 IOobject
 (
  "class0_N",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("zero", dimless, 0.0)
 );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volVectorField class1_U
(
 IOobject
 (
  "class1_U",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class1_U", dimensionSet(0, 1, -1, 0, 0),vector::zero)
 );

volVectorField class1_mU
(
 IOobject
 (
  "class1_mU",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class1_mU", dimensionSet(1, 1, -1, 0, 0),vector::zero)
 );

volScalarField class1_T
(
 IOobject
 (
  "class1_T",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class1_T", dimensionSet(0, 0, 0, 1, 0), 0.0)
 );

volScalarField class1_m
(
 IOobject
 (
  "class1_m",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class1_m", dimensionSet(1, 0, 0, 0, 0), 0.0)
 );

volScalarField class1_N
(
 IOobject
 (
  "class1_N",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("zero", dimless, 0.0)
 );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volVectorField class2_U
(
 IOobject
 (
  "class2_U",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class2_U", dimensionSet(0, 1, -1, 0, 0),vector::zero)
 );
 
volVectorField class2_mU
(
 IOobject
 (
  "class2_mU",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class2_mU", dimensionSet(1, 1, -1, 0, 0),vector::zero)
 );

volScalarField class2_T
(
 IOobject
 (
  "class2_T",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class2_T", dimensionSet(0, 0, 0, 1, 0), 0.0)
 );

volScalarField class2_m
(
 IOobject
 (
  "class2_m",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class2_m", dimensionSet(1, 0, 0, 0, 0), 0.0)
 );

volScalarField class2_N
(
 IOobject
 (
  "class2_N",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("zero", dimless, 0.0)
 );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volVectorField class3_U
(
 IOobject
 (
  "class3_U",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class3_U", dimensionSet(0, 1, -1, 0, 0),vector::zero)
 );

volVectorField class3_mU
(
 IOobject
 (
  "class3_mU",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("class3_mU", dimensionSet(1, 1, -1, 0, 0),vector::zero)
 );

volScalarField class3_T
(
 IOobject
 (
  "class3_T",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class3_T", dimensionSet(0, 0, 0, 1, 0), 0.0)
 );

volScalarField class3_m
(
 IOobject
 (
  "class3_m",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("class3_m", dimensionSet(1, 0, 0, 0, 0), 0.0)
 );

volScalarField class3_N
(
 IOobject
 (
  "class3_N",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("zero", dimless, 0.0)
 );
 
 volScalarField rDeltaT
    (
        IOobject
        (
            "rDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimensionSet(0, 0, -1, 0, 0), 1),
        zeroGradientFvPatchScalarField::typeName
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
