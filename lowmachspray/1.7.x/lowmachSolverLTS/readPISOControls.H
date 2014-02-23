    dictionary piso = mesh.solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr =
        piso.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        piso.lookupOrDefault<Switch>("momentumPredictor", true);

    bool transonic =
        piso.lookupOrDefault<Switch>("transonic", false);

    int nOuterCorr =
        piso.lookupOrDefault<int>("nOuterCorrectors", 1);

	scalar maxCo(readScalar(piso.lookup("maxCo")));

// Maximum time scale
	scalar maxDeltaT(piso.lookupOrDefault<scalar>("maxDeltaT", GREAT));

// Smoothing parameter (0-1) when smoothing iterations > 0
	scalar rDeltaTSmoothingCoeff
	(
		piso.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
	);

// Maximum change in cell temperature per iteration (relative to previous value)
	scalar alphaTemp(piso.lookupOrDefault("alphaTemp", 0.05));

Info<< "    maxCo..................." << maxCo << nl
    << "    maxDeltaT..............." << maxDeltaT << nl
    << "    alphaTemp..............." << alphaTemp << nl
    << "    rDeltaTSmoothingCoeff..." << rDeltaTSmoothingCoeff << nl
    << endl;
