    label Nparcels = dieselSpray.size();
    reduce(Nparcels, sumOp<label>());

    Info<< "\nNumber of parcels in system.... | "
        << Nparcels << endl
        << "Injected liquid mass........... | "
        << 1e6*dieselSpray.injectedMass(runTime.value()) << " mg" << endl
        << "Liquid Mass in system.......... | "
        << 1e6*dieselSpray.liquidMass() << " mg" << endl
        << "SMD, Dmax...................... | "
        << dieselSpray.smd()*1e6 << " mu, "
        << dieselSpray.maxD()*1e6 << " mu"
        << endl;
