    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimple.lookupOrDefault<int>("nOuterCorrectors", 1);
        
    const int nCorr =
        pimple.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimple.lookupOrDefault("momentumPredictor", true);
        
//    const bool transonic =
//        pimple.lookupOrDefault("transonic", false);
        
//    const bool consistent =
//        pimple.lookupOrDefault("consistent", false);
        

Switch updateEnergyInPISO (false);
word fluxScheme (word::null);

    pimple.lookup("updateEnergyInPISO") >> updateEnergyInPISO;
        
    pimple.lookup("fluxScheme") >> fluxScheme;

    if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
    {
    }
    else
    {
        FatalErrorIn
        (
        "readAdditionalPimpleControl.H "
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Tadmor, Kurganov"
            << abort(FatalError);
    }

