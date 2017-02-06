 //   fvSolution solutionDict(runTime);

    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    const int nCorr =
        pimple.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const int nOuterCorr =
        pimple.lookupOrDefault<int>("nOuterCorrectors", 1);

    const bool momentumPredictor =
        pimple.lookupOrDefault("momentumPredictor", true);
        
//    const bool transonic =
//        pimple.lookupOrDefault("transonic", false);
        
//    const bool consistent =
//        pimple.lookupOrDefault("consistent", false);
        

    Info << "Marking cells with bad quality to limit pressure gradient" << endl;
    
    badQualityCells.resize(0);
    
    const scalar nonOrthogonalityThreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("nonOrthogonalityThreshold")
        )
    );
    
    const scalar skewnessThreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("skewnessThreshold")
        )
    );

    cellQuality cq (mesh);
    
    scalarField cellNonOrtho (cq.nonOrthogonality());
    
    scalarField cellSkewness (cq.skewness());
    
    forAll(cellNonOrtho, iCell)
    {
        if ( (cellNonOrtho[iCell] >= nonOrthogonalityThreshold) ||
            (cellSkewness[iCell] >= skewnessThreshold) )
        {
            badQualityCells.append(iCell);
            
        }
    }


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

