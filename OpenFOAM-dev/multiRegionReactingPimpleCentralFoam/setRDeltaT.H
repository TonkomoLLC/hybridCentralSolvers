{

    fvMesh& mesh = fluidRegions[0];
    combustionModels::psiCombustionModel& reaction = reactionFluid[0];
    surfaceScalarField& phi = phiFluid[0];
    volScalarField& rho = rhoFluid[0];    
    PtrList<volScalarField>& Y = reactionFluid[0].thermo().composition().Y();
    basicMultiComponentMixture& composition = reactionFluid[0].thermo().composition();	
    psiReactionThermo& thermo = reaction.thermo();
    surfaceScalarField& own = ownFluid[0];
    surfaceScalarField& nei = neiFluid[0];
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    const volScalarField& psi = thermo.psi();
    volVectorField& U = UFluid[0];

        surfaceScalarField rho_own
    (
        "rho_own",
        fvc::interpolate(thermo.rho()(), own, "reconstruct(rho)")
    );

    surfaceScalarField rho_nei
    (
        "rho_nei",
        fvc::interpolate(thermo.rho()(), nei, "reconstruct(rho)")
    );

    surfaceScalarField psi_own
    (
        "psi_own",
        fvc::interpolate(psi, own, "reconstruct(psi)")
    );
    surfaceScalarField psi_nei
    (
        "psi_nei",
        fvc::interpolate(psi, nei, "reconstruct(psi)")
    );

    surfaceVectorField psiU_own
    (
        "rhoU_own",
        fvc::interpolate(psi*U, own, "reconstruct(U)")
    );
    surfaceVectorField psiU_nei
    (
        "rhoU_nei",
        fvc::interpolate(psi*U, nei, "reconstruct(U)")
    );


    surfaceScalarField phiv_own
    (
        "phiv_own",
        (psiU_own / psi_own)  & mesh.Sf()
    );
    surfaceScalarField phiv_nei
    (
        "phiv_nei",
        (psiU_nei / psi_nei) & mesh.Sf()
    );

volScalarField c
(
    "cSound",
    sqrt(thermo.Cp()/thermo.Cv() / psi)
);

    surfaceScalarField cf_own
    (
        "cf_own",
        fvc::interpolate(c, own, "reconstruct(psi)")
    );

    surfaceScalarField cf_nei
    (
        "cf_nei",
        fvc::interpolate(c, nei, "reconstruct(psi)")
    );
    correctCentralACMIInterpolation(cf_nei);

    surfaceScalarField cSf_own
    (
        "cSf_own",
        cf_own * mesh.magSf()
    );
    surfaceScalarField cSf_nei
    (
        "cSf_nei",
        cf_nei * mesh.magSf()
    );


    surfaceScalarField ap
    (
        "ap",
        max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero)
    );



        surfaceScalarField alpha_own
    (
        "alpha_own",
        ap/(ap - am)
    );
    surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
    surfaceScalarField aSf
    (
        "aSf",
        am*alpha_own*0.0
    );
    surfaceScalarField alpha_nei
    (
        "alpha_nei",
        1.0 - alpha_own
    );

//OpenFOAM-dev
/*
    const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!composition.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie << " not found in available species "
            << composition.species()
            << exit(FatalIOError);
    }
    const label inertIndex(composition.species()[inertSpecie]);   
*/

    volScalarField& rDeltaT = trDeltaT.ref();

    scalar maxCo
    (
        readScalar(runTime.controlDict().lookup("maxCo")) //0.8 is default
    );

    scalar rDeltaTSmoothingCoeff
    (
        readScalar(runTime.controlDict().lookup("rDeltaTSmoothingCoeff")) // 0.02 is default
    );

    scalar rDeltaTDampingCoeff
    (
        readScalar(runTime.controlDict().lookup("rDeltaTDampingCoeff")) // 1.0 is default - no damping
    );

    scalar maxDeltaT
    (
        readScalar(runTime.controlDict().lookup("maxDeltaT")) // GREAT - is a default
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);
    
    if (initialDeltaT < 0.0)
    {
        runTime.controlDict().lookup("initialDeltaT") >> initialDeltaT;
        rDeltaT.primitiveFieldRef() = 1.0 / initialDeltaT;
        rDeltaT0.primitiveFieldRef() = rDeltaT.primitiveField();
    }
    else
    {
        if ( CourantType == "cellCourant" )
        {
            // Set the reciprocal time-step from the local Courant number
    
            rDeltaT.ref() = max
            (
                1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
                (
                    0.5 * fvc::surfaceSum(mag(phi))().internalField() 
                    / mesh.V() / rho.internalField()
                ) / maxCo
            );
        }
        else if ( CourantType == "faceCourant" )
        {
            if (mesh.nInternalFaces())
            {
                surfaceScalarField phiOverRho =
                mag
                (
                    phi / (rho_own * alpha_own + rho_nei * alpha_nei)
                );
            
                surfaceScalarField SfUfbyDelta =
                    mesh.surfaceInterpolation::deltaCoeffs()*phiOverRho*2.0;
            
                forAll(mesh.V(), iCell)
                {
                    const labelList& cellFaces = mesh.cells()[iCell];
                    
                    scalar maxCellRDeltaT = -1.0;
                    scalar faceRDeltaT = maxCellRDeltaT;
                    label faceId = -1;
                    forAll(cellFaces, jFace)
                    {
                        faceId = cellFaces[jFace];
                        if (mesh.isInternalFace(faceId))
                        {
                            faceRDeltaT = SfUfbyDelta.internalField()[faceId] / mesh.magSf().internalField()[faceId] / maxCo;
                        }
                        if (faceRDeltaT > maxCellRDeltaT)
                        {
                            maxCellRDeltaT = faceRDeltaT;
                        }
                    }
                    rDeltaT.primitiveFieldRef()[iCell] = maxCellRDeltaT;
                }
                forAll(SfUfbyDelta.boundaryField(), kPatch)
                {
                    if(polyPatch::constraintType(mesh.boundaryMesh().types()[kPatch]))
                    {
                    }
                    else
                    {
                        const fvsPatchField<scalar>& pSfUfbyDelta = SfUfbyDelta.boundaryField()[kPatch];
                        const labelList& patchCells = mesh.boundary()[kPatch].faceCells();
                        scalar faceRDeltaT = -1;
                        forAll(pSfUfbyDelta, qFace)
                        {
                            faceRDeltaT = pSfUfbyDelta[qFace] / mesh.magSf().boundaryField()[kPatch][qFace] / maxCo;
                            if (faceRDeltaT > rDeltaT.primitiveField()[patchCells[qFace]])
                            {
                                rDeltaT.primitiveFieldRef()[patchCells[qFace]] = faceRDeltaT;
                            }
                        }
                    }
                }
           }
        
        }
        else
        {
            FatalErrorIn
            (
                "setRDeltaT.H:"
            )
            << "Wrong type of Courant criterion: " << CourantType
            << endl << " must be one of:" 
            << endl << "1) cellCourant"
            << endl << "2) faceCourant"
            << endl << abort(FatalError);
        }
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/(rDeltaT.primitiveField() + 1.0e-12))
        << ", " << gMax(1/(rDeltaT.primitiveField() + 1.0e-12)) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
