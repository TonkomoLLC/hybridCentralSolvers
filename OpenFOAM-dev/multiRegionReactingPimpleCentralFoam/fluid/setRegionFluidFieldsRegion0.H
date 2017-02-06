    fvMesh& mesh = fluidRegions[0];
 
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields = fieldsFluid[0];
 
    combustionModels::psiCombustionModel& reaction = reactionFluid[0];
    psiReactionThermo& thermo = reaction.thermo();
    PtrList<volScalarField>& Y = reactionFluid[0].thermo().composition().Y(); 
    thermo.validate(args.executable(), "h", "e");
    volScalarField& p = thermo.p();
    volScalarField& h = thermo.he();
    volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    volScalarField& rho = rhoFluid[0];
    volVectorField& U = UFluid[0];
    surfaceScalarField& phi = phiFluid[0];
    surfaceScalarField& kappa = kappaFluid[0];
    surfaceScalarField& own = ownFluid[0];
    surfaceScalarField& nei = neiFluid[0];
    surfaceScalarField& tuMagSf = tuMagSfFluid[0];
    surfaceVectorField& tuSf = tuSfFluid[0];

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

    surfaceScalarField p_own
    (
        "p_own",
        fvc::interpolate(p, own, "reconstruct(p)")
    );
    surfaceScalarField p_nei
    (
        "p_nei",
        fvc::interpolate(p, nei, "reconstruct(p)")
    );

    compressible::turbulenceModel& turb = turbulence[0];
    fv::kappaFunction& kappaFunc = kappaFuncPtr[0];    

    reaction.setTurbulence(turb);
    
    volScalarField& K = KFluid[0];
    volScalarField& dpdt = dpdtFluid[0];


    volScalarField& Qdot = QdotFluid[0];

    radiation::radiationModel& rad = radiation[0];

    fv::options& fvOptions = fluidFvOptions[0];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[0]
    );

    basicMultiComponentMixture& composition = reactionFluid[0].thermo().composition();	

    const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!composition.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie << " not found in available species "
            << composition.species()
            << exit(FatalIOError);
    }

dimensionedScalar v_one("v_one", pow(dimTime,-1.0), 1.0);
volScalarField rAU ("rAU", 1.0/rho/v_one);
rAU.primitiveFieldRef() *= runTime.deltaTValue();
rAU.primitiveFieldRef() /= mesh.V().field();
volVectorField HbyA ("HbyA", U);
 dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
volScalarField KChange
(
    "KChange",
    rho*K / runTime.deltaT() * 0.0
);

volScalarField KSource
(
    "KSource",
    p / runTime.deltaT() * 0.0
);


volVectorField gradp
(
    "gradp",
    fvc::grad(p)
);
volVectorField limGradP 
(
    "limGradP",
    fvc::grad(p)
);

DynamicList<label> badQualityCells(0);

Switch updateEnergyInPISO (false);
word fluxScheme (word::null);

    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(rho_nei);

    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);

    
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

    if (fluxScheme == "Kurganov")
    {
        //do nothing, field is updated at each time step
    }
    else if (fluxScheme == "Tadmor")
    {
        //fields are equal
        alpha_own = 0.5;
        alpha_nei = 1.0 - alpha_nei;
    }
    else
    {
        
    }
    
    phiv_own *= alpha_own;
    phiv_nei *= alpha_nei;
    
    surfaceScalarField aphiv_own("aphiv_own", phiv_own - aSf);
    surfaceScalarField aphiv_nei("aphiv_nei", phiv_nei + aSf);
    
    surfaceScalarField phid_own
    (
        "phid_own",
        aphiv_own * psi_own
    );
    surfaceScalarField phid_nei
    (
        "phid_nei",
        aphiv_nei * psi_nei
    );

    surfaceScalarField phi_own
    (
        "phi_own",
        phi * 1.0
    );
    surfaceScalarField phi_nei
    (
        "phi_nei",
        phi * 0.0
    );

    surfaceScalarField Dp_own
    (
        "Dp_own",
        alpha_own * dimensionedScalar ("zeroA", dimTime, 0.0)
    );
    surfaceScalarField Dp_nei
    (
        "Dp_nei",
        Dp_own
    );




    forAll(kappa, iFace)
    {
        kappa[iFace] = 0.0;
    }
    forAll(kappa.boundaryField(), iPatch)
    {
        kappa.boundaryFieldRef()[iPatch] = 0.0;
    }


    const label inertIndex(composition.species()[inertSpecie]);   
    const bool frozenFlow = frozenFlowFluid[0];
    
    
    PtrList<surfaceScalarField> mDCf(Y.size());
    PtrList<surfaceScalarField> diffusiveFlux(Y.size());
    PtrList<surfaceScalarField> phiY(Y.size());
    for (label cmptIdx=0; cmptIdx < Y.size(); cmptIdx++)
    {
        mDCf.set
        (
            cmptIdx,
            new surfaceScalarField
            (
                "mDCf",
                fvc::interpolate(turb.muEff())
            )
        );
        
        diffusiveFlux.set
        (
            cmptIdx,
            new surfaceScalarField
            (
                "phiD",
                mDCf[cmptIdx] * fvc::snGrad(Y[cmptIdx]) * mesh.magSf() * 0.0
            )
        );
        
        if (Y[cmptIdx].name() == inertSpecie)
        {
            phiY.set
            (
                cmptIdx,
                new surfaceScalarField
                (
                    "rhoiPhi",
                    phi * 0.0
                )
            );
        }
    }

