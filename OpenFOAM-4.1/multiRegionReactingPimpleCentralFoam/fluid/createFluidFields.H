// Initialise fluid field pointer lists
PtrList<combustionModels::psiCombustionModel> reactionFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volScalarField> CpFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> kappaFluid(fluidRegions.size());

PtrList<surfaceScalarField> ownFluid(fluidRegions.size());
PtrList<surfaceScalarField> neiFluid(fluidRegions.size());
PtrList<surfaceScalarField> tuMagSfFluid(fluidRegions.size());
PtrList<surfaceVectorField> tuSfFluid(fluidRegions.size());

PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
PtrList<fv::kappaFunction> kappaFuncPtr(fluidRegions.size());    
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
PtrList<volScalarField> dQFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable> fieldsFluid(fluidRegions.size());

PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());


#include "createRDeltaT.H"

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{

    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;
        
        
    Info<< "Creating reaction model and thermophysical properties\n" << endl;

    reactionFluid.set
    (
        i,
          autoPtr<combustionModels::psiCombustionModel>(combustionModels::psiCombustionModel::New(fluidRegions[i])).ptr()
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            reactionFluid[i].thermo().rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );


    kappaFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
               "kappa",
                runTime.timeName(),
                fluidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimless
        )
    );


    CpFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Cp",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
           reactionFluid[i].thermo().Cp()
        )
    );
    
    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            reactionFluid[i].thermo()
        ).ptr()
    );
    
    reactionFluid[i].setTurbulence(turbulence[i]);

       
    radiation.set
    (
        i,
        radiation::radiationModel::New(reactionFluid[i].thermo().T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                "dpdt",
                reactionFluid[i].thermo().p().dimensions()/dimTime,
                0
            )
        )
    );

   Info<< "    Adding to dQ\n" << endl;
   dQFluid.set
   (
        i,
        new volScalarField
        (
            IOobject
            (
                "dQ",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
                fluidRegions[i],
                dimensionedScalar
            (
                "dQ",
                dimEnergy/dimTime,
                0
            )
        )
    );


    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    rhoMax.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMax",
                pimpleDict,
                dimDensity,
                GREAT
            )
        )
    );

    rhoMin.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMin",
                pimpleDict,
                dimDensity,
                0
            )
        )
    );


    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;

    fluidFvOptions.set
     (
        i,
        new fv::options(fluidRegions[i])
     );
 
 
    fieldsFluid.set
    	(
    	    i,
			new multivariateSurfaceInterpolationScheme<scalar>::fieldTable 
		);

    basicMultiComponentMixture& composition = reactionFluid[i].thermo().composition();	
		
	
	forAll(composition.Y(), j)
	{
		    fieldsFluid[i].add(composition.Y()[i]);
	}



//Central
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    ownFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "own",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("own", dimless, 1.0)
        )
    );

    neiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "nei",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("nei", dimless, -1.0)
        )
    );


    tuMagSfFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "uMagSf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i].magSf() * 1.0
        )
    );
    

    tuSfFluid.set
    (
        i,
        new surfaceVectorField
        (
            IOobject
            (
                "uSf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i].Sf() * 1.0
        )
    );



    correctCentralACMIInterpolation(tuMagSfFluid[i]);
    correctCentralACMIInterpolation(tuSfFluid[i]);


    

    
    fieldsFluid[i].add(reactionFluid[i].thermo().he());

    turbulence[i].validate();    

    //
     //
     // Flux scheme switching surface field
     //
     //
    kappaFuncPtr.set
    (
        i,
        fv::kappaFunction::New
        (
            "kappaFunction",
            fluidRegions[i].solutionDict().subDict("PIMPLE").subDict("kappaFunction"),
            fluidRegions[i]
        ).ptr()
    );
    


}




