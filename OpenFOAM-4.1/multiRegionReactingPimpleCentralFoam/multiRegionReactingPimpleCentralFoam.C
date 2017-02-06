/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiRegionReactingCentralPimpleFoam

Description
    Pressure-based semi implicit compressible flow solver based on central
    -upwind schemes of Kurganov and Tadmor for combustion with chemical 
    reactions and LTS support for steady-state calculations
    
    MultiRegion version includes solid heat conduction with conjugate heat 
    transfer between solid and fluid regions. It handles secondary fluid or 
    solid circuits which can be coupled thermally with the main fluid region. 
    i.e radiators, etc.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "psiCombustionModel.H"
#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "pimpleControl.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "coordinateSystem.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "multivariateScheme.H"
#include "fvcSmooth.H"
#include "localEulerDdtScheme.H"

//reactingCentralPimpleFoam includes
#include "gaussConvectionScheme.H"
#include "coupledFvsPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "kappaFunction.H"
#include "correctCentralACMIInterpolation.H"
#include "centralMULES.H"
#include "cellQuality.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"

    pimpleControl pimple(fluidRegions[0]);

    #include "createTimeControls.H"
    scalar initialDeltaT = -VGREAT;

    #include "createFields.H"
// reactingCentralPimpleFoam - added to createFluidFields and setFluidRegionFields
//    #include "createRDeltaT.H"
//    #include "readAdditionalPimpleControl.H"
//    #include "createCommonCentralFields.H"
//    #include "createMulticomponentSurfaceFields.H"

    #include "createControl.H"
    #include "readCourantType.H"
    #include "createMRF.H"

// use multiRegion version
//    #include "initContinuityErrs.H"
    #include "initContinuityErrs.H"
    #include "readSolidTimeControls.H"
    
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
// REVISIT FOR MULTIREGION
//    #include "setInitialMultiRegionDeltaT.H"

        #include "setRegionFluidFieldsRegion0.H"
        #include "markBadQualityCells.H"
        #include "psiUpdateCentralFields.H"
        #include "updateKappa.H"
        #include "updateCentralWeights.H"

    if (!LTS)
    {
        #include "setRegionFluidFieldsRegion0.H"
        #include "psiUpdateCentralFields.H"
        #include "updateKappa.H"
        #include "updateCentralWeights.H"
        #include "createCentralCourantNo.H"
        #include "centralCompressibleCourantNo.H"
        #include "setInitialDeltaT.H"
// REVISIT FOR MULTIREGION
//        #include "compressibleMultiRegionCourantNo.H"
//        #include "setMultiRegionDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
        #include "createTimeControls.H"
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
// REVISIT FOR MULTIREGION
//            #include "setRDeltaT.H"
        }
        else
 
          {
              #include "setRegionFluidFieldsRegion0.H"
              #include "psiUpdateCentralFields.H"
              #include "updateKappa.H"
              #include "updateCentralWeights.H"
              #include "createCentralCourantNo.H"
              #include "acousticCourantNo.H"
              #include "centralCompressibleCourantNo.H"
              #include "readTimeControls.H"
              #include "setDeltaT.H"
// REVISIT FOR MULTIREGION
//              #include "readFluidMultiRegionPIMPLEControls.H"
//              #include "compressibleMultiRegionCourantNo.H"
//              #include "solidRegionDiffusionNo.H"
//              #include "setMultiRegionDeltaT.H"
          }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
                // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            bool finalIter = oCorr == nOuterCorr-1;

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
// REVISIT FOR MULTIREGION
//                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
// REVISIT FOR MULTIREGION
//                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveSolid.H"
            }

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
