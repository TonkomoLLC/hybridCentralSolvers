        psi.oldTime();
        rho.oldTime();
        p.oldTime();
        U.oldTime();
        h.oldTime();
        K.oldTime();
        
if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
      #include "hEqnFrozenFlow.H"
}
else
{

        // --- Predict density
        #include "massEqn.H"

        // --- Solve momentum
        #include "UEqn.H"

       // --- Solve energy
        if (!updateEnergyInPISO)
        {
          #include "YEqn.H"
          #include "hEqn.H"
        }


    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {

            turb.correct();

          if (updateEnergyInPISO) //update each iteration before pressure
          {
              #include "YEqn.H"
              #include "hEqn.H"
          }


               #include "pEqn.H"

          if (updateEnergyInPISO)
                {
                    #define PISOCENTRALFOAM_LTS
                    
                    //// --- update weightings for central scheme
                    //#include "updateCentralWeights.H"
                    
                    // --- update blending function
                    #include "updateKappa.H"
                    
                    // --- update mechanical fields
                    #include "updateMechanicalFields.H"
                }

    }

    if (!updateEnergyInPISO)
    {
           //// --- update weightings for central scheme
           //#include "updateCentralWeights.H"
                
           // --- update blending function
           #include "updateKappa.H"
                
           // --- update mechanical fields
           #include "updateMechanicalFields.H"
    }


}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}

