MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho,U)
    +
    fvm::div(phi_own,U)
    +
    fvm::div(phi_nei,U)
    +
    MRF.DDt(rho,U)
    +
//multiRegionReactingPimpleCentralFoam
    turb.divDevRhoReff(U)
    ==
    fvOptions(rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

fvOptions.constrain(UEqn);

//multiRegionReactingPimpleCentralFoam
//if (pimple.momentumPredictor())
if (momentumPredictor)
{
    solve
    (
        UEqn == -gradp
    );
    
    fvOptions.correct(U);
}
