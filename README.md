# hybridCentralSolvers
United collection of hybrid  Central solvers - one-phase, two-phase and multicomponent versions.

Latest OpenFOAM version: 4.1

OpenFOAM-dev checked with commit f5b91be 6-Feb-17. Note OpenFOAM 4.1 tutorials work with OpenFOAM-dev as of commit f5b91be.

1. **pimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
Kurganov and Tadmor with LTS support for steady-state calculations

2. **pimpleCentralDyMFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
Kurganov and Tadmor with mesh motion and LTS support for steady-state calculations

3. **reactingPimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of 
Kurganov and Tadmor for combustion with chemical reactions and LTS support for steady-state calculations

4. **multiRegionReactingPimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind
schemes of Kurganov and Tadmor for combustion with chemical reactions and LTS support for steady-state calculations. Multi-
Region version includes solid heat conduction with conjugate heat transfer between solid and fluid regions.
