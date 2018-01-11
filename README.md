# PimpleMeshCoFoam

A modified version of OpenFOAM-v1706's pimpleDyMFoam. 
Takes up to four arguments to calculate time step size:

- maxCo
- maxMeshCo
- maxDeltaT
- minDeltaT

Time step will be chosen so that 

- Courant number < maxCo
- Mesh Courant number < maxMeshCo
- minDeltaT < deltaT < maxDeltaT

## Install
- Clone this repository
- cd pimpleDyMFoam
- source {path-to-OpenFOAM-v1706}/etc/bash
- wmake

Solver will be installed to $FOAM_USER_APPBIN

## Use
In a moving mesh case, add following keywords:

adjustTimeStep          yes;

maxCo                   value;
maxMeshCo               value;
maxDeltaT               value;
minDeltaT               value;


maxMeshCo should be < 1 to prevent negative cell volumes.
