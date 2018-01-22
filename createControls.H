#include "createControl.H"
#include "createTimeControls.H"

scalar maxMeshCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxMeshCo", 1.0);

scalar minDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);
