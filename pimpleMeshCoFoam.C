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
    pimpleDyMFoam.C

Group
    grpIncompressibleSolvers grpMovingMeshSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControls.H" // local file, here we source maxMeshCo and minDeltaT
    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    Info<< "Max Co: " << maxCo << endl;
    Info<< "Max MeshCo: " << maxMeshCo << endl;


    while (runTime.run())
    {
        #include "readControls.H"
        #include "CourantNo.H"


        if (adjustTimeStep && mesh.changing())
        {

            scalar meshCoNum = 0.0;
            scalar meanMeshCoNum = 0.0;

            if (mesh.nInternalFaces() && mesh.moving())
            {
                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
                );

                meshCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

                meanMeshCoNum =
                    0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

                Info<< "Mesh Courant Number mean: " << meanMeshCoNum
                    << " max: " << meshCoNum << endl;
            }

            // evaluate delta T for Mesh Courant no
            scalar maxDeltaTFact = 
                min(maxMeshCo/(meshCoNum + SMALL), maxCo/(CoNum + SMALL));

            scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

            // Debug statement
            Info << "DeltaTFact = " << deltaTFact << endl;

            // set delta T to be smaller than maxDeltaT but bigger than minDeltaT
            // Currently, minDeltaT is not checked for
            runTime.setDeltaT
            (
                min
                (
                    deltaTFact*runTime.deltaTValue(),
                    maxDeltaT
                )
            );

            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }


        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

        // Calculate absolute flux from the mapped surface velocity
        phi = mesh.Sf() & Uf;

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
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
