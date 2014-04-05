/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "immersedBoundary.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

    immersedBoundary ib(mesh);

#   include "createIbFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

//     volVectorField gradP = ib.gamma()*ib.grad(p);
    volVectorField gradP = ib.gamma()*fvc::grad(p);

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        for (int outerCorr=0; outerCorr<nOuterCorr; outerCorr++)
        {
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              - fvm::laplacian(nu, U)
            );

            ib.imposeDirichletBoundaryCondition(UEqn, ibVelocity);

            solve(UEqn == -gradP);

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rUA = 1.0/UEqn.A();
                
                scalarField ibPhi = ib.flux(U);

                U = rUA*UEqn.H();
                phi = fvc::interpolate(U) & mesh.Sf();

                ib.setFlux(phi, ibPhi);

                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian
                        (
                            ib.sGamma()*fvc::interpolate(rUA),
                            p
                        )
                     == fvc::div(phi)
                    );

                    ib.correctEqn(pEqn);

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve();

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrsIb.H"

                ibPressure = 
                    ib.imposeNeumannBoundaryCondition(p, ibSnGradP);
                p.correctBoundaryConditions();

//                 gradP = ib.gamma()*ib.grad(p);
                gradP = ib.gamma()*fvc::grad(p);

                U -= rUA*gradP;
                ibSnGradU = 
                    ib.imposeDirichletBoundaryCondition(UEqn, ibVelocity);
                U.correctBoundaryConditions();
            }
        }

        ibPressure = 
            ib.imposeNeumannBoundaryCondition(p, ibSnGradP);
        p.correctBoundaryConditions();

        runTime.write();

        scalarField uTau = sqrt(nu.value()*mag(ibSnGradU));
        scalarField yPlus = uTau*ib.delta()/nu.value();

        Info << "yPlus, max: " << max(yPlus) 
            << ", avg: " << average(yPlus) << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
