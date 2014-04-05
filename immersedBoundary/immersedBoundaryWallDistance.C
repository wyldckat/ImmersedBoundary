/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "immersedBoundary.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void immersedBoundary::makeWallDistance() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeWallDistance() const : "
            << "create wall (and IB) distance field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (wallDistancePtr_)
    {
        FatalErrorIn("immersedBoundary::makeWallDistance() const")
            << "wall distance field already exists"
            << abort(FatalError);
    }

    wallDistancePtr_ = 
        new volScalarField
        (
            IOobject
            (
                "ibWallDist",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("dist", dimLength, GREAT)
        );

    wordList patchFieldTypes
    (
        mesh_.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    bool needWallDist = false;

    forAll(mesh_.boundary(), patchI)
    {
        if
        (
            mesh_.boundary()[patchI].type()
         == wallFvPatch::typeName
        )
        {
            patchFieldTypes[patchI] = 
                fixedValueFvPatchScalarField::typeName;

            needWallDist = true;
        }
    }

    if (cells().size()>0)
    {
        needWallDist = true;
    }

    if (needWallDist)
    {
        volScalarField psi
        (
            IOobject
            (
                "psi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dist", dimLength*dimLength, 0),
            patchFieldTypes
        );

        label iCorr = 0;
        scalar initialResidual = 0;
        lduMatrix::solverPerformance solverPerf;

        lduMatrix::debug = 0;

        do
        {
            fvScalarMatrix psiEqn
            (
                fvm::laplacian(psi)
             == -dimensionedScalar("1", dimless, 1)
            );

            scalarField ibPsi(points().size(), 0);

            imposeDirichletBoundaryCondition(psiEqn, ibPsi);

            solverPerf = psiEqn.solve();

            if(iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }
        }
        while
        (
            solverPerf.initialResidual() > 1e-8
         && ++iCorr < 100
        );

        lduMatrix::debug = 1;

        Info << solverPerf.solverName() << ": Solving for " << psi.name()
             << ", Initial residula = " << initialResidual
             << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations = " << iCorr << endl; 

        volVectorField gradPsi = fvc::grad(psi);

        (*wallDistancePtr_) = 
            gamma()
           *(
               sqrt((gradPsi&gradPsi) + 2*psi)
             - mag(gradPsi)
            );

        scalarField ibCellsDistance =
            mag
            (
                vectorField(mesh_.cellCentres(), cells())
              - points()
            );

        forAll(cells(), cellI)
        {
            (*wallDistancePtr_)[cells()[cellI]] = ibCellsDistance[cellI];
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const volScalarField& immersedBoundary::wallDistance() const
{
    if (!wallDistancePtr_)
    {
        makeWallDistance();
    }

    return *wallDistancePtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
