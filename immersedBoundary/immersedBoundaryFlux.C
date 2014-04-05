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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *  //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<scalarField> immersedBoundary::flux
(
    const volVectorField& U
) const
{
    tmp<scalarField> tIbPhi(new scalarField(faces().size(), 0));
    scalarField& ibPhi = tIbPhi();

    const vectorField& UI = U.internalField();

    const vectorField& Sf = mesh_.faceAreas();
    const unallocLabelList& owner = mesh_.owner(); 
    const unallocLabelList& neighbour = mesh_.neighbour(); 
    const surfaceScalarField& w = mesh_.weights();
    
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        vector ibU = vector::zero;

        if (curFace < mesh_.nInternalFaces())
        {
            ibU = w[curFace]*UI[owner[curFace]] 
              + (1.0 - w[curFace])*UI[neighbour[curFace]];
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(curFace);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = faces()[faceI] - start;
                
            const unallocLabelList& fCells = 
                mesh_.boundary()[patchID].faceCells();

            const fvsPatchScalarField& pw = w.boundaryField()[patchID];

            if (U.boundaryField()[patchID].coupled())
            {
                ibU = pw[localIndex]*UI[fCells[localIndex]]
                  + (1.0 - pw[localIndex])
                   *U.boundaryField()[patchID][localIndex];
            }
            else
            {
                ibU = U.boundaryField()[patchID][localIndex];
            }
        }

        ibPhi[faceI] = (ibU&Sf[curFace]);
    }

    // Set flux sign
    forAll(ibPhi, faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            if (gamma()[owner[curFace]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(curFace);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = faces()[faceI] - start;
                
            const unallocLabelList& fCells = 
                mesh_.boundary()[patchID].faceCells();

            if (gamma().internalField()[fCells[localIndex]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
    }

    // Scale flux
    scalarField weights = mag(ibPhi);

    scalar sumWeights = gSum(weights);

    forAll(weights, faceI)
    {
        label curFace = faces()[faceI];

        if (curFace >= mesh_.nInternalFaces())
        {
            sumWeights -= weights[faceI];
        }
    }

    if(sumWeights > SMALL)
    {
        weights /= sumWeights;
    }

    scalar sumIbPhi = gSum(ibPhi);

    forAll(ibPhi, faceI)
    {
        label curFace = faces()[faceI];

        if (curFace >= mesh_.nInternalFaces())
        {
            sumIbPhi -= ibPhi[faceI];
        }
    }

    ibPhi -= weights*sumIbPhi;

    // Reset flux sign
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            if (gamma()[owner[curFace]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(curFace);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = faces()[faceI] - start;

            const unallocLabelList& fCells = 
                mesh_.boundary()[patchID].faceCells();

            if (gamma()[fCells[localIndex]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
    }

    return tIbPhi;
}


tmp<scalarField> immersedBoundary::flux
(
    const volVectorField& U,
    const vectorField& ibVelocity
) const
{
    tmp<scalarField> tIbPhi(new scalarField(faces().size(), 0));
    scalarField& ibPhi = tIbPhi();

    const vectorField& UI = U.internalField();

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const unallocLabelList& owner = mesh_.owner(); 

    label nCoeffs = 5;
    if(mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

#   include "sendAndReceiveVelocity.H"

    forAll(cells(), cellI)
    {
        label curCell = cells()[cellI];

        const labelList& curCells = cellCells()[cellI];
        const List<labelPair>& curProcCells = cellProcCells()[cellI];

        const scalarRectangularMatrix& curInvMatrix =
            invDirichletMatrices()[cellI];

        Field<vector> coeffs(nCoeffs, pTraits<vector>::zero);
        Field<vector> source
        (
            curCells.size() + curProcCells.size(), 
            pTraits<vector>::zero
        );

        label pointID = 0;
        for (label i=0; i<curCells.size(); i++)
        {
            source[pointID++] = UI[curCells[i]] - ibVelocity[cellI];
        }

        for (label i=0; i<curProcCells.size(); i++)
        {
            source[pointID++] = 
                procVelocities
                [
                    curProcCells[i].first()
                ]
                [
                    curProcCells[i].second()
                ]
              - ibVelocity[cellI];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        const labelList& curCellFaces = mesh_.cells()[curCell];
        forAll(curCellFaces, faceI)
        {
            const label& curFace = curCellFaces[faceI];

            label index = findIndex(faces(), curFace);

            if (index != -1)
            {
                vector curCf = Cf[curFace] - points()[cellI];

                vector curUf =
                    ibVelocity[cellI]
                  + coeffs[0]*curCf.x()
                  + coeffs[1]*curCf.y()
                  + coeffs[2]*curCf.x()*curCf.y()
                  + coeffs[3]*sqr(curCf.x())
                  + coeffs[4]*sqr(curCf.y());

                if(mesh_.nGeometricD() == 3)
                {
                    curUf += 
                        coeffs[5]*curCf.z()
                      + coeffs[6]*curCf.x()*curCf.z()
                      + coeffs[7]*curCf.y()*curCf.z()
                      + coeffs[8]*sqr(curCf.z());
                }

                ibPhi[index] = (curUf&Sf[curFace]);
            }
        }
    }

    forAll(ibPhi, faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            if (gamma()[owner[curFace]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(curFace);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = faces()[faceI] - start;
                
            const unallocLabelList& fCells = 
                mesh_.boundary()[patchID].faceCells();

            if (gamma().internalField()[fCells[localIndex]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
    }

    scalarField weights = mag(ibPhi);

    scalar sumWeights = gSum(weights);

    forAll(weights, faceI)
    {
        label curFace = faces()[faceI];

        if (curFace >= mesh_.nInternalFaces())
        {
            sumWeights -= weights[faceI];
        }
    }

    if(sumWeights > SMALL)
    {
        weights /= sumWeights;
    }

    scalar sumIbPhi = gSum(ibPhi);

    forAll(ibPhi, faceI)
    {
        label curFace = faces()[faceI];

        if (curFace >= mesh_.nInternalFaces())
        {
            sumIbPhi -= ibPhi[faceI];
        }
    }

    ibPhi -= weights*sumIbPhi;

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            if (gamma()[owner[curFace]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(curFace);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = faces()[faceI] - start;

            const unallocLabelList& fCells = 
                mesh_.boundary()[patchID].faceCells();

            if (gamma()[fCells[localIndex]] < SMALL)
            {
                ibPhi[faceI] *= -1;
            }
        }
    }

    return tIbPhi;
}


void immersedBoundary::setFlux
(
    surfaceScalarField& phi,
    const scalarField& ibPhi
) const
{
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            phi.internalField()[curFace] = ibPhi[faceI];
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(curFace);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = curFace - start;

            phi.boundaryField()[patchID][localIndex] = ibPhi[faceI];
        }
    }
}

void immersedBoundary::scaleFlux
(
    surfaceScalarField& phi, 
    const surfaceScalarField& phiU
) const
{
//     // Set IB internal faces flux
//     forAll(internalFaces(), faceI)
//     {
//         label curFace = internalFaces()[faceI];

//         if (curFace < mesh_.nInternalFaces())
//         {
//             phi.internalField()[curFace] =
//                 phiU.internalField()[curFace];
                
//         }
//         else
//         {
//             label patchID = mesh_.boundaryMesh().whichPatch(curFace);
//             label start = mesh_.boundaryMesh()[patchID].start();
//             label localIndex = curFace - start;

//             phi.boundaryField()[patchID][localIndex] =
//                 phiU.boundaryField()[patchID][localIndex];
//         }
//     }


    // Scale flux for IB faces
    {
        const scalarField& gammaI = gamma().internalField();
        const unallocLabelList& owner = mesh_.owner(); 

        scalarField ibPhi(faces().size(), 0);

        forAll(faces(), faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh_.nInternalFaces())
            {
                ibPhi[faceI] = phiU.internalField()[curFace];
            }
            else
            {
                label patchID = mesh_.boundaryMesh().whichPatch(curFace);
                label start = mesh_.boundaryMesh()[patchID].start();
                label localIndex = curFace - start;
                
                ibPhi[faceI] = phiU.boundaryField()[patchID][localIndex];
            }
        }

        forAll(ibPhi, faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh_.nInternalFaces())
            {
                if (gammaI[owner[curFace]] < SMALL)
                {
                    ibPhi[faceI] *= -1;
                }
            }
            else
            {
                label patchID = mesh_.boundaryMesh().whichPatch(curFace);
                label start = mesh_.boundaryMesh()[patchID].start();
                label localIndex = faces()[faceI] - start;
                
                const unallocLabelList& fCells = 
                    mesh_.boundary()[patchID].faceCells();

                if (gamma().internalField()[fCells[localIndex]] < SMALL)
                {
                    ibPhi[faceI] *= -1;
                }
            }
        }

        scalarField weights = mag(ibPhi);

        scalar sumWeights = gSum(weights);

        forAll(weights, faceI)
        {
            label curFace = faces()[faceI];

            if (curFace >= mesh_.nInternalFaces())
            {
                sumWeights -= weights[faceI];
            }
        }

        if(sumWeights > SMALL)
        {
            weights /= sumWeights;
        }

        scalar sumIbPhi = gSum(ibPhi);

        forAll(ibPhi, faceI)
        {
            label curFace = faces()[faceI];

            if (curFace >= mesh_.nInternalFaces())
            {
                sumIbPhi -= ibPhi[faceI];
            }
        }

        ibPhi -= weights*sumIbPhi;


        // Set flux
        forAll(faces(), faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh_.nInternalFaces())
            {
                if (gammaI[owner[curFace]] > SMALL)
                {
                    phi.internalField()[curFace] = ibPhi[faceI];
                }
                else
                {
                    phi.internalField()[curFace] = -ibPhi[faceI];
                }
            }
            else
            {
                label patchID = mesh_.boundaryMesh().whichPatch(curFace);
                label start = mesh_.boundaryMesh()[patchID].start();
                label localIndex = faces()[faceI] - start;

                const unallocLabelList& fCells = 
                    mesh_.boundary()[patchID].faceCells();

                if (gamma().internalField()[fCells[localIndex]] < SMALL)
                {
                    phi.boundaryField()[patchID][localIndex] = -ibPhi[faceI];
                }
                else
                {
                    phi.boundaryField()[patchID][localIndex] = ibPhi[faceI];
                }
            }
        }
    }


//     // Scale flux for inside IB faces
//     {
//         const scalarField& gammaExtI = gammaExt().internalField();
//         const unallocLabelList& owner = mesh_.owner(); 

//         scalarField ibPhi(insideFaces().size(), 0);
//         forAll(insideFaces(), faceI)
//         {
//             if (insideFaces()[faceI] < mesh_.nInternalFaces())
//             {
//                 ibPhi[faceI] = 
//                     phiU.internalField()[insideFaces()[faceI]];
//             }
//             else
//             {
//                 label patchID = 
//                     mesh_.boundaryMesh().whichPatch(insideFaces()[faceI]);
//                 label start = mesh_.boundaryMesh()[patchID].start();
//                 label localIndex = insideFaces()[faceI] - start;

//                 ibPhi[faceI] = phiU.boundaryField()[patchID][localIndex];
//             }
//         }


//         forAll(ibPhi, faceI)
//         {
//             label curFace = insideFaces()[faceI];

//             if (curFace < mesh_.nInternalFaces())
//             {
//                 if (gammaExtI[owner[curFace]] < SMALL)
//                 {
//                     ibPhi[faceI] *= -1;
//                 }
//             }
//             else
//             {
//                 label patchID = mesh_.boundaryMesh().whichPatch(curFace);
//                 label start = mesh_.boundaryMesh()[patchID].start();
//                 label localIndex = insideFaces()[faceI] - start;
                
//                 const unallocLabelList& fCells = 
//                     mesh_.boundary()[patchID].faceCells();

//                 if (gammaExt().internalField()[fCells[localIndex]] < SMALL)
//                 {
//                     ibPhi[faceI] *= -1;
//                 }
//             }
//         }

//         scalarField weights = mag(ibPhi);

//         scalar sumWeights = gSum(weights);

//         forAll(weights, faceI)
//         {
//             label curFace = insideFaces()[faceI];

//             if (curFace >= mesh_.nInternalFaces())
//             {
//                 sumWeights -= weights[faceI];
//             }
//         }

//         if(sumWeights > SMALL)
//         {
//             weights /= sumWeights;
//         }

//         scalar sumIbPhi = gSum(ibPhi);

//         forAll(ibPhi, faceI)
//         {
//             label curFace = insideFaces()[faceI];

//             if (curFace >= mesh_.nInternalFaces())
//             {
//                 sumIbPhi -= ibPhi[faceI];
//             }
//         }

//         ibPhi -= weights*sumIbPhi;


//         // Set flux
//         forAll(insideFaces(), faceI)
//         {
//             label curFace = insideFaces()[faceI];

//             if (curFace < mesh_.nInternalFaces())
//             {
//                 if (gammaExtI[owner[curFace]] > SMALL)
//                 {
//                     phi.internalField()[curFace] = ibPhi[faceI];
//                 }
//                 else
//                 {
//                     phi.internalField()[curFace] = -ibPhi[faceI];
//                 }
//             }
//             else
//             {
//                 label patchID = mesh_.boundaryMesh().whichPatch(curFace);
//                 label start = mesh_.boundaryMesh()[patchID].start();
//                 label localIndex = insideFaces()[faceI] - start;

//                 const unallocLabelList& fCells = 
//                     mesh_.boundary()[patchID].faceCells();

//                 if (gammaExt().internalField()[fCells[localIndex]] < SMALL)
//                 {
//                     phi.boundaryField()[patchID][localIndex] = -ibPhi[faceI];
//                 }
//                 else
//                 {
//                     phi.boundaryField()[patchID][localIndex] = ibPhi[faceI];
//                 }
//             }
//         }
//     }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
