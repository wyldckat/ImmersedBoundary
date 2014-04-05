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
#include "zeroGradientFvPatchFields.H"
#include "volMesh.H"
#include "SortableList.H"
#include "skewCorrectionVectors.H"
#include "leastSquaresGrad.H"

#include "OPstream.H"
#include "IPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(immersedBoundary, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void immersedBoundary::makeTriSurfSearch() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeTriSurfSearch() const : "
            << "creating triSurface search algorithm"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already
    if (triSurfSearchPtr_)
    {
        FatalErrorIn("immersedBoundary::makeTriSurfSearch() const")
            << "triSurface search algorithm already exist"
            << abort(FatalError);
    }

    triSurfSearchPtr_ = new triSurfaceSearch(*this);
}


void immersedBoundary::makeGamma() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeGamma() const : "
            << "creating fluid cells indicator"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaPtr_)
    {
        FatalErrorIn("immersedBoundary::makeGamma() const")
            << "fluid cells indicator already exist"
            << abort(FatalError);
    }

    gammaPtr_ =
        new volScalarField
        (
            IOobject
            (
                "ibGamma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1),
            zeroGradientFvPatchScalarField::typeName
        );

    scalarField& gammaI = gammaPtr_->internalField();

    gammaI = gammaExt().internalField();

    forAll(cells(), cellI)
    {
        gammaI[cells()[cellI]] = 0.0;
    }

    gammaPtr_->correctBoundaryConditions();
}


void immersedBoundary::makeGammaExt() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeGammaExtended() const : "
            << "creating extended fluid cells indicator"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaExtPtr_)
    {
        FatalErrorIn("immersedBoundary::makeGammaExtended() const")
            << "extended fluid cells indicator already exist"
            << abort(FatalError);
    }

    gammaExtPtr_ =
        new volScalarField
        (
            IOobject
            (
                "ibGammaExt",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1),
            zeroGradientFvPatchScalarField::typeName
        );

    scalarField& gammaExtI = gammaExtPtr_->internalField();

    const vectorField& C = mesh_.cellCentres();

    boolList inside = triSurfSearch().calcInside(C);

    if (!internalFlow_)
    {
        forAll(gammaExtI, cellI)
        {
            if(inside[cellI])
            {
                gammaExtI[cellI] = 0;
            }
        }
    }
    else
    {
        forAll(gammaExtI, cellI)
        {
            if(!inside[cellI])
            {
                gammaExtI[cellI] = 0;
            }
        }
    }

    gammaExtPtr_->correctBoundaryConditions();
}


void immersedBoundary::makeSGamma() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeSGamma() const : "
            << "creating fluid faces indicator"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (sGammaPtr_)
    {
        FatalErrorIn("immersedBoundary::makeSGamma() const")
            << "fluid faces indicator already exist"
            << abort(FatalError);
    }

    sGammaPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "sGamma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1),
            calculatedFvsPatchScalarField::typeName
        );

    scalarField& sGammaI = sGammaPtr_->internalField();

    forAll(faces(), faceI)
    {
        if (faces()[faceI] < mesh_.nInternalFaces())
        {
            sGammaI[faces()[faceI]] = 0.0;
        }
        else
        {
            label patchID = mesh_.boundaryMesh().whichPatch(faces()[faceI]);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = faces()[faceI] - start;

            sGammaPtr_->boundaryField()[patchID][localIndex] = 0.0;
        }
    }

    forAll(insideFaces(), faceI)
    {
        if (insideFaces()[faceI] < mesh_.nInternalFaces())
        {
            sGammaI[insideFaces()[faceI]] = 0.0;
        }
        else
        {
            label patchID = 
                mesh_.boundaryMesh().whichPatch(insideFaces()[faceI]);
            label start = mesh_.boundaryMesh()[patchID].start();
            label localIndex = insideFaces()[faceI] - start;

            sGammaPtr_->boundaryField()[patchID][localIndex] = 0.0;
        }
    }

//     forAll(internalFaces(), faceI)
//     {
//         if (internalFaces()[faceI] < mesh_.nInternalFaces())
//         {
//             sGammaI[internalFaces()[faceI]] = 0.0;
//         }
//         else
//         {
//             label patchID = 
//                 mesh_.boundaryMesh().whichPatch(internalFaces()[faceI]);
//             label start = mesh_.boundaryMesh()[patchID].start();
//             label localIndex = internalFaces()[faceI] - start;

//             sGammaPtr_->boundaryField()[patchID][localIndex] = 0.0;
//         }
//     }
}


void immersedBoundary::makeCells() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeIbCells() const : "
            << "create list of cells next to immersed boundary"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeIbCells() const")
            << "list of cells next to immersed boundary already exist"
            << abort(FatalError);
    }

    labelHashSet ibCellSet;

    const unallocLabelList& owner = mesh_.owner(); 
    const unallocLabelList& neighbour = mesh_.neighbour();

    const scalarField& gammaExtI = gammaExt().internalField();

    forAll(neighbour, faceI) 
    {
        if(mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            if(gammaExtI[owner[faceI]] > SMALL)
            {
                if(!ibCellSet.found(owner[faceI]))
                {
                    ibCellSet.insert(owner[faceI]);
                }
            }
            else
            {
                if(!ibCellSet.found(neighbour[faceI]))
                {
                    ibCellSet.insert(neighbour[faceI]);
                }
            }
        }
    }

    forAll(gammaExt().boundaryField(), patchI)
    {
        if (gammaExt().boundaryField()[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gammaExt().boundaryField()[patchI].patchInternalField();

            scalarField gammaExtNgb =
                gammaExt().boundaryField()[patchI].patchNeighbourField();

            const unallocLabelList& fCells = 
                mesh_.boundary()[patchI].faceCells();
            

            forAll(gammaExtOwn, faceI)
            {
                if
                (
                    mag(gammaExtNgb[faceI] - gammaExtOwn[faceI])
                  > SMALL
                )
                {
                    if(gammaExtOwn[faceI] > SMALL)
                    {
                        if(!ibCellSet.found(fCells[faceI]))
                        {
                            ibCellSet.insert(fCells[faceI]);
                        }
                    }
                    else if (2*gammaExtOwn.size() == fCells.size())
                    {
                        if
                        (
                           !ibCellSet.found
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            )
                        )
                        {
                            ibCellSet.insert
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            );
                        }
                    }
                }
            }
        }
    }

    cellsPtr_ = new labelList(ibCellSet.toc());

    Pout << "Num of IB cells: " << cellsPtr_->size() << endl;
}


void immersedBoundary::addCornerCells() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::addCornerIbCells() const : "
            << "add cells next to sharp corner"
            << endl;
    }

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const unallocLabelList& own = mesh_.owner(); 
    const unallocLabelList& ngb = mesh_.neighbour();

    label nCornerCells = 0;
    labelList cornerCells;

    do
    {
        const scalarField& gammaI = gamma().internalField();

        labelHashSet cornerIbCellSet;

        forAll(faces(), faceI)
        {
            const label& ownCell = own[faces()[faceI]];
            const label& ngbCell = ngb[faces()[faceI]];

            label ibCell = -1;
            label liveCell = -1;

            if (gammaI[ownCell] > SMALL)
            {
                liveCell = ownCell;
                ibCell = ngbCell;
            }
            else
            {
                liveCell = ngbCell;
                ibCell = ownCell;
            }

            scalar delta = cellSize(liveCell);
//             delta *= 2;
            vector span(delta, delta, delta);

            pointIndexHit pih = 
                triSurfSearch().nearest(C[liveCell], span);

            if(pih.hit())
            {
                vector p = pih.hitPoint();
                vector n =
                    triSurfaceTools::surfaceNormal
                    (
                        *this,
                        pih.index(),
                        pih.hitPoint()
                    );
            
                scalar totalArea = 0;
                {
                    const labelList& cellFaces = mesh_.cells()[liveCell];

                    forAll(cellFaces, faceI)
                    {
                        label curFace = cellFaces[faceI];
                    
                        vector curSf = Sf[curFace];
                    
                        if ((curSf&(Cf[curFace] - C[liveCell])) < 0)
                        {
                            curSf *= -1;
                        }

                        if ((curSf&n) > 0)
                        {
                            totalArea += (curSf&n);
                        }
                    }
                }

                scalar area = 0;
                {
                    const labelList& cellFaces = mesh_.cells()[liveCell];

                    forAll(cellFaces, faceI)
                    {
                        label curFace = cellFaces[faceI];
                    
                        vector curSf = Sf[curFace];

                        label ngbCell = -1;
                        {
                            if (own[curFace] == liveCell)
                            {
                                ngbCell = ngb[curFace];
                            }
                            else
                            {
                                ngbCell = own[curFace];
                            }
                        }

                        label ibCell = findIndex(cells(), ngbCell);

                        if (ibCell != -1)
                        {
                            if ((normals()[ibCell]&n) > 0)
                            {
                                area += mag(curSf&n);
                            }
                        }
                    }
                }

                if (area/totalArea < 0.5)
                {
                    if (!cornerIbCellSet.found(liveCell))
                    {
                        cornerIbCellSet.insert(liveCell);
                    }
                }                
            }
        }

        cornerCells = cornerIbCellSet.toc();
        cellsPtr_->append(cornerCells);
        nCornerCells += cornerCells.size();

        deleteDemandDrivenData(gammaPtr_);
        deleteDemandDrivenData(facesPtr_);
        deleteDemandDrivenData(pointsPtr_);
        deleteDemandDrivenData(normalsPtr_);        
        deleteDemandDrivenData(hitFacesPtr_);
    }
    while(cornerCells.size() > 0);

    Info << "Num of corner cells: " << nCornerCells << endl;
}


void immersedBoundary::makeFaces() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeFaces() const : "
            << "create list of faces next to immersed boundary"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (facesPtr_)
    {
        FatalErrorIn("immersedBoundary::makeFaces() const")
            << "list of faces next to immersed boundary already exist"
            << abort(FatalError);
    }

    labelHashSet ibFacesSet;

    const unallocLabelList& owner = mesh_.owner(); 
    const unallocLabelList& neighbour = mesh_.neighbour();

    const scalarField& gammaI = gamma().internalField();

    forAll(neighbour, faceI) 
    { 
        if (mag(gammaI[neighbour[faceI]] - gammaI[owner[faceI]]) > SMALL)
        {
            ibFacesSet.insert(faceI);
        }
    }

    forAll(gamma().boundaryField(), patchI)
    {
        if (gamma().boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gamma().boundaryField()[patchI].patchInternalField();

            scalarField gammaNgb =
                gamma().boundaryField()[patchI].patchNeighbourField();

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll(gammaOwn, faceI)
            {
                if
                (
                    mag(gammaNgb[faceI] - gammaOwn[faceI]) > SMALL
                )
                {
                    if (!ibFacesSet.found(start + faceI))
                    {
                        ibFacesSet.insert(start + faceI);
                    }
                    
                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibFacesSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibFacesSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    facesPtr_ = new labelList(ibFacesSet.toc());

//     Info << "nIbFaces: " << facesPtr_->size() << endl;
}


void immersedBoundary::makeInsideFaces() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeIbInsideFaces() const : "
            << "create list of faces next to immersed boundary"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (insideFacesPtr_)
    {
        FatalErrorIn("immersedBoundary::makeInsideFaces() const")
            << "list of faces next to immersed boundary already exist"
            << abort(FatalError);
    }

    labelHashSet ibInsideFacesSet;

    const unallocLabelList& owner = mesh_.owner(); 
    const unallocLabelList& neighbour = mesh_.neighbour();

    const scalarField& gammaExtI = gammaExt().internalField();

    forAll(neighbour, faceI) 
    { 
        if(mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            ibInsideFacesSet.insert(faceI);
        }
    }

    forAll(gammaExt().boundaryField(), patchI)
    {
        if (gammaExt().boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gammaExt().boundaryField()[patchI].patchInternalField();

            scalarField gammaNgb =
                gammaExt().boundaryField()[patchI].patchNeighbourField();

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll(gammaOwn, faceI)
            {
                if
                (
                    mag(gammaNgb[faceI] - gammaOwn[faceI]) > SMALL
                )
                {
                    if (!ibInsideFacesSet.found(start + faceI))
                    {
                        ibInsideFacesSet.insert(start + faceI);
                    }
                    
                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibInsideFacesSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibInsideFacesSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    insideFacesPtr_ = new labelList(ibInsideFacesSet.toc());

//     Info << "nIbInsideFaces: " << insideFacesPtr_->size() << endl;
}


void immersedBoundary::makeInternalFaces() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeInternalFaces() const : "
            << "create list of faces next to immersed boundary"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (internalFacesPtr_)
    {
        FatalErrorIn("immersedBoundary::makeInternalFaces() const")
            << "list of faces next to immersed boundary already exist"
            << abort(FatalError);
    }

    labelHashSet ibInternalFacesSet;

    const unallocLabelList& owner = mesh_.owner(); 
    const unallocLabelList& neighbour = mesh_.neighbour();

    volScalarField gammaTmp = gammaExt() - gamma();
    gammaTmp.correctBoundaryConditions();
    const scalarField& gammaTmpI = gammaTmp.internalField();

    forAll(neighbour, faceI)
    { 
        if
        (
            (gammaTmpI[neighbour[faceI]] > SMALL)
         && (gammaTmpI[owner[faceI]] > SMALL)
        )
        {
            ibInternalFacesSet.insert(faceI);
        }
    }

    forAll(gammaTmp.boundaryField(), patchI)
    {
        if (gammaTmp.boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gammaTmp.boundaryField()[patchI].patchInternalField();

            scalarField gammaNgb =
                gammaTmp.boundaryField()[patchI].patchNeighbourField();

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll(gammaOwn, faceI)
            {
                if
                (
                    (gammaNgb[faceI] > SMALL)
                 && (gammaOwn[faceI] > SMALL)
                )
                {
                    if (!ibInternalFacesSet.found(start + faceI))
                    {
                        ibInternalFacesSet.insert(start + faceI);
                    }
                    
                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibInternalFacesSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibInternalFacesSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    internalFacesPtr_ = new labelList(ibInternalFacesSet.toc());

//     Info << "nIbInternalFaces: " << internalFacesPtr_->size() << endl;
}


void immersedBoundary::makePointsAndNormals() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeIbPointsAndNormals() const : "
            << "create immersed  boundary points and normals"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointsPtr_ || normalsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeIbPointsAndNormals() const")
            << "immersed boundary points and normals already exist"
            << abort(FatalError);
    }

    // Find average cell dimension
    scalarField delta(cells().size(), 0.0);

    forAll(delta, cellI)
    {
        delta[cellI] = cellSize(cells()[cellI]);
    }

    // Find nearest triSurface point for each interface cell centre
    pointsPtr_ = new vectorField(cells().size(), vector::zero);
    normalsPtr_ = new vectorField(cells().size(), vector::zero);
    hitFacesPtr_ = new labelList(cells().size(), -1);

    vectorField& ibPoints = *pointsPtr_;
    vectorField& ibNormals = *normalsPtr_;
    labelList& ibHitFaces = *hitFacesPtr_;

    const vectorField& C = mesh_.cellCentres();

    forAll(cells(), cellI)
    {
        vector span
        (
            2*delta[cellI], 
            2*delta[cellI], 
            2*delta[cellI]
        );

        pointIndexHit pih = 
            triSurfSearch().nearest(C[cells()[cellI]], span);

        if(pih.hit())
        {
            ibPoints[cellI] = pih.hitPoint();
            ibNormals[cellI] =
                triSurfaceTools::surfaceNormal
                (
                    *this,
                    pih.index(),
                    pih.hitPoint()
                );

            if (internalFlow_)
            {
                ibNormals[cellI] *= -1;
            }

            ibHitFaces[cellI] = pih.index();
        }
        else
        {
            FatalErrorIn("immersedBoundary::makeIbPointsAndNormals() const")
                << "Can't find nearest triSurface point for cell "
                    << cells()[cellI] << ", " 
                    << mesh_.cellCentres()[cells()[cellI]]
                    << abort(FatalError);;
        }

//         Info << cellI << ", " << ibPoints[cellI] << endl;

        if
        (
            mesh_.nGeometricD() < 3 
         && mag(C[cells()[cellI]].z() - ibPoints[cellI].z()) > SMALL
        )
        {
            Info << C[cells()[cellI]].z() << endl;
            Info << ibPoints[cellI].z() << endl;
            Info << ibNormals[cellI] << endl;

            Info << "deltaZ = " 
                << mag(C[cells()[cellI]].z() - ibPoints[cellI].z())
                << endl;
            FatalErrorIn("immersedBoundary::makeIbPointsAndNormals() const")
                << "Intersection point is not on simmetry plane "
                    << "for 2-D geometry" << abort(FatalError);
        }
    }
}


void immersedBoundary::makeCellCells() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeCellCellsExt() const : "
            << "create neighbour cells for ib cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellCellsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeCellCells() const")
            << "cell-cell addressing already exists"
            << abort(FatalError);
    }

    cellCellsPtr_ = new labelListList(cells().size());
    labelListList& cellCells = *cellCellsPtr_;

    if (procCentresPtr_)
    {
        FatalErrorIn("immersedBoundary::makeCellCellsExt() const")
            << "procCentres already exists"
                << abort(FatalError);
    }
    procCentresPtr_ = new FieldField<Field, vector>(Pstream::nProcs());
    FieldField<Field, vector>& procCentres = *procCentresPtr_;
    forAll(procCentres, procI)
    {
        procCentres.set(procI, new vectorField(0));
    }


    if (cellProcCellsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeCellCellsExt() const")
            << "cellProcCells addressing already exists"
                << abort(FatalError);
    }
    cellProcCellsPtr_ = new List<List<labelPair> >(cells().size());
    List<List<labelPair> >& cellProcCells = *cellProcCellsPtr_;

    // 
    const vectorField& C = mesh_.cellCentres();

    scalarField rM = cellSizes();
    rM = 3.5*rM;

//     forAll(cellCellsExt, cellI) 
//     { 
//         findCellCells 
//         ( 
//             C[cells()[cellI]],  
//             cells()[cellI],  
//             cellCellsExt[cellI] 
//         ); 
 
//         for (label i = 0; i < cellCellsExt[cellI].size(); i++) 
//         { 
//             label curCell = cellCellsExt[cellI][i]; 
//             scalar r = mag(C[curCell] - C[cells()[cellI]]); 
//             if (r > rM[cellI]) 
//             { 
//                 cellCellsExt[cellI].setSize(i); 
//                 break; 
//             } 
//         } 
//     } 
 
    forAll(cellCells, cellI)
    {
        labelList curCells;

        findCellCells
        (
            C[cells()[cellI]], 
            cells()[cellI],
            curCells
        );

        cellCells[cellI] = labelList(curCells.size(), -1);
        label cI = 0;
        for (label i = 0; i < curCells.size(); i++)
        {
            label curCell = curCells[i];
            scalar r = mag(C[curCell] - C[cells()[cellI]]);
            if (r <= rM[cellI])
            {
                scalar limit = -::cos(80*M_PI/180.0);
                vector dir = (C[curCell]-points()[cellI]);
                dir /= mag(dir) + SMALL;
                if ((normals()[cellI]&dir)>=limit)
                {
                    cellCells[cellI][cI++] = curCell;
                }
            }
        }
        cellCells[cellI].setSize(cI);
    }

    if (Pstream::parRun())
    {
        // Find immersed boundary cells next to processor boundaries
        labelHashSet procIbCellsSet;

        forAll(cells(), cellI)
        {
            const labelList& curCellCells = cellCells[cellI];

            if (curCellCells.size())
            {
                forAll(curCellCells, cI)
                {
                    const labelList& faces = mesh_.cells()[curCellCells[cI]];

                    bool foundProcessorFace = false;

                    forAll(faces, faceI)
                    {
                        label patchID =
                            mesh_.boundaryMesh().whichPatch(faces[faceI]);

                        if (patchID != -1)
                        {
                            if
                            (
                                mesh_.boundaryMesh()[patchID].type()
                             == processorPolyPatch::typeName
                            )
                            {
                                foundProcessorFace = true;
                            }
                        }
                    }

                    if (foundProcessorFace)
                    {
                        procIbCellsSet.insert(cellI);
                        break;
                    }
                }
            }
            else
            {
                const labelList& faces = mesh_.cells()[cells()[cellI]];

                bool foundProcessorFace = false;

                forAll(faces, faceI)
                {
                    label patchID =
                        mesh_.boundaryMesh().whichPatch(faces[faceI]);

                    if (patchID != -1)
                    {
                        if
                        (
                            mesh_.boundaryMesh()[patchID].type()
                         == processorPolyPatch::typeName
                        )
                        {
                            foundProcessorFace = true;
                        }
                    }
                }

                if (foundProcessorFace)
                {
                    procIbCellsSet.insert(cellI);
                }
            }
        }
        labelList procIbCells = procIbCellsSet.toc();

        // Send and receive num of immersed boundary cells 
        // next to processor boundaries
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        sizeof(label)
                    );

                    toProc << procIbCells.size();
                }
            }
        }

        labelList sizes(Pstream::nProcs(), 0);
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI, 
                        sizeof(label)
                    );

                    fromProc >> sizes[procI];
                }
            }
        }

        // Send and receive ibCells centres and radii
        vectorField centres(procIbCells.size(), vector::zero);
        forAll(centres, cellI)
        {
            centres[cellI] = C[cells()[procIbCells[cellI]]];
        }
        scalarField rMax(rM, procIbCells);

        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        centres.size()*sizeof(vector)
                      + rMax.size()*sizeof(scalar)
                    );

                    toProc << centres << rMax;
                }
            }
        }

        FieldField<Field, vector> centres_(Pstream::nProcs());
        FieldField<Field, scalar> rMax_(Pstream::nProcs());

        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                centres_.set
                (
                    procI, 
                    new vectorField(sizes[procI], vector::zero)
                );

                rMax_.set
                (
                    procI, 
                    new scalarField(sizes[procI], 0)
                );

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI, 
                        sizes[procI]*sizeof(vector)
                      + sizes[procI]*sizeof(scalar)
                    );

                    fromProc >> centres_[procI] >> rMax_[procI];
                }
            }
            else
            {
                centres_.set(procI, new vectorField(0));
                rMax_.set(procI, new scalarField(0));
            }
        }

        // Find cells needed by other processors
        if (procCellsPtr_)
        {
            FatalErrorIn("immersedBoundary::makeCellCells() const")
                << "procCells addressing already exists"
                    << abort(FatalError);
        }
        procCellsPtr_ = new labelListList(Pstream::nProcs());
        labelListList& procCells = *procCellsPtr_;

        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            labelHashSet procCellSet;

            if (procI != Pstream::myProcNo())
            {
                forAll(centres_[procI], cellI)
                {
                    label nearestCellID =
                        findNearestCell(centres_[procI][cellI]);

                    if (nearestCellID == -1)
                    {
                        FatalErrorIn
                        (
                            "immersedBoundary::makeCellCells() const"
                        ) << "Can't find nearest cell."
                            << abort(FatalError);
                    }

                    scalar R = mag(C[nearestCellID] - centres_[procI][cellI]);

                    if (R < rMax_[procI][cellI])
                    {
                        if (!procCellSet.found(nearestCellID))
                        {
                            procCellSet.insert(nearestCellID);
                        }

                        labelList tmpCellList;
                        findCellCells
                        (
                            centres_[procI][cellI],
                            nearestCellID, 
                            tmpCellList
                        );
                        forAll(tmpCellList, cI)
                        {
                            scalar r = 
                                mag
                                (
                                    C[tmpCellList[cI]]
                                  - centres_[procI][cellI]
                                );

                            if (r <= rMax_[procI][cellI])
                            {
                                if (!procCellSet.found(tmpCellList[cI]))
                                {
                                    procCellSet.insert(tmpCellList[cI]);
                                }
                            }
                        }
                    }
                }
            }

            procCells[procI] = procCellSet.toc();
        }


        // Send and receive sizes
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        sizeof(label)
                    );

                    toProc << procCells[procI].size();
                }
            }
        }

        labelList procSizes(Pstream::nProcs(), 0);
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI, 
                        sizeof(label)
                    );

                    fromProc >> procSizes[procI];
                }
            }
        }
        

        // Send cell centres
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField centres(C, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        centres.size()*sizeof(vector)
                    );

                    toProc << centres;
                }
            }
        }

        // Receive cell centres
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procCentres.set
                (
                    procI,
                    new vectorField(procSizes[procI], vector::zero)
                );

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI, 
                        procSizes[procI]*sizeof(vector)
                    );

                    fromProc >> procCentres[procI];
                }
            }
            else
            {
                procCentres.set(procI, new vectorField(0));
            }
        }

        // Cell-procCells addressing
        forAll(cellProcCells, cellI)
        {
            scalar rMax = rM[cellI];

            cellProcCells[cellI].setSize(100);

            label index = 0;
            forAll(procCentres, procI)
            {
                forAll(procCentres[procI], pointI)
                {
                    scalar r = 
                        mag
                        (
                            procCentres[procI][pointI]
                          - C[cells()[cellI]]
                        );

                    if (r <= rMax)
                    {
                        cellProcCells[cellI][index].first() = procI;
                        cellProcCells[cellI][index].second() = pointI;
                        index++;
                    }
                }
            }

            cellProcCells[cellI].setSize(index);
        }
    }
}


void immersedBoundary::makeDeadCells() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeDeadCells() const : "
            << "create list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadCellsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeDeadCells() const")
            << "list of dead cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaExtI = gammaExt().internalField();

    deadCellsPtr_ = new labelList(label(sum(1.0 - gammaExtI)), -1);
    labelList& deadCells = *deadCellsPtr_;

    label counter = 0;
    forAll(gammaExtI, cellI)
    {
        if(gammaExtI[cellI] < SMALL)
        {
            deadCells[counter++] = cellI;
        }
    }
}


void immersedBoundary::makeDeadCellsExt() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeDeadCellsExt() const : "
            << "create extended list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadCellsExtPtr_)
    {
        FatalErrorIn("immersedBoundary::makeDeadCellsExt() const")
            << "extended list of dead cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaI = gamma().internalField();

    deadCellsExtPtr_ = new labelList(label(sum(1.0 - gammaI)), -1);
    labelList& deadCellsExt = *deadCellsExtPtr_;

    label counter = 0;
    forAll(gammaI, cellI)
    {
        if(gammaI[cellI] < SMALL)
        {
            deadCellsExt[counter++] = cellI;
        }
    }
}


void immersedBoundary::makeLiveCells() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeLiveCells() const : "
            << "create list of live cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (liveCellsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeLiveCells() const")
            << "list of live cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaI = gamma().internalField();

    liveCellsPtr_ = new labelList(label(sum(gammaI)), -1);
    labelList& liveCells = *liveCellsPtr_;

    label counter = 0;
    forAll(gammaI, cellI)
    {
        if(gammaI[cellI] > (1.0 - SMALL))
        {
            liveCells[counter++] = cellI;
        }
    }
}


void immersedBoundary::makeCellSizes() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeIbCellsSize() const : "
            << "create average sizes of immersed boundary cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellSizesPtr_)
    {
        FatalErrorIn("immersedBoundary::makeIbCellsSize() const")
            << "average sizes of immersed boundary cells already exist"
            << abort(FatalError);
    }

    cellSizesPtr_ = new scalarField(points().size(), 0.0);
    scalarField& delta = *cellSizesPtr_;

    if(mesh_.nGeometricD() == 3)
    {
        scalarField V(mesh_.V().field(), cells());

        delta = Foam::pow(V, 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = sqrt(scalarField(mesh_.V().field(), cells())/thickness);
    }
}


void immersedBoundary::makeDelta() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeDelta() const : "
            << "creating delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deltaPtr_)
    {
        FatalErrorIn("immersedBoundary::makeDelta() const")
            << "delta field already exist"
            << abort(FatalError);
    }

    const vectorField C = mesh_.cellCentres();

    deltaPtr_ = new scalarField(mag(vectorField(C, cells()) - points()));
}


// Find the cell with the nearest cell centre
label immersedBoundary::findNearestCell(const point& location) const
{
    const vectorField& centres = mesh_.cellCentres();
    const scalarField& gammaExtI = gammaExt().internalField();

    label nearestCellI = -1;
    scalar minProximity = GREAT;

    for (label cellI = 0; cellI < centres.size(); cellI++)
    {    
        if (gammaExtI[cellI] > SMALL)
        {
            scalar proximity = magSqr(centres[cellI] - location);

            if (proximity < minProximity)
            {
                nearestCellI = cellI;
                minProximity = proximity;
            }
        }
    }

    return nearestCellI;
}


void immersedBoundary::findCellCells
(
    const vector& pt,
    const label cellID, 
    labelList& cellCells
) const
{
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();

    const scalarField& gammaExtI = gammaExt().internalField();

    labelHashSet cellSet;
    cellSet.insert(cellID);

    // First row
    const labelList& curCellPoints = cellPoints[cellID];
    forAll(curCellPoints, pointI)
    {
        label curPoint = curCellPoints[pointI];
        const labelList& curPointCells = pointCells[curPoint];

        forAll(curPointCells, cI)
        {
            if (gammaExtI[curPointCells[cI]] > SMALL)
            {
                if (!cellSet.found(curPointCells[cI]))
                {
                    cellSet.insert(curPointCells[cI]);
                }
            }
        }
    }

    // Second row
    labelList curCells = cellSet.toc();
    forAll(curCells, cellI)
    {
        label curCell = curCells[cellI];
        const labelList& curCellPoints = cellPoints[curCell];

        forAll(curCellPoints, pointI)
        {
            label curPoint = curCellPoints[pointI];
            const labelList& curPointCells = pointCells[curPoint];

            forAll(curPointCells, cI)
            {
                if (gammaExtI[curPointCells[cI]] > SMALL)
                {
                    if (!cellSet.found(curPointCells[cI]))
                    {
                        cellSet.insert(curPointCells[cI]);
                    }            
                }
            }
        }
    }

    // Third row
    curCells = cellSet.toc();
    forAll(curCells, cellI)
    {
        label curCell = curCells[cellI];
        const labelList& curCellPoints = cellPoints[curCell];

        forAll(curCellPoints, pointI)
        {
            label curPoint = curCellPoints[pointI];
            const labelList& curPointCells = pointCells[curPoint];

            forAll(curPointCells, cI)
            {
                if (gammaExtI[curPointCells[cI]] > SMALL)
                {
                    if (!cellSet.found(curPointCells[cI]))
                    {
                        cellSet.insert(curPointCells[cI]);
                    }            
                }
            }
        }
    }

//     // Fourth row
//     curCells = cellSet.toc();
//     forAll(curCells, cellI)
//     {
//         label curCell = curCells[cellI];
//         const labelList& curCellPoints = cellPoints[curCell];

//         forAll(curCellPoints, pointI)
//         {
//             label curPoint = curCellPoints[pointI];
//             const labelList& curPointCells = pointCells[curPoint];

//             forAll(curPointCells, cI)
//             {
//                 if (gammaExtI[curPointCells[cI]] > SMALL)
//                 {
//                     if (!cellSet.found(curPointCells[cI]))
//                     {
//                         cellSet.insert(curPointCells[cI]);
//                     }            
//                 }
//             }
//         }
//     }

    // Erase current cell
    cellSet.erase(cellID);

    // Sorting cells
    const vectorField& C = mesh_.cellCentres();
    curCells = cellSet.toc();
    scalarField distances(curCells.size(), 0);
    forAll(distances, cI)
    {
        distances[cI] =
            mag(C[curCells[cI]] - pt);
    }

    SortableList<scalar> sortedDistances(distances);
    labelList sortedCells(curCells.size(), -1);
    for (label i=0; i<sortedCells.size(); i++)
    {
        sortedCells[i] =
            curCells[sortedDistances.indices()[i]];
    }

    cellCells = sortedCells;
}


scalar immersedBoundary::cellSize(label cellID) const
{
    scalar delta;

    if(mesh_.nGeometricD() == 3)
    {
        delta = Foam::pow(mesh_.V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = sqrt(mesh_.V().field()[cellID]/thickness);
    }

    return delta;
}


scalar immersedBoundary::cellProjection
(
    label cellID,
    const vector& dir
) const
{
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const labelList& cellFaces = mesh_.cells()[cellID];

    scalar area = 0;

    forAll(cellFaces, faceI)
    {
        label curFace = cellFaces[faceI];

        vector curSf = Sf[curFace];

        if ((curSf&(Cf[curFace] - C[cellID])) < 0)
        {
            curSf *= -1;
        }

        if ((curSf&dir) > 1)
        {
            area += (curSf&dir);
        }
    }

    return area;
}


void immersedBoundary::clearOut()
{
    deleteDemandDrivenData(triSurfSearchPtr_);
    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(gammaExtPtr_);
    deleteDemandDrivenData(sGammaPtr_);
    deleteDemandDrivenData(cellsPtr_);
    deleteDemandDrivenData(facesPtr_);
    deleteDemandDrivenData(insideFacesPtr_);
    deleteDemandDrivenData(internalFacesPtr_);
    deleteDemandDrivenData(pointsPtr_);
    deleteDemandDrivenData(normalsPtr_);
    deleteDemandDrivenData(hitFacesPtr_);
    deleteDemandDrivenData(cellCellsPtr_);
    deleteDemandDrivenData(procCellsPtr_);
    deleteDemandDrivenData(procCentresPtr_);
    deleteDemandDrivenData(cellProcCellsPtr_);
    deleteDemandDrivenData(deadCellsPtr_);
    deleteDemandDrivenData(deadCellsExtPtr_);
    deleteDemandDrivenData(liveCellsPtr_);
    deleteDemandDrivenData(cellSizesPtr_);
    invDirichletMatrices_.clear();
    invNeumannMatrices_.clear();
    deleteDemandDrivenData(refinementCellsPtr_);
    deleteDemandDrivenData(wallDistancePtr_);
    deleteDemandDrivenData(deltaPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundary::immersedBoundary
(
    const fvMesh& mesh
)
:
    triSurfaceMesh
    (
        IOobject
        (
            triSurface::typeName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        triSurface
        (
            mesh.time().path()/mesh.time().caseConstant()
           /triSurface::typeName/"immersedBoundary.ftr"
        )
    ),
    mesh_(mesh),
    internalFlow_(false),
    triSurfSearchPtr_(NULL),
    gammaPtr_(NULL),
    gammaExtPtr_(NULL),
    sGammaPtr_(NULL),
    cellsPtr_(NULL),
    facesPtr_(NULL),
    insideFacesPtr_(NULL),
    internalFacesPtr_(NULL),
    pointsPtr_(NULL),
    normalsPtr_(NULL),
    hitFacesPtr_(NULL),
    cellCellsPtr_(NULL),
    procCellsPtr_(NULL),
    procCentresPtr_(NULL),
    cellProcCellsPtr_(NULL),
    deadCellsPtr_(NULL),
    deadCellsExtPtr_(NULL),
    liveCellsPtr_(NULL),
    cellSizesPtr_(NULL),
    invDirichletMatrices_(),
    invNeumannMatrices_(),
    refinementCellsPtr_(NULL),
    wallDistancePtr_(NULL),
    deltaPtr_(NULL)
{
    IOdictionary ibProperties
    (
        IOobject
        (
            "immersedBoundaryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    internalFlow_ = Switch(ibProperties.lookup("internalFlow"));
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

immersedBoundary::~immersedBoundary()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::immersedBoundary::movePoints(const pointField& newPoints)
{
    clearOut();

    triSurfaceMesh::movePoints(newPoints);
}


tmp<volVectorField> immersedBoundary::grad(const volScalarField& p) const
{
    tmp<volVectorField> tGradP
    (
        new volVectorField
        (
            IOobject
            (
                "grad(" + p.name() + ')',
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", p.dimensions()/dimLength, vector::zero)
        )
    );
    volVectorField& gradP = tGradP();

    surfaceScalarField pf = fvc::interpolate(p);

//     const volVectorField& prevGradP = 
//         mesh_.lookupObject<volVectorField>("grad(" + p.name() + ')');

    // Skew-correction
    if (skewCorrectionVectors::New(this->mesh_).skew())
    {
        const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh_);

        pf += 
        (
            scv()
          & linear<vector>(mesh_).interpolate
            (
                fv::leastSquaresGrad<scalar>(mesh_).grad(p)
            )
        );
    }

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();

    label iCorr = 0;

    do
    {
        iCorr++;

        // Interface correction
        forAll(faces(), faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh_.nInternalFaces())
            {
                label liveCell = -1;
                if (gamma()[owner[curFace]] > SMALL)
                {
                    liveCell = owner[curFace];
                }
                else
                {
                    liveCell = neighbour[curFace];
                }

                vector delta = Cf[curFace] - C[liveCell];

                pf.internalField()[curFace] =
                    p[liveCell] + (delta&gradP[liveCell]);
            }
            else
            {
                label patchID = mesh_.boundaryMesh().whichPatch(curFace);
                label start = mesh_.boundaryMesh()[patchID].start();
                label localIndex = curFace - start;

                const unallocLabelList& faceCells =
                    mesh_.boundary()[patchID].faceCells();

                if (gamma()[faceCells[localIndex]] > SMALL)
                {
                    vector delta = Cf[curFace] - C[faceCells[localIndex]];

                    pf.boundaryField()[patchID][localIndex] =
                        p[faceCells[localIndex]] 
                      + (delta&gradP[faceCells[localIndex]]); 
                }
            }
        }
        
        // Gradient calculation using Gauss method
        gradP = fv::gaussGrad<scalar>(mesh_).grad(pf);
        fv::gaussGrad<scalar>(mesh_).correctBoundaryConditions(p, gradP);
    }
    while(iCorr < 3);

    return tGradP;
}


const triSurfaceSearch& immersedBoundary::triSurfSearch() const
{
    if (!triSurfSearchPtr_)
    {
        makeTriSurfSearch();
    }

    return *triSurfSearchPtr_;
}


const volScalarField& immersedBoundary::gamma() const
{
    if (!gammaPtr_)
    {
        makeGamma();
    }

    return *gammaPtr_;
}


const volScalarField& immersedBoundary::gammaExt() const
{
    if (!gammaExtPtr_)
    {
        makeGammaExt();
    }

    return *gammaExtPtr_;
}


const surfaceScalarField& immersedBoundary::sGamma() const
{
    if (!sGammaPtr_)
    {
        makeSGamma();
    }

    return *sGammaPtr_;
}


const labelList& immersedBoundary::cells() const
{
    if (!cellsPtr_)
    {
        makeCells();

//         addCornerCells();
    }

    return *cellsPtr_;
}


const labelList& immersedBoundary::faces() const
{
    if (!facesPtr_)
    {
        makeFaces();
    }

    return *facesPtr_;
}


const labelList& immersedBoundary::insideFaces() const
{
    if (!insideFacesPtr_)
    {
        makeInsideFaces();
    }

    return *insideFacesPtr_;
}


const labelList& immersedBoundary::internalFaces() const
{
    if (!internalFacesPtr_)
    {
        makeInternalFaces();
    }

    return *internalFacesPtr_;
}


const vectorField& immersedBoundary::points() const
{
    if (!pointsPtr_)
    {
        makePointsAndNormals();
    }

    return *pointsPtr_;
}


const vectorField& immersedBoundary::normals() const
{
    if (!normalsPtr_)
    {
        makePointsAndNormals();
    }

    return *normalsPtr_;
}


const labelList& immersedBoundary::hitFaces() const
{
    if (!hitFacesPtr_)
    {
        makePointsAndNormals();
    }

    return *hitFacesPtr_;
}


const labelListList& immersedBoundary::cellCells() const
{
    if (!cellCellsPtr_)
    {
        makeCellCells();
    }

    return *cellCellsPtr_;
}


const FieldField<Field, vector>& immersedBoundary::procCentres() const
{
    if (!procCentresPtr_)
    {
        makeCellCells();
    }

    return *procCentresPtr_;
}

const List<List<labelPair> >& immersedBoundary::cellProcCells() const
{
    if (!cellProcCellsPtr_)
    {
        makeCellCells();
    }

    return *cellProcCellsPtr_;
}

const labelListList& immersedBoundary::procCells() const
{
    if (!procCellsPtr_)
    {
        makeCellCells();
    }

    return *procCellsPtr_;
}

const labelList& immersedBoundary::deadCells() const
{
    if (!deadCellsPtr_)
    {
        makeDeadCells();
    }

    return *deadCellsPtr_;
}

const labelList& immersedBoundary::deadCellsExt() const
{
    if (!deadCellsExtPtr_)
    {
        makeDeadCellsExt();
    }

    return *deadCellsExtPtr_;
}

const labelList& immersedBoundary::liveCells() const
{
    if (!liveCellsPtr_)
    {
        makeLiveCells();
    }

    return *liveCellsPtr_;
}

const scalarField& immersedBoundary::cellSizes() const
{
    if (!cellSizesPtr_)
    {
        makeCellSizes();
    }

    return *cellSizesPtr_;
}

const PtrList<scalarRectangularMatrix>& 
immersedBoundary::invDirichletMatrices() const
{
    label size = invDirichletMatrices_.size();

    reduce(size, maxOp<label>());

    if (size == 0)
    {
        makeInvDirichletMatrices();
    }

    return invDirichletMatrices_;
}

const PtrList<scalarRectangularMatrix>& 
immersedBoundary::invNeumannMatrices() const
{
    label size = invNeumannMatrices_.size();

    reduce(size, maxOp<label>());

    if (size == 0)
    {
        makeInvNeumannMatrices();
    }

    return invNeumannMatrices_;
}

const scalarField& immersedBoundary::delta() const
{
    if (!deltaPtr_)
    {
        makeDelta();
    }

    return *deltaPtr_;
}

bool immersedBoundary::write() const
{
    if (mesh_.time().outputTime())
    {
        if(!isDir(mesh_.time().timePath()/triSurface::typeName))
        {
            mkDir(mesh_.time().timePath()/triSurface::typeName);
        }

        fileName triSurfFileName  =
            mesh_.time().timePath()/triSurface::typeName
           /("immersedBoundary.ftr");

        triSurface::write
        (
            triSurfFileName, 
            false
        );
    }

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
