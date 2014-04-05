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
#include "plane.H"
#include "wedgePolyPatch.H"
#include "multiDirRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void immersedBoundary::makeRefinementCells() const
{
    if (debug)
    {
        Info<< "void immersedBoundary::makeRefinementCells() const : "
            << "create list of cells for refinement"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (refinementCellsPtr_)
    {
        FatalErrorIn("immersedBoundary::makeRefinementCells() const")
            << "list of refinement cells already exist"
            << abort(FatalError);
    }

    labelHashSet refCellSet;

//     const labelListList& cellCells = cellCellsExt();

    forAll(cells(), cellI)
    {
        refCellSet.insert(cells()[cellI]);
    }

    forAll(cellCells(), cellI)
    {
        const labelList& curCells = cellCells()[cellI];

        forAll(curCells, cI)
        {
            if (!refCellSet.found(curCells[cI]))
            {
                refCellSet.insert(curCells[cI]);
            }            
        }
    }

    const unallocLabelList& owner = mesh_.owner(); 
    const unallocLabelList& neighbour = mesh_.neighbour();

//     forAll(faces(), faceI)
//     {
//         const label& ownCell = owner[faces()[faceI]];
//         const label& ngbCell = neighbour[faces()[faceI]];

//         if(!refCellSet.found(ownCell))
//         {
//             refCellSet.insert(ownCell);
//         }

//         if(!refCellSet.found(ngbCell))
//         {
//             refCellSet.insert(ngbCell);
//         }
//     }

    const scalarField& gammaExtI = gammaExt().internalField();

    forAll(insideFaces(), faceI)
    {
        const label& ownCell = owner[insideFaces()[faceI]];
        const label& ngbCell = neighbour[insideFaces()[faceI]];

        if (gammaExtI[ownCell] < SMALL)
        {
            if (!refCellSet.found(ownCell))
            {
                refCellSet.insert(ownCell);
            }
        }
        else
        {
            if (!refCellSet.found(ngbCell))
            {
                refCellSet.insert(ngbCell);
            }
        }
    }

    refinementCellsPtr_ = new labelList(refCellSet.toc());
}

// Return index of coordinate axis.
label immersedBoundary::axis(const vector& normal) const
{
    const scalar edgeTol = 1E-3;

    label axisIndex = -1;

    if (mag(normal & point(1, 0, 0)) > (1-edgeTol))
    {
        axisIndex = 0;
    }
    else if (mag(normal & point(0, 1, 0)) > (1-edgeTol))
    {
        axisIndex = 1;
    }
    else if (mag(normal & point(0, 0, 1)) > (1-edgeTol))
    {
        axisIndex = 2;
    }

    return axisIndex;
}


//- Returns -1 or cartesian coordinate component (0=x, 1=y, 2=z) of normal
//  in case of 2D mesh
label immersedBoundary::twoDNess(const polyMesh& mesh) const
{
    const pointField& ctrs = mesh.cellCentres();

    if (ctrs.size() < 2)
    {
        return -1;
    }


    //
    // 1. All cell centres on single plane aligned with x, y or z
    //

    // Determine 3 points to base plane on.
    vector vec10 = ctrs[1] - ctrs[0];
    vec10 /= mag(vec10);

    label otherCellI = -1;

    for (label cellI = 2; cellI < ctrs.size(); cellI++)
    {
        vector vec(ctrs[cellI] - ctrs[0]);
        vec /= mag(vec);

        if (mag(vec & vec10) < 0.9)
        {
            // ctrs[cellI] not in line with n
            otherCellI = cellI;

            break;
        }
    }

    if (otherCellI == -1)
    {
        // Cannot find cell to make decent angle with cell0-cell1 vector.
        // Note: what to do here? All cells (almost) in one line. Maybe 1D case?
        return -1;
    }

    plane cellPlane(ctrs[0], ctrs[1], ctrs[otherCellI]);


    forAll(ctrs, cellI)
    {
        const labelList& cEdges = mesh.cellEdges()[cellI];

        scalar minLen = GREAT;

        forAll(cEdges, i)
        {
            minLen = min(minLen, mesh.edges()[cEdges[i]].mag(mesh.points()));
        }

        if (cellPlane.distance(ctrs[cellI]) > 1E-6*minLen)
        {
            // Centres not in plane
            return  -1;
        }
    }

    label axisIndex = axis(cellPlane.normal());

    if (axisIndex == -1)
    {
        return axisIndex;
    }


    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    //
    // 2. No edges without points on boundary
    //

    // Mark boundary points
    boolList boundaryPoint(mesh.allPoints().size(), false);

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            forAll(f, fp)
            {
                boundaryPoint[f[fp]] = true;
            }
        }
    }


    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (!boundaryPoint[e.start()] && !boundaryPoint[e.end()])
        {
            // Edge has no point on boundary.
            return -1;
        }
    }


    // 3. For all non-wedge patches: all faces either perp or aligned with
    //    cell-plane normal. (wedge patches already checked upon construction)

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!isA<wedgePolyPatch>(patch))
        {
            const vectorField& n = patch.faceAreas();

            scalarField cosAngle = mag(n/mag(n) & cellPlane.normal());

            if (mag(min(cosAngle) - max(cosAngle)) > 1E-6)
            {
                // cosAngle should be either ~1 over all faces (2D front and
                // back) or ~0 (all other patches perp to 2D)
                return -1;
            }
        }
    }

    return axisIndex;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void immersedBoundary::refineMesh
(
    polyMesh& pMesh, 
    const labelList& refCells
) const
{
    Info << "nRefCells = " << refCells.size() << "\n" << endl;

    // Dictionary to control refinement
    dictionary refineDict;

    // Set refinement directions based on 2D/3D
    label axisIndex = twoDNess(pMesh);

    if (axisIndex == -1)
    {
        Info<< "3D case; refining all directions" << nl << endl;

        wordList directions(3);
        directions[0] = "tan1";
        directions[1] = "tan2";
        directions[2] = "normal";
        refineDict.add("directions", directions);

        // Use hex cutter
        refineDict.add("useHexTopology", "true");
    }
    else
    {
        wordList directions(2);

        if (axisIndex == 0)
        {
            Info<< "2D case; refining in directions y,z\n" << endl;
            directions[0] = "tan2";
            directions[1] = "normal";
        }
        else if (axisIndex == 1)
        {
            Info<< "2D case; refining in directions x,z\n" << endl;
            directions[0] = "tan1";
            directions[1] = "normal";
        }
        else
        {
            Info<< "2D case; refining in directions x,y\n" << endl;
            directions[0] = "tan1";
            directions[1] = "tan2";
        }

        refineDict.add("directions", directions);

        // Use standard cutter
        refineDict.add("useHexTopology", "false");
    }

    refineDict.add("coordinateSystem", "global");

    dictionary coeffsDict;
    coeffsDict.add("tan1", vector(1, 0, 0));
    coeffsDict.add("tan2", vector(0, 1, 0));
    refineDict.add("globalCoeffs", coeffsDict);

    refineDict.add("geometricCut", "false");
    refineDict.add("writeMesh", "false");

    // Multi-directional refinement (does multiple iterations)
    multiDirRefinement multiRef(pMesh, refCells, refineDict);
}

const labelList& immersedBoundary::refinementCells() const
{
    if (!refinementCellsPtr_)
    {
        makeRefinementCells();
    }

    return *refinementCellsPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
