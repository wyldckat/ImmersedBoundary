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

\*---------------------------------------------------------------------------*/

#include "refineImmersedBoundaryMesh.H"
#include "immersedBoundaryFvPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "plane.H"
#include "wedgePolyPatch.H"
#include "multiDirRefinement.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refineImmersedBoundaryMesh::makeRefinementCells() const
{
    if (fvMesh::debug)
    {
        InfoIn("void refineImmersedBoundaryMesh::makeRefinementCells() const")
            << "create list of cells for refinement"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (refinementCellsPtr_)
    {
        FatalErrorIn
        (
            "void refineImmersedBoundaryMesh::makeRefinementCells() const"
        )   << "list of refinement cells already exist"
            << abort(FatalError);
    }

    labelHashSet refCellSet;

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Insert immersed boundary cells from all immersed boundary patches
    forAll (mesh_.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh_.boundary()[patchI]))
        {
            Info<< "Found immersed boundary patch " << patchI
                << " named " << mesh_.boundary()[patchI].name()
                << endl;

            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh_.boundary()[patchI]
                );

            const labelList& c = ibPatch.ibCells();

            forAll (c, cI)
            {
                if (!refCellSet.found(c[cI]))
                {
                    refCellSet.insert(c[cI]);
                }
            }

            const labelListList& cc = ibPatch.ibCellCells();

            forAll (cc, cellI)
            {
                const labelList& curCells = cc[cellI];

                forAll (curCells, cI)
                {
                    if (!refCellSet.found(curCells[cI]))
                    {
                        refCellSet.insert(curCells[cI]);
                    }
                }
            }

            const scalarField& gammaExtI = ibPatch.gammaExt().internalField();

            const labelList& ibInsideFaces = ibPatch.ibInsideFaces();

            forAll (ibInsideFaces, faceI)
            {
                const label ownCell = owner[ibInsideFaces[faceI]];
                const label ngbCell = neighbour[ibInsideFaces[faceI]];

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
        }
    }

    refinementCellsPtr_ = new labelList(refCellSet.toc());
}


// Return index of coordinate axis.
Foam::label Foam::refineImmersedBoundaryMesh::axis(const vector& normal) const
{
    const scalar edgeTol = 1e-3;

    label axisIndex = -1;

    if (mag(normal & vector(1, 0, 0)) > (1 - edgeTol))
    {
        axisIndex = 0;
    }
    else if (mag(normal & vector(0, 1, 0)) > (1 - edgeTol))
    {
        axisIndex = 1;
    }
    else if (mag(normal & vector(0, 0, 1)) > (1 - edgeTol))
    {
        axisIndex = 2;
    }

    return axisIndex;
}


//- Returns -1 or cartesian coordinate component (0=x, 1=y, 2=z) of normal
//  in case of 2D mesh
Foam::label Foam::refineImmersedBoundaryMesh::twoDNess() const
{
    const pointField& ctrs = mesh_.cellCentres();

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
        // Note: what to do here? All cells (almost) in one line.
        // Maybe 1D case?
        return -1;
    }

    plane cellPlane(ctrs[0], ctrs[1], ctrs[otherCellI]);


    forAll (ctrs, cellI)
    {
        const labelList& cEdges = mesh_.cellEdges()[cellI];

        scalar minLen = GREAT;

        forAll (cEdges, i)
        {
            minLen = min(minLen, mesh_.edges()[cEdges[i]].mag(mesh_.points()));
        }

        if (cellPlane.distance(ctrs[cellI]) > 1e-6*minLen)
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


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    //
    // 2. No edges without points on boundary
    //

    // Mark boundary points
    boolList boundaryPoint(mesh_.allPoints().size(), false);

    forAll (patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        forAll (patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            forAll (f, fp)
            {
                boundaryPoint[f[fp]] = true;
            }
        }
    }


    const edgeList& edges = mesh_.edges();

    forAll (edges, edgeI)
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

    forAll (patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!isA<wedgePolyPatch>(patch))
        {
            const vectorField& n = patch.faceAreas();

            scalarField cosAngle = mag(n/mag(n) & cellPlane.normal());

            if (mag(min(cosAngle) - max(cosAngle)) > 1e-6)
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

Foam::refineImmersedBoundaryMesh::refineImmersedBoundaryMesh
(
    fvMesh& mesh
)
:
    mesh_(mesh),
    refinementCellsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::refineImmersedBoundaryMesh::~refineImmersedBoundaryMesh()
{
    deleteDemandDrivenData(refinementCellsPtr_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::labelList&
Foam::refineImmersedBoundaryMesh::refinementCells() const
{
    if (!refinementCellsPtr_)
    {
        makeRefinementCells();
    }

    return *refinementCellsPtr_;
}


void Foam::refineImmersedBoundaryMesh::refineMesh
(
    const labelList& refCells
) const
{
    Info << "nRefCells = " << refCells.size() << "\n" << endl;

    // Dictionary to control refinement
    dictionary refineDict;

    // Set refinement directions based on 2D/3D
    label axisIndex = twoDNess();

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
    multiDirRefinement multiRef(mesh_, refCells, refineDict);
}


// ************************************************************************* //
