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

#include "triSurfaceHandler.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "boundBox.H"
#include "SortableList.H"
#include "PackedBoolList.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceHandler::triSurfaceHandler()
:
    triSurface()
{}



Foam::triSurfaceHandler::triSurfaceHandler
(
    const List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const pointField& points
)
:
    triSurface(triangles, patches, points)
{}


Foam::triSurfaceHandler::triSurfaceHandler
(
    List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    pointField& points,
    const bool reUse
)
:
    triSurface(triangles, patches, points, reUse)
{}


Foam::triSurfaceHandler::triSurfaceHandler
(
    const List<labelledTri>& triangles,
    const pointField& points
)
:
    triSurface(triangles, points)
{
}


Foam::triSurfaceHandler::triSurfaceHandler
(
    const triFaceList& triangles,
    const pointField& points
)
:
    triSurface(triangles, points)
{
}


Foam::triSurfaceHandler::triSurfaceHandler(const fileName& name)
:
    triSurface(name)
{
}


Foam::triSurfaceHandler::triSurfaceHandler(Istream& is)
:
    triSurface(is)
{
}


Foam::triSurfaceHandler::triSurfaceHandler(const Time& d)
:
    triSurface(d)
{
}


Foam::triSurfaceHandler::triSurfaceHandler(const triSurfaceHandler& ts)
:
    triSurface(ts)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceHandler::~triSurfaceHandler()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Adapted from OpenFOAM 1.6-ext
//"src/OpenFOAM/meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
void
Foam::triSurfaceHandler::
writeVTKNormals
(
    const fileName& name,
    const List<labelledTri>& faces,
    const Field<point>& points
)
{
    // Write patch and points into VTK
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << faces.size() << " float" << nl;

    // Write points
    List<float> mlPointBuffer(3*faces.size());

    label counter = 0;
    forAll (faces, i)
    {
        const vector c = faces[i].centre(points);

        mlPointBuffer[counter++] = float(c.x());
        mlPointBuffer[counter++] = float(c.y());
        mlPointBuffer[counter++] = float(c.z());
    }

    forAll (mlPointBuffer, i)
    {
        mps << mlPointBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;

    // Write normals
    mps << "POINT_DATA " << faces.size() << nl
        << "FIELD attributes " << 1 << nl
        << "normals" << " 3 "
        << faces.size() << " float" << nl;

    List<float> mlNormalBuffer(3*faces.size());

    counter = 0;
    forAll (faces, i)
    {
        const vector n = faces[i].normal(points);

        mlNormalBuffer[counter++] = float(n.x());
        mlNormalBuffer[counter++] = float(n.y());
        mlNormalBuffer[counter++] = float(n.z());
    }

    forAll (mlNormalBuffer, i)
    {
        mps << mlNormalBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;
}

// ************************************************************************* //
