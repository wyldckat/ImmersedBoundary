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

Class
    Foam::triSurfaceHandler

Description
    Triangulated surface description with patch information.

SourceFiles
    triSurfaceHandler.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceHandler_H
#define triSurfaceHandler_H

#include "PrimitivePatch.H"
#include "triSurface.H"
#include "pointField.H"
#include "labelledTri.H"
#include "boolList.H"
#include "geometricSurfacePatchList.H"
#include "surfacePatchList.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class IFstream;

/*---------------------------------------------------------------------------*\
                           Class triSurfaceHandler Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceHandler
:
    public triSurface
{
public:

    // Static

        // Adapted from OpenFOAM 1.6-ext
        //"src/OpenFOAM/meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
        static void writeVTKNormals
        (
            const fileName& name,
            const List<labelledTri>& faces,
            const Field<point>& points
        );


    // Constructors

        //- Construct null
        triSurfaceHandler();

        //- Construct from triangles, patches, points.
        triSurfaceHandler
        (
            const List<labelledTri>&,
            const geometricSurfacePatchList&,
            const pointField&
        );

        //- Construct from triangles, patches, points. Reuse storage.
        triSurfaceHandler
        (
            List<labelledTri>&,
            const geometricSurfacePatchList&,
            pointField&,
            const bool reUse
        );

        //- Construct from triangles, points. Set patchnames to default.
        triSurfaceHandler(const List<labelledTri>&, const pointField&);

        //- Construct from triangles, points. Set region to 0 and default
        //  patchName.
        triSurfaceHandler(const triFaceList&, const pointField&);

        //- Construct from file name (uses extension to determine type)
        triSurfaceHandler(const fileName&);

        //- Construct from Istream
        triSurfaceHandler(Istream&);

        //- Construct from objectRegistry
        triSurfaceHandler(const Time& d);

        //- Construct as copy
        triSurfaceHandler(const triSurfaceHandler&);



    // Destructor

        ~triSurfaceHandler();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
