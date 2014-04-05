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
    makeTriSurfaceMesh

Description
    A mesh generator for triSurface mesh from fvMesh patches.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
// #include "labelHashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    labelHashSet includedPatches;

    for(label patchI=0; patchI<mesh.boundaryMesh().size(); patchI++)
    {
        includedPatches.insert(patchI);
    }

    triSurface triSurf = 
        triSurfaceTools::triangulate
        (
            mesh.boundaryMesh(), 
            includedPatches
        );

    // Create triSurface folder
//     mkDir
//     (
//         runTime.path()
//        /triSurf.triSurfInstance(runTime)
//        /triSurface::typeName
//     );

    // Writing triSurface mesh
    Info << "Write triSurface mesh ... ";
//     triSurf.write(runTime);
    triSurf.write("immersedBoundary.stl");
    triSurf.write("immersedBoundary.ftr");


    Info << "Done" << endl;

    return(0);
}

// ************************************************************************* //
