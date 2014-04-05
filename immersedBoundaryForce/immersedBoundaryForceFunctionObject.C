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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryForceFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryForceFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        immersedBoundaryForceFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryForceFunctionObject::
immersedBoundaryForceFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    forcesFilePtr_(NULL),
    Uref_(readScalar(dict.lookup("Uref"))),
    Aref_(readScalar(dict.lookup("Aref")))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info << "Creating immersed boundary force function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::immersedBoundaryForceFunctionObject::start()
{
    write();

    return true;
}


bool Foam::immersedBoundaryForceFunctionObject::execute()
{
    write();

    return true;
}


bool Foam::immersedBoundaryForceFunctionObject::read
(
    const dictionary& dict
)
{
    return false;
}

void Foam::immersedBoundaryForceFunctionObject::makeFile()
{
    // Create the forces file if not already created
    if (forcesFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating forces file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName forcesDir;
            word startTimeName =
                time_.timeName(time_.startTime().value());

            word name_("ibForce");

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                forcesDir = time_.path()/".."/name_/startTimeName;
            }
            else
            {
                forcesDir = time_.path()/name_/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(forcesDir);

            // Open new file at start up
            forcesFilePtr_.reset(new OFstream(forcesDir/(type() + ".dat")));

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::immersedBoundaryForceFunctionObject::writeFileHeader()
{
    if (forcesFilePtr_.valid())
    {
        forcesFilePtr_()
            << "# Time" << tab << "Fx" << tab << "Fy" << tab << "Fz" << tab
            << "Cx" << tab << "Cy" << tab << "Cz" << endl;
    }
}


Foam::vector
Foam::immersedBoundaryForceFunctionObject::calcForce() const
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volVectorField& U = 
        mesh.lookupObject<volVectorField>("U");

    const volScalarField& p = 
        mesh.lookupObject<volScalarField>("p");

    const surfaceScalarField& phi = 
        mesh.lookupObject<surfaceScalarField>("phi");

    const volScalarField& gammaExt = 
        mesh.lookupObject<volScalarField>("ibGammaExt");

    const dictionary& transportProperties =
        mesh.lookupObject<dictionary>("transportProperties");

    dimensionedScalar nu(transportProperties.lookup("nu"));

    vector force = fvc::domainIntegrate(gammaExt*fvc::ddt(U)).value();

    {
        volTensorField gradU = fvc::grad(U);

        forAll(mesh.boundary(), patchI)
        {
            if (!mesh.boundary()[patchI].coupled())
            {
                force += gSum
                (
                    phi.boundaryField()[patchI]
                   *U.boundaryField()[patchI]
                );
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            if (!mesh.boundary()[patchI].coupled())
            {
                force += gSum
                (
                    mesh.Sf().boundaryField()[patchI]
                   *p.boundaryField()[patchI]
                );
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            if (!mesh.boundary()[patchI].coupled())
            {
                force -= gSum
                (
                    mesh.Sf().boundaryField()[patchI]
                   &(
                       nu.value()*gradU.boundaryField()[patchI] 
                     + nu.value()*gradU.boundaryField()[patchI].T()
                    )                    
                );
            }
        }
    }

    return force;
}


void Foam::immersedBoundaryForceFunctionObject::write()
{
    makeFile();

    vector force = calcForce();

    vector coeff = force/(0.5*Aref_*sqr(Uref_));

    if (Pstream::master())
    {
        forcesFilePtr_() << time_.value() << tab 
            << force.x() << tab
            << force.y() << tab
            << force.z() << tab
            << coeff.x() << tab
            << coeff.y() << tab
            << coeff.z() << tab << endl;

        Info << "Immersed boundary force coefficients (Cx Cy Cz): "
            << coeff << nl << nl << endl;
    }
}

// ************************************************************************* //
