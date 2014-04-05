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
    immersedBoundaryForceFunctionObject

Description
    FunctionObject setting mesh motion 

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    meshCheckFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryForceFunctionObject_H
#define immersedBoundaryForceFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "Switch.H"
#include "OFstream.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class immersedBoundaryForceFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryForceFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Force file ptr
        autoPtr<OFstream> forcesFilePtr_;

        //- Reference velocity for force coefficient calculation
        scalar Uref_;

        //- Reference area for force coefficient calculation
        scalar Aref_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        immersedBoundaryForceFunctionObject
        (
            const immersedBoundaryForceFunctionObject&
        );

        //- Disallow default bitwise assignment
        void operator=(const immersedBoundaryForceFunctionObject&);

public:

    //- Runtime type information
    TypeName("immersedBoundaryForce");


    // Constructors

        //- Construct from components
        immersedBoundaryForceFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- If the forces file has not been created create it
        void makeFile();

        //- Output file header information
        void writeFileHeader();

        //- Output forces
        void write();
 
        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Calculate and return forces 
        vector calcForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
