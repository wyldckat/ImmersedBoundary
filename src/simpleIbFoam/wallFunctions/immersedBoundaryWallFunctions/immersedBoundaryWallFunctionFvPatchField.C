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

#include "immersedBoundaryWallFunctionFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryFvPatchField<Type>(p, iF)
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& tkqrwfpf
)
:
    immersedBoundaryFvPatchField<Type>(tkqrwfpf)
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& tkqrwfpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryFvPatchField<Type>(tkqrwfpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    immersedBoundaryFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::write(Ostream& os) const
{
    immersedBoundaryFvPatchField<Type>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
