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

template<class Type>
tmp<Field<Type> >  immersedBoundary::imposeDirichletBoundaryCondition
(
    fvMatrix<Type>& eqn, 
    const Field<Type>& ibValues
) const
{
    Field<Type>& psiI = eqn.psi().internalField();

    tmp<Field<Type> > tIbSnGradPsi
    (
        new Field<Type>(cells().size(), pTraits<Type>::zero)
    );
    Field<Type>& ibSnGradPsi = tIbSnGradPsi();

    Field<Type> polyPsi(psiI, cells());

    const vectorField& C = mesh_.cellCentres();

    label nCoeffs = 5;
    if(mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    label counter = 0;

    scalarField error(cells().size(), 0);

    const PtrList<scalarRectangularMatrix>& invMat = invDirichletMatrices();

    do
    {
        counter++;

#       include "sendAndReceive.H"

        forAll(cells(), cellI)
        {
            label curCell = cells()[cellI];

            const labelList& curCells = cellCells()[cellI];

            const List<labelPair>& curProcCells = cellProcCells()[cellI];

            const scalarRectangularMatrix& curInvMatrix = invMat[cellI];

            Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
            Field<Type> source
            (
                curCells.size() + curProcCells.size(), 
                pTraits<Type>::zero
            );

            label pointID = 0;
            for (label i=0; i<curCells.size(); i++)
            {
                source[pointID++] = psiI[curCells[i]] - ibValues[cellI];
            }

            for (label i=0; i<curProcCells.size(); i++)
            {
                source[pointID++] = 
                    procPsi
                    [
                        curProcCells[i].first()
                    ]
                    [
                        curProcCells[i].second()
                    ] 
                  - ibValues[cellI];
            }

            for (label i=0; i<nCoeffs; i++)
            {
                for (label j=0; j<source.size(); j++)
                {
                    coeffs[i] += curInvMatrix[i][j]*source[j];
                }
            }

            Type oldPolyPsi = polyPsi[cellI];

            vector R =  C[curCell] - points()[cellI];

            polyPsi[cellI] = 
                ibValues[cellI]
              + coeffs[0]*R.x()
              + coeffs[1]*R.y()
              + coeffs[2]*R.x()*R.y()
              + coeffs[3]*sqr(R.x())
              + coeffs[4]*sqr(R.y());

            if (mesh_.nGeometricD() == 3)
            {
                polyPsi[cellI] +=
                    coeffs[5]*R.z()
                  + coeffs[6]*R.x()*R.z()
                  + coeffs[7]*R.y()*R.z()
                  + coeffs[8]*sqr(R.z());
            }

            ibSnGradPsi[cellI] =
                coeffs[0]*normals()[cellI].x()
              + coeffs[1]*normals()[cellI].y();

            if (mesh_.nGeometricD() == 3)
            {
                ibSnGradPsi[cellI] +=
                    coeffs[5]*normals()[cellI].z();
            }

            error[cellI] = mag(polyPsi[cellI] - oldPolyPsi)
               /(mag(oldPolyPsi) + SMALL);
        }

        forAll(polyPsi, cellI)
        {
            psiI[cells()[cellI]] = polyPsi[cellI];
        }
    }
    while(gMax(error) > 1e-6 && counter < 200);

    if (counter==200)
    {
        Info << eqn.psi().name() << ", error, max: " << gMax(error) 
            << ", min: " << gMin(error) 
            << ", avg: "  << gAverage(error) << endl;
    }

    // Dead cells
    Field<Type> deadCellsPsi(deadCells().size(), pTraits<Type>::zero);
    eqn.setValues(deadCells(), deadCellsPsi);

    // IB cells
    setValues(eqn, cells(), polyPsi);

    return tIbSnGradPsi;
}


template<class Type>
tmp<Field<Type> > immersedBoundary::imposeNeumannBoundaryCondition
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    const Field<Type>& ibSnGradPsi
) const
{
    Field<Type>& psiI = psi.internalField();

    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(cells().size(), pTraits<Type>::zero)
    );
    scalarField& ibPsi = tIbPsi();

    const vectorField& C = mesh_.cellCentres();

    label nCoeffs = 6;
    if(mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMat = invNeumannMatrices();

    label counter = 0;
    scalarField error(cells().size(), 0);

    do
    {
        counter++;

#       include "sendAndReceive.H"

        forAll(cells(), cellI)
        {
            label curCell = cells()[cellI];

            const labelList& curCells = cellCells()[cellI];
            const List<labelPair>& curProcCells = cellProcCells()[cellI];

            const scalarRectangularMatrix& curInvMatrix = invMat[cellI];

            Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
            Field<Type> source
            (
                curCells.size() + 1 + curProcCells.size(), 
                pTraits<Type>::zero
            );

            label pointID = 0;
            for (label i=0; i<curCells.size(); i++)
            {
                source[pointID++] = psiI[curCells[i]];
            }

            source[pointID++] = ibSnGradPsi[cellI];

            for (label i=0; i<curProcCells.size(); i++)
            {
                source[pointID++] = 
                    procPsi
                    [
                        curProcCells[i].first()
                    ]
                    [
                        curProcCells[i].second()
                    ];
            }

            for (label i=0; i<nCoeffs; i++)
            {
                for (label j=0; j<source.size(); j++)
                {
                    coeffs[i] += curInvMatrix[i][j]*source[j];
                }
            }

            Type oldPsi = psiI[curCell];

            vector ibR =  C[curCell] - points()[cellI];

            psiI[curCell] =
                coeffs[0]
              + coeffs[1]*ibR.x()
              + coeffs[2]*ibR.y()
              + coeffs[3]*ibR.x()*ibR.y()
              + coeffs[4]*sqr(ibR.x())
              + coeffs[5]*sqr(ibR.y());

            if (mesh_.nGeometricD() == 3)
            {
                psiI[curCell] +=
                    coeffs[6]*ibR.z()
                  + coeffs[7]*ibR.x()*ibR.z()
                  + coeffs[8]*ibR.y()*ibR.z()
                  + coeffs[9]*sqr(ibR.z());
            }

            ibPsi[cellI] = coeffs[0];

            error[cellI] = mag(psiI[curCell] - oldPsi)
               /(mag(oldPsi) + SMALL);
        }
    }
    while(gMax(error) > 1e-6 && counter < 200);

    if (counter==200)
    {
        Info << psi.name() << ", error, max: " << gMax(error) 
            << ", min: " << gMin(error) 
            << ", avg: "  << gAverage(error) << endl;
    }

    return tIbPsi;
}


template<class Type>
void immersedBoundary::correctEqn
(
    fvMatrix<Type>& psiEqn
) const
{
    scalarField& Diag = psiEqn.diag();

    forAll(deadCellsExt(), cellI)
    {
        if (mag(Diag[deadCellsExt()[cellI]]) < SMALL)
        {
            Diag[deadCellsExt()[cellI]] = 1.0;
        }
    }

    // Dead cells
    Field<Type> deadCellsPsi(deadCellsExt().size(), pTraits<Type>::zero);
    psiEqn.setValues(deadCellsExt(), deadCellsPsi);
}


template<class Type>
void immersedBoundary::setValues
(
    fvMatrix<Type>& eqn,
    const labelList& cellLabels,
    const Field<Type>& values
) const
{
    const fvMesh& mesh = eqn.psi().mesh();

    const cellList& cells = mesh.cells();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    scalarField& Diag = eqn.diag();
    Field<Type>& psi_ = eqn.psi();
    Field<Type>& source_ = eqn.source();

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        psi_[celli] = values[i];
        source_[celli] = values[i]*Diag[celli];

        if (eqn.symmetric() || eqn.asymmetric())
        {
            const cell& c = cells[celli];

            forAll(c, j)
            {
                label facei = c[j];

                if (mesh.isInternalFace(facei))
                {
                    if (eqn.symmetric())
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= 
                                eqn.upper()[facei]*values[i];
                        }
                        else
                        {
                            source_[own[facei]] -= 
                                eqn.upper()[facei]*values[i];
                        }

                        eqn.upper()[facei] = 0.0;
                    }
                    else
                    {
                        if (celli == own[facei])
                        {
//                             source_[nei[facei]] -= lower()[facei]*values[i];
                            eqn.upper()[facei] = 0.0;
                        }
                        else
                        {
//                             source_[own[facei]] -= upper()[facei]*values[i];
                            eqn.lower()[facei] = 0.0;
                        }

//                         upper()[facei] = 0.0;
//                         lower()[facei] = 0.0;
                    }
                }
                else
                {
                    label patchi = mesh.boundaryMesh().whichPatch(facei);

                    if (eqn.internalCoeffs()[patchi].size())
                    {
                        label patchFacei =
                            mesh.boundaryMesh()[patchi].whichFace(facei);

                        eqn.internalCoeffs()[patchi][patchFacei] =
                            pTraits<Type>::zero;

                        eqn.boundaryCoeffs()[patchi][patchFacei] =
                            pTraits<Type>::zero;
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
