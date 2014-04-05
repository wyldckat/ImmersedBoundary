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

#include "OPstream.H"
#include "IPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void immersedBoundary::makeInvDirichletMatrices() const
{
    if (debug)
    {
        Info<< "immersedBoundary::makeInvDirichletMatrices() : "
            << "making immersed boundary inverse matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invDirichletMatrices_.size() != 0)
    {
        FatalErrorIn("immersedBoundary::makeInvDirichletMatrices()")
            << "immersed boundary inverse least squares matrices already exist"
            << abort(FatalError);
    }

    invDirichletMatrices_.setSize(cells().size());

    const vectorField& C = mesh_.cellCentres();

    scalarField conditionNumber(cells().size(), 0.0);

    const FieldField<Field, vector>& procC = procCentres();

    label nCoeffs = 5;
    if(mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    forAll(invDirichletMatrices_, cellI)
    {
        const labelList& interpCells = cellCells()[cellI];
        const List<labelPair>& interpProcCells = cellProcCells()[cellI];

        vectorField allPoints
        (
            interpCells.size()
          + interpProcCells.size(),
            vector::zero
        );

        if (allPoints.size() < nCoeffs)
        {
            FatalErrorIn("immersedBoundary::makeInvDirichletMatrices()")
                << "allPoints.size() < " << nCoeffs << " : " 
                    << allPoints.size() << abort(FatalError);
        }

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Processor cells
        for (label i=0; i<interpProcCells.size(); i++)
        {
            allPoints[pointID++] = 
                procC
                [
                    interpProcCells[i].first()
                ]
                [
                    interpProcCells[i].second()
                ];
        }

        // Weights 
        scalarField W(allPoints.size(), 1.0);
        vector origin = C[cells()[cellI]];

        scalar maxR = 0;
        for (label i=0; i<allPoints.size(); i++)
        {
            if(mag(allPoints[i] - origin) > maxR)
            {
                maxR = mag(allPoints[i] - origin);
            }
        }
        maxR *= 1.1;
        
        for (label i=0; i<allPoints.size(); i++)
        {
            scalar curR =  mag(allPoints[i] - origin);

            W[i] = 0.5*(1 + ::cos(M_PI*curR/maxR));
        }

        invDirichletMatrices_.set
        (
            cellI, 
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size(),
                0.0
            ) 
        );
        scalarRectangularMatrix& curMatrix = invDirichletMatrices_[cellI];

        scalarRectangularMatrix M
        (
            allPoints.size(),
            nCoeffs,
            0.0
        );

        origin = points()[cellI];
        for(label i=0; i<allPoints.size(); i++)
        {
            scalar X = allPoints[i].x() - origin.x();
            scalar Y = allPoints[i].y() - origin.y();

            label coeff = 0;
            M[i][coeff++] = X;
            M[i][coeff++] = Y;
            M[i][coeff++] = X*Y;
            M[i][coeff++] = sqr(X);
            M[i][coeff++] = sqr(Y);
            
            if(mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[i].z() - origin.z();
                M[i][coeff++] = Z;
                M[i][coeff++] = X*Z;
                M[i][coeff++] = Y*Z;
                M[i][coeff++] = sqr(Z);
            }
        }

        for (label i=0; i<M.n(); i++)
        {
            for (label j=0; j<M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        scalarSquareMatrix lsM(nCoeffs, 0.0);

        for (label i=0; i<lsM.n(); i++)
        {
            for (label j=0; j<lsM.m(); j++)
            {
                for (label k=0; k<M.n(); k++)
                {
                    lsM[i][j] += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate matrix norm
        scalar maxRowSum = 0.0;
        for (label i=0; i<lsM.n(); i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<lsM.m(); j++)
            {
                curRowSum += lsM[i][j];
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }
        conditionNumber[cellI] = maxRowSum;

        // Calculate inverse
        scalarSquareMatrix invLsM = lsM.LUinvert();

        for (label i=0; i<lsM.n(); i++)
        {
            for (label j=0; j<M.n(); j++)
            {
                for (label k=0; k<lsM.n(); k++)
                {
                    curMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
                }
            }
        }

        // Calculate condition number
        maxRowSum = 0.0;
        for (label i=0; i<lsM.n(); i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<lsM.m(); j++)
            {
                curRowSum += invLsM[i][j];
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }
        conditionNumber[cellI] *= maxRowSum;
    }

    Info << "Max dirichlet matrix condition number: " 
        << gMax(conditionNumber) << endl;
}


void immersedBoundary::makeInvNeumannMatrices() const
{
    if (debug)
    {
        Info<< "immersedBoundary::makeInvNeumannMatrices() : "
            << "making immersed boundary inverse least sqares matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invNeumannMatrices_.size() != 0)
    {
        FatalErrorIn("immersedBoundary::makeInvNeumannMatrices()")
            << "immersed boundary inverse least squares matrices already exist"
            << abort(FatalError);
    }

    invNeumannMatrices_.setSize(cells().size());

    const vectorField& C = mesh_.cellCentres();

    scalarField conditionNumber(cells().size(), 0.0);

    const FieldField<Field, vector>& procC = procCentres();

    label nCoeffs = 6;
    if(mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    forAll(invNeumannMatrices_, cellI)
    {
        const labelList& interpCells = cellCells()[cellI];
        const List<labelPair>& interpProcCells = cellProcCells()[cellI];

        vectorField allPoints
        (
            interpCells.size() + 1
          + interpProcCells.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // IB point
        allPoints[pointID++] = points()[cellI];

        // Processor cells
        for (label i=0; i<interpProcCells.size(); i++)
        {
            allPoints[pointID++] = 
                procC
                [
                    interpProcCells[i].first()
                ]
                [
                    interpProcCells[i].second()
                ];
        }

        // Weights 
        scalarField W(allPoints.size(), 1.0);

        vector origin = C[cells()[cellI]];

        scalar maxR = 0;
        for (label i=0; i<allPoints.size(); i++)
        {
            if(mag(allPoints[i] - origin) > maxR)
            {
                maxR = mag(allPoints[i] - origin);
            }
        }

        for (label i=0; i<allPoints.size(); i++)
        {
            scalar curR =  mag(allPoints[i] - origin);

            W[i] = 0.5*(1 + ::cos(M_PI*curR/maxR));
        }

        invNeumannMatrices_.set
        (
            cellI, 
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size(),  
                0.0
            ) 
        );
        scalarRectangularMatrix& curMatrix = invNeumannMatrices_[cellI];

        scalarRectangularMatrix M
        (
            allPoints.size(),
            nCoeffs, 
            0.0
        );

        pointID = 0;
        origin = points()[cellI];
        for(label i=0; i<interpCells.size(); i++)
        {
            scalar X = allPoints[pointID].x() - origin.x();
            scalar Y = allPoints[pointID].y() - origin.y();

            M[pointID][0] = 1.0;
            M[pointID][1] = X;
            M[pointID][2] = Y;
            M[pointID][3] = X*Y;
            M[pointID][4] = sqr(X);
            M[pointID][5] = sqr(Y);

            if(mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[pointID].z() - origin.z();

                M[pointID][6] = Z;
                M[pointID][7] = X*Z;
                M[pointID][8] = Y*Z;
                M[pointID][9] = sqr(Z);
            }

            pointID++;
        }

        scalar X = allPoints[pointID].x() - origin.x();
        scalar Y = allPoints[pointID].y() - origin.y();

        M[pointID][0] = 0;
        M[pointID][1] = normals()[cellI].x();
        M[pointID][2] = normals()[cellI].y();
        M[pointID][3] = 
        (
            normals()[cellI].x()*Y
          + normals()[cellI].y()*X
        );
        M[pointID][4] = 
            2*normals()[cellI].x()*X;
        M[pointID][5] = 
            2*normals()[cellI].y()*Y;

        if(mesh_.nGeometricD() == 3)
        {
            scalar Z = allPoints[pointID].z() - origin.z();

            M[pointID][6] = normals()[cellI].z();
            M[pointID][7] =
            (
                normals()[cellI].x()*Z
              + normals()[cellI].z()*X
            );
            M[pointID][8] =
            (
                normals()[cellI].y()*Z
              + normals()[cellI].z()*Y
            );
            M[pointID][9] = 2*normals()[cellI].z()*Z;
        }

        pointID++;

        for(label i=0; i<interpProcCells.size(); i++)
        {
            scalar X = allPoints[pointID].x() - origin.x();
            scalar Y = allPoints[pointID].y() - origin.y();

            M[pointID][0] = 1.0;
            M[pointID][1] = X;
            M[pointID][2] = Y;
            M[pointID][3] = X*Y;
            M[pointID][4] = sqr(X);
            M[pointID][5] = sqr(Y);

            if(mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[pointID].z() - origin.z();

                M[pointID][6] = Z;
                M[pointID][7] = X*Z;
                M[pointID][8] = Y*Z;
                M[pointID][9] = sqr(Z);
            }

            pointID++;
        }

        for (label i=0; i<M.n(); i++)
        {
            for (label j=0; j<M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        scalarSquareMatrix lsM(nCoeffs, 0.0);

        for (label i=0; i<lsM.n(); i++)
        {
            for (label j=0; j<lsM.m(); j++)
            {
                for (label k=0; k<M.n(); k++)
                {
                    lsM[i][j] += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate matrix norm
        scalar maxRowSum = 0.0;
        for (label i=0; i<lsM.n(); i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<lsM.m(); j++)
            {
                curRowSum += lsM[i][j];
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }
        conditionNumber[cellI] = maxRowSum;
        
        // Calculate inverse
        scalarSquareMatrix invLsM = lsM.LUinvert();

        for (label i=0; i<lsM.n(); i++)
        {
            for (label j=0; j<M.n(); j++)
            {
                for (label k=0; k<lsM.n(); k++)
                {
                    curMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
                }
            }
        }

        // Calculate condition number
        maxRowSum = 0.0;
        for (label i=0; i<lsM.n(); i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<lsM.m(); j++)
            {
                curRowSum += invLsM[i][j];
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }
        conditionNumber[cellI] *= maxRowSum;
    }

    Info << "Max neumann matrix condition number: " 
        << gMax(conditionNumber) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * *  //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
