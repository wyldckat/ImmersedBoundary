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

Class
    immersedBoundary

Description

SourceFiles
    immersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef ImmersedBoundary_H
#define ImmersedBoundary_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "octree.H"
#include "octreeDataTriSurface.H"
#include "triSurfaceSearch.H"
#include "triSurfaceFields.H"

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundary
:
    public triSurfaceMesh
{
    // Private data

        //- Finite volume mesh reference
        const fvMesh& mesh_;

        //- Internal or external flow calculation
        Switch internalFlow_;
    
    // Demand-driven data

        //- triSurface search algorithm
        mutable triSurfaceSearch* triSurfSearchPtr_;

        //- Fluid cells indicator
        mutable volScalarField* gammaPtr_;

        //- Fluid cells indicator
        mutable volScalarField* gammaExtPtr_;

        //- Fluid faces indicator
        mutable surfaceScalarField* sGammaPtr_;

        //- List of fluid cells next to immersed boundary (IB cells)
        mutable labelList* cellsPtr_;

        //- List of faces for which one neighbour is an IB cell
        //  and another neighbour is a live fluid cell (IB faces)
        mutable labelList* facesPtr_;

        //- List of fluid faces for which one neighbour is an IB cell
        //  and another neighbour is a dead cell (inside IB faces)
        mutable labelList* insideFacesPtr_;

        //- List of internal faces in the region bounded by IB faces
        //  inside IB faces
        mutable labelList* internalFacesPtr_;

        //- Points at the immersed boundary (IB points) 
        //  nearest to the IB cell centres
        mutable vectorField* pointsPtr_;

        //- Normals at IB points 
        mutable vectorField* normalsPtr_;

        //- List of faces (triangles) which are part of IB mesh
        //  nearest to the IB cell centres
        mutable labelList* hitFacesPtr_;

        //- Neighbour cells for immersed boundary cells
        //  (extended stencil) 
        mutable labelListList* cellCellsPtr_;

        //- List of cells needed by neighbour processors
        mutable labelListList* procCellsPtr_;

        //- Centres of cells from neighbour processors
        mutable FieldField<Field, vector>* procCentresPtr_;

        //- Cell-proc-cell addressing
        mutable List<List<labelPair> >* cellProcCellsPtr_;

        //- Dead cells list
        mutable labelList* deadCellsPtr_;

        //- Extended dead cells list (dead cells + IB cells)
        mutable labelList* deadCellsExtPtr_;

        //- List of live cells
        mutable labelList* liveCellsPtr_;

        //- Average IB cell sizes
        mutable scalarField* cellSizesPtr_;

        //- Inverse interpolation matrices for Dirichlet BC at the IB
        mutable PtrList<scalarRectangularMatrix> invDirichletMatrices_;

        //- Inverse interpolation matrices for Neumann BC at the IB
        mutable PtrList<scalarRectangularMatrix> invNeumannMatrices_;

        //- Cells for refinement
        mutable labelList* refinementCellsPtr_;

        //- Wall distance
        mutable volScalarField* wallDistancePtr_;

        //- IB cell centre distances to IB
        mutable scalarField* deltaPtr_;

    // Private Member Functions

        // Make demand-driven data

            //- Make triSurface search algorithm
            void makeTriSurfSearch() const;

            //- Make fluid cells indicator
            void makeGamma() const;

            //- Make extended fluid cells indicator
            void makeGammaExt() const;

            //- Make fluid faces indicator
            void makeSGamma() const;

            //- Make list of cells next to immersed boundary
            void makeCells() const;

            //- Add corner points to IB cells list
            void addCornerCells() const;

            //- Make IB faces
            void makeFaces() const;

            //- Make inside IB faces
            void makeInsideFaces() const;

            //- Make internal IB faces
            void makeInternalFaces() const;

            //- Make immersed boundary points and normals
            void makePointsAndNormals() const;

            //- Make extended IB cells stencils
            void makeCellCells() const;

            //- Make list of dead cells
            void makeDeadCells() const;

            //- Make extended list of dead cells
            void makeDeadCellsExt() const;

            //- Make list of live cells
            void makeLiveCells() const;

            //- Make immersed boundary cell sizes
            void makeCellSizes() const;

            //- Make inverse Dirichlet interpolation matrices
            void makeInvDirichletMatrices() const;

            //- Make inverse Neumann interpolation matrices
            void makeInvNeumannMatrices() const;

            //- Make list of cells for refinement
            void makeRefinementCells() const;

            //- Make wall and IB distance
            void makeWallDistance() const;

            //- Make distance between IB cell centres 
            //  and corresponding IB points
            void makeDelta() const;

            //- Find nearest cell
            label findNearestCell(const point& location) const;

            //- Return extended cell-cell addressing
            void findCellCells
            (
                const vector& pt,
                const label cellID, 
                labelList& cellCells
            ) const;

            //- Set value at IB cells
            template<class Type>
            void setValues
            (
                fvMatrix<Type>& eqn, 
                const labelList& cellLabels,
                const Field<Type>& values
            ) const;

            //- From refineMesh application
            label axis(const vector& normal) const;

            //- From refineMesh application
            label twoDNess(const polyMesh& mesh) const;

            //- Calc cell size
            scalar cellSize(label cellID) const;

            //- Calc cell projection area
            scalar cellProjection(label cellID, const vector& dir) const;

        // Storage management

            //- Clear all demand-driven data
            void clearOut();


        //- Disallow default bitwise copy construct
        immersedBoundary(const immersedBoundary&);

        //- Disallow default bitwise assignment
        void operator=(const immersedBoundary&);

public:

    // Declare name of the class and it's debug switch
    ClassName("immersedBoundary");


    // Constructors
    immersedBoundary(const fvMesh& mesh);


    // Destructor
    ~immersedBoundary();


    // Member Functions

        //- Move points
        virtual void movePoints(const pointField& newPoints);

        //- Impose Dirichlet BC at IB cells
        template<class Type>
        tmp<Field<Type> > imposeDirichletBoundaryCondition
        (
            fvMatrix<Type>& eqn, 
            const Field<Type>& ibValues
        ) const;

        //- Impose Neuman BC at IB cells
        template<class Type>
        tmp<Field<Type> > imposeNeumannBoundaryCondition
        (
            GeometricField<Type, fvPatchField, volMesh>& psi,
            const Field<Type>& ibSnGradPsi
        ) const;

        //- Correct eqn by setting zero field for dead cells
        template<class Type>
        void correctEqn
        (
            fvMatrix<Type>& psiEqn
        ) const;

        //- Calc flux for IB faces
        tmp<scalarField> flux
        (
            const volVectorField& U
        ) const;

        //- Calc flux for IB faces
        tmp<scalarField> flux
        (
            const volVectorField& U,
            const vectorField& ibVelocity
        ) const;

        //- Set flux for IB faces
        void setFlux
        (
            surfaceScalarField& phi,
            const scalarField& ibPhi
        ) const;

        //- Scale and set flux for IB faces
        void scaleFlux
        (
            surfaceScalarField& phi, 
            const surfaceScalarField& phiU
        ) const;

        //- Calc pressure gradient using extrapolation near IB
        tmp<volVectorField> grad(const volScalarField& p) const;

        // Access

            //- Get triSurface search object
            const triSurfaceSearch& triSurfSearch() const;

            //- Get flud cells indicator
            const volScalarField& gamma() const;

            //- Get extended flud cells indicator
            const volScalarField& gammaExt() const;

            //- Get fluid faces indicator
            const surfaceScalarField& sGamma() const;

            //- Get list of IB cells
            const labelList& cells() const;

            //- Get list of IB faces
            const labelList& faces() const;

            //- Get list of inside IB faces
            const labelList& insideFaces() const;

            //- Get list of internal IB faces
            const labelList& internalFaces() const;

            //- Get IB points
            const vectorField& points() const;

            //- Get IB normals
            const vectorField& normals() const;

            //- Get list of triangles in IB mesh nearest
            //  nearest to IB cell centres
            const labelList& hitFaces() const;

            //- Get IB cell extended stencil
            const labelListList& cellCells() const;

            //- Get neighbour proc centres
            const FieldField<Field, vector>& procCentres() const;

            //- Get nighbour centres addressing
            const List<List<labelPair> >& cellProcCells() const;

            //- Get neighbour proc cells
            const labelListList& procCells() const;

            //- Get dead cells
            const labelList& deadCells() const;

            //- Get extended dead cells
            const labelList& deadCellsExt() const;

            //- Get live cells
            const labelList& liveCells() const;

            //- Get immersed boundary cell sizes
            const scalarField& cellSizes() const;

            //- Get inverse Dirichlet interpolation matrix
            const PtrList<scalarRectangularMatrix>& 
            invDirichletMatrices() const;

            //- Get inverse Neumann interpolation matrix
            const PtrList<scalarRectangularMatrix>& 
            invNeumannMatrices() const;

            //- Get list of cells for refinement
            const labelList& refinementCells() const;

            //- Get wall (IB) distance
            const volScalarField& wallDistance() const;

            //- Get distance to IB
            const scalarField& delta() const;

            //- Refine mesh near immersed boundary
            void refineMesh(polyMesh& pMesh, const labelList& refCells) const;

        // Write

        bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
