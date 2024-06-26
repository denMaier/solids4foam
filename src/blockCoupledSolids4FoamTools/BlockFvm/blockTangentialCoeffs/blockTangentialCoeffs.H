/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Description
    To avoid code repetition, the implementation of inserting the coefficients
    for the block laplacian, laplacianTranspose and laplacianTrace operators are
    given here.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef logVolSymmTensorField_H
#define logVolSymmTensorField_H

#include "solidPolyMesh.H"
#include "newLeastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{
    // Natural log of a symmTensor field
    //tmp<volSymmTensorField> log(const volSymmTensorField& vf);

    //- Insert coefficients for the block discretisation of the tangential
    //  derivative terms
    void blockFvmInsertCoeffsTang
    (
        const solidPolyMesh& solidMesh,
        const surfaceScalarField& muf,
        const surfaceScalarField& lambdaf,
        const GeometricField<vector, fvPatchField, volMesh>& vf,
        Field<vector>& blockB,
        BlockLduMatrix<vector>& blockM,
        const newLeastSquaresVolPointInterpolation& volToPointInterp,
        const int op
    );

    //- Insert the coefficients for the current face
    void blockFvmInsertCoeffsTangForFace
    (
        tensorField& d,
        tensorField& u,
        tensorField& l,
        vector& blockB,
        Map<tensorField>& pointProcFacesCoeffs,
        Map<tensorField>& pointProcBndFacesCoeffs,
        Map<tensorField>& pointProcCellsCoeffs,
        Map<tensorField>& gPtNgbProcBndFaceCoeffs,
        Map<tensorField>& gPtNgbProcCellCoeffs,
        const label cellI,
        const face& curFace,
        const bool flipFace,
        const bool nonOrthogonalMesh,
        const scalar faceMu,
        const scalar faceLambda,
        const vector& faceN,
        const vector& faceK,
        const vectorField& points,
        const vectorField& origins,
        const scalarField& refL,
        const PtrList<scalarRectangularMatrix>& invMatrices,
        const FieldField<Field, scalar>& weights,
        const bool enforcePointConstraints,
        const Map<vector>& pointFixedComp,
        const Map<symmTensor>& pointFixedDir,
        const labelListList& pointCells,
        const List< Tuple2<vector, tensor> >& mirrorPlaneTrans,
        const solidPolyMesh& solidMesh,
        const labelListList& pointBndFaces,
        const labelListList& pointCyclicFaces,
        const labelListList& pointProcFaces,
        const Map<Field<vector> >& gPtNgbProcBndFaceFieldData,
        const Map<Field<vector> >& gPtNgbProcCellFieldData,
        const Map< List<labelPair> >& pointProcCells,
        const List< List<labelPair> >& pointProcBndFaces,
        const int op
    );

    //- Calculate the coefficient
    //  If fixed point is enforced then there can be a blockB contribution
    void blockFvmCalculateCoeff
    (
        tensor& coeff,
        vector& blockB,
        const scalar wCell,
        const vector& m,
        const vector& dr,
        const scalar faceMu,
        const scalar faceLambda,
        const vector& faceN,
        const tensor& LeFaceN,
        const scalar kDotLe,
        const bool nonOrthogonalMesh,
        const Tuple2<vector, tensor>& mirrorPlaneTrans,
        const bool pointHasFixedComp,
        const vector& pointFixedComp,
        const symmTensor& pointFixedDir,
        const int op
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
