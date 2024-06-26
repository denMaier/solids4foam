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

Class
    analyticalPlateHoleTractionPointPatchVectorField

Description
    The point displacement version of analyticalPlateHoleTractionFvPatchField.

    Used to enforce the analytical traction on the hole in a plate case.
SourceFiles
    analyticalPlateHoleTractionPointPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef analyticalPlateHoleTractionPointPatchVectorField_H
#define analyticalPlateHoleTractionPointPatchVectorField_H

#include "solidTractionPointPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
#ifdef OPENFOAM_NOT_EXTEND
    typedef pointPatchFieldMapper PointPatchFieldMapper;
#endif

/*---------------------------------------------------------------------------*\
            Class analyticalPlateHoleTractionPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class analyticalPlateHoleTractionPointPatchVectorField
:
    public solidTractionPointPatchVectorField
{

    // Private data

        //- Far field traction
        //  It is assumed that it acts in the positive X direction
        const scalar T_;

        //- Hole radius
        const scalar holeR_;

    // Private Member Functions

        //- Analytical solution
        symmTensor plateHoleSolution(const vector& C);


public:

    //- Runtime type information
    TypeName("pointAnalyticalPlateHoleTraction");


    // Constructors

        //- Construct from patch and internal field
        analyticalPlateHoleTractionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticalPlateHoleTractionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        analyticalPlateHoleTractionPointPatchVectorField
        (
            const analyticalPlateHoleTractionPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

#ifndef OPENFOAM_ORG
        //- Construct as copy
        analyticalPlateHoleTractionPointPatchVectorField
        (
            const analyticalPlateHoleTractionPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new analyticalPlateHoleTractionPointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        analyticalPlateHoleTractionPointPatchVectorField
        (
            const analyticalPlateHoleTractionPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new analyticalPlateHoleTractionPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto
            //  this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        //- Insert boundary value into the internal field
        //  Fixed values need to be inserted first. In the second go,
        //  there's nothing to do
        virtual void initEvaluate
        (
#ifdef OPENFOAM_NOT_EXTEND
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
