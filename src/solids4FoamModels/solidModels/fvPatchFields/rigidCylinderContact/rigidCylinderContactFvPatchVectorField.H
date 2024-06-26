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
    rigidCylinderContactFvPatchVectorField

Description
    A rigid cylinder in contact with the patch.

Author
    Philip Cardiff, UCD.

SourceFiles
    rigidCylinderContactFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rigidCylinderContactFvPatchVectorField_H
#define rigidCylinderContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class rigidCylinderContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class rigidCylinderContactFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{

    // Private Data

        //- Cylinder radius
        const scalar radius_;

        //- Time-varying cylinder centre
        const interpolationTable<vector> cylinderCentre_;

        //- Stiffness of the penalty springs
        const scalar penaltyStiffness_;

        //- Under-relaxation factor to imrpove convergence
        const scalar relaxFactor_;

public:

    //- Runtime type information
    TypeName("rigidCylinderContact");


    // Constructors

        //- Construct from patch and internal field
        rigidCylinderContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rigidCylinderContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  rigidCylinderContactFvPatchVectorField onto a new patch
        rigidCylinderContactFvPatchVectorField
        (
            const rigidCylinderContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rigidCylinderContactFvPatchVectorField
        (
            const rigidCylinderContactFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new rigidCylinderContactFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        rigidCylinderContactFvPatchVectorField
        (
            const rigidCylinderContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new rigidCylinderContactFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
