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
    cantileverTractionPointPatchVectorField

Description
    The boundary condition applies the analytical point traction to the
    fixed end of a simple bending cantilever test case, as described in:

    C.E. Augarde, A.J. Deeks, The use of Timoshenko’s exact solution for a
    cantilever beam in adaptive analysis. Finite Elements in Analysis and
    Design, 44, 2008, 595–601, 10.1016/j.finel.2008.01.010.

    It is assumed that the origin is at the centre of the fixed patch.

    The beam is L long and D deep (thick).

SourceFiles
    cantileverTractionPointPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef cantileverTractionPointPatchVectorField_H
#define cantileverTractionPointPatchVectorField_H

#include "solidTractionPointPatchVectorField.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "pointPatchFieldMapper.H"
#else
    #include "PointPatchFieldMapper.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
#ifdef OPENFOAM_NOT_EXTEND
    typedef pointPatchFieldMapper PointPatchFieldMapper;
#endif

/*---------------------------------------------------------------------------*\
            Class cantileverTractionPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class cantileverTractionPointPatchVectorField
:
    public solidTractionPointPatchVectorField
{

    // Private data

        //- Load applied in the minus y direction at the other end of the beam
        const scalar P_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Length of the beam
        const scalar L_;

        //- Depth of the beam
        const scalar D_;

        //- Second moment of area
        const scalar I_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("pointCantileverTraction");


    // Constructors

        //- Construct from patch and internal field
        cantileverTractionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cantileverTractionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        cantileverTractionPointPatchVectorField
        (
            const cantileverTractionPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

#ifndef OPENFOAM_ORG
        //- Construct as copy
        cantileverTractionPointPatchVectorField
        (
            const cantileverTractionPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new cantileverTractionPointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        cantileverTractionPointPatchVectorField
        (
            const cantileverTractionPointPatchVectorField&,
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
                new cantileverTractionPointPatchVectorField(*this, iF)
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
