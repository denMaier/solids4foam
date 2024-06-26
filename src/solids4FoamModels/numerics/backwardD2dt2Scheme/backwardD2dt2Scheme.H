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
    Foam::fv::backwardD2dt2Scheme

Description
    Second-order backward implicit d2dt2 using the current and four previous
    time-step values.

SourceFiles
    backwardD2dt2Scheme.C

Author
    Zeljko Tukovic
    Some additions by Philip Cardiff

\*---------------------------------------------------------------------------*/

#ifndef backwardD2dt2Scheme_H
#define backwardD2dt2Scheme_H

#include "d2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class backwardD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class backwardD2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        scalar deltaT0_
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Disallow default bitwise copy construct
        backwardD2dt2Scheme(const backwardD2dt2Scheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const backwardD2dt2Scheme&) = delete;


public:

    //- Runtime type information
    TypeName("backward");


    // Constructors

        //- Construct from mesh
        backwardD2dt2Scheme(const fvMesh& mesh)
        :
            d2dt2Scheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        backwardD2dt2Scheme(const fvMesh& mesh, Istream& is)
        :
            d2dt2Scheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
