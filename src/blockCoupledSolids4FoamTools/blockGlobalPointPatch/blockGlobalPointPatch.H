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
    blockGlobalPointPatch

Description
    Block coupled global patch for solid solver

Author
    Philip Cardiff UCD

SourceFiles
    blockGlobalPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef blockGlobalPointPatch_H
#define blockGlobalPointPatch_H

#include "facePointPatch.H"
#include "blockGlobalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class blockGlobalPointPatch Declaration
\*---------------------------------------------------------------------------*/

class blockGlobalPointPatch
:
    public facePointPatch
{

public:

    //- Runtime type information
    TypeName(blockGlobalPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        blockGlobalPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            facePointPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
