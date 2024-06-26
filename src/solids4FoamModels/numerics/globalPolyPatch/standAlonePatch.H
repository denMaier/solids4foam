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

Typedef
    Foam::standAlonePatch

Description
    standAlonePatch is a primiteve patch that holds its point and face list.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef standAlonePatch_H
#define standAlonePatch_H

#ifdef OPENFOAM_NOT_EXTEND
    #include "PrimitivePatch.H"
#else
    #include "PrimitivePatchTemplate.H"
#endif
#include "face.H"
#include "List.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
#ifdef OPENFOAM_ORG
    typedef PrimitivePatch<List<face>, pointField> standAlonePatch;
#elif OPENFOAM_COM
    typedef PrimitivePatch<faceList, pointField> standAlonePatch;
#else
    typedef PrimitivePatch<face, List, pointField> standAlonePatch;
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
