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

\*---------------------------------------------------------------------------*/

#include "globalPolyPatch.H"
#include "FieldSumOp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::patchPointToGlobal
(
    const Field<Type>& pField
) const
{
    if (pField.size() != patch().nPoints())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > globalPolyPatch::patchPointToGlobal\n"
            "(\n"
            "    const Field<Type>& pField\n"
            ") const"
        )   << "Patch field does not correspond to patch points.  Patch size: "
            << patch().nPoints() << " field size: " << pField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tgField
    (
        new Field<Type>(globalPatch().nPoints(), pTraits<Type>::zero)
    );
#ifdef OPENFOAM_NOT_EXTEND
    Field<Type>& gField = tgField.ref();
#else
    Field<Type>& gField = tgField();
#endif

    if (Pstream::parRun())
    {
        // PC, 16/12/17
        // We have removed duplicate points so multiple local processor points
        // may map to the same global point, which we will account for using
        // the nPoints field
        scalarField nPoints(gField.size(), 0.0);

        const labelList& addr = pointToGlobalAddr();

        forAll(addr, i)
        {
            const label globalPointID = addr[i];
            gField[globalPointID] = pField[i];
            nPoints[globalPointID] += 1.0;
        }

        // Global comm
        reduce(gField, FieldSumOp<Type>());
        reduce(nPoints, FieldSumOp<scalar>());
        gField /= nPoints;
    }
    else
    {
        gField = pField;
    }

    return tgField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::globalPointToPatch
(
    const Field<Type>& gField
) const
{
    if (gField.size() != globalPatch().nPoints())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > globalPolyPatch::globalPointToPatch\n"
            "(\n"
            "    const Field<Type>& gField\n"
            ") const"
        )   << "Patch field does not correspond to global patch points.  "
            << "Global patch size: " << globalPatch().nPoints()
            << " field size: " << gField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tpField
    (
        new Field<Type>(patch().nPoints(), pTraits<Type>::zero)
    );
#ifdef OPENFOAM_NOT_EXTEND
    Field<Type>& pField = tpField.ref();
#else
    Field<Type>& pField = tpField();
#endif

    if (Pstream::parRun())
    {
        const labelList& addr = pointToGlobalAddr();

        forAll (addr, i)
        {
            pField[i] = gField[addr[i]];
        }
    }
    else
    {
        pField = gField;
    }

    return tpField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::patchFaceToGlobal
(
    const Field<Type>& pField
) const
{
    if (pField.size() != patch().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > globalPolyPatch::patchFaceToGlobal\n"
            "(\n"
            "    const Field<Type>& pField\n"
            ") const"
        )   << "Patch field does not correspond to patch faces.  Patch size: "
            << patch().size() << " field size: " << pField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tgField
    (
        new Field<Type>(globalPatch().size(), pTraits<Type>::zero)
    );
#ifdef OPENFOAM_NOT_EXTEND
    Field<Type>& gField = tgField.ref();
#else
    Field<Type>& gField = tgField();
#endif

    if (Pstream::parRun())
    {
        const labelList& addr = faceToGlobalAddr();

        forAll (addr, i)
        {
            gField[addr[i]] = pField[i];
        }

        // Global comm
        reduce(gField, FieldSumOp<Type>());
    }
    else
    {
        gField = pField;
    }

    return tgField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::globalPolyPatch::globalFaceToPatch
(
    const Field<Type>& gField
) const
{
    if (gField.size() != globalPatch().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > globalPolyPatch::globalFaceToPatch\n"
            "(\n"
            "    const Field<Type>& gField\n"
            ") const"
        )   << "Patch field does not correspond to global patch faces.  "
            << "Global patch size: " << globalPatch().size()
            << " field size: " << gField.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tpField
    (
        new Field<Type>(patch().size(), pTraits<Type>::zero)
    );
#ifdef OPENFOAM_NOT_EXTEND
    Field<Type>& pField = tpField.ref();
#else
    Field<Type>& pField = tpField();
#endif

    if (Pstream::parRun())
    {
        const labelList& addr = faceToGlobalAddr();

        forAll(addr, i)
        {
            pField[i] = gField[addr[i]];
        }
    }
    else
    {
        pField = gField;
    }

    return tpField;
}


// ************************************************************************* //
