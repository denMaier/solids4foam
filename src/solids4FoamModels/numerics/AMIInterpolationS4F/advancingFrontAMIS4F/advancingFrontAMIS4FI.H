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

inline const Foam::standAlonePatch& Foam::advancingFrontAMIS4F::srcPatch() const
{
    if (!tsrcPatch0_.valid())
    {
        FatalErrorInFunction
            << "tsrcPatch0_ not allocated"
            << abort(FatalError);
    }

    return tsrcPatch0_();
}


inline const Foam::standAlonePatch& Foam::advancingFrontAMIS4F::tgtPatch() const
{
    if (extendedTgtPatchPtr_)
    {
        return extendedTgtPatchPtr_();
    }

    if (!ttgtPatch0_.valid())
    {
        FatalErrorInFunction
            << "srcPatch0Ptr not allocated"
            << abort(FatalError);
    }

    return ttgtPatch0_();
}


inline const Foam::labelList& Foam::advancingFrontAMIS4F::srcNonOverlap() const
{
    return srcNonOverlap_;
}


// ************************************************************************* //
