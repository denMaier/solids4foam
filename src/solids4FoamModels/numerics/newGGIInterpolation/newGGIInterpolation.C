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
    Interpolation class dealing with transfer of data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "newGGIInterpolationTemplate.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
label newGGIInterpolation<MasterPatch, SlavePatch>::parMasterStart() const
{
    if (globalData())
    {
        // Integer division intended
        return Foam::min
        (
            masterPatch_.size(),
            Pstream::myProcNo()*(masterPatch_.size()/Pstream::nProcs() + 1)
        );
    }
    else
    {
        // No parallel search: do complete patch
        return 0;
    }
}


template<class MasterPatch, class SlavePatch>
label newGGIInterpolation<MasterPatch, SlavePatch>::parMasterEnd() const
{
    if (globalData())
    {
        // Integer division intended
        return Foam::min
        (
            masterPatch_.size(),
            (Pstream::myProcNo() + 1)*
            (masterPatch_.size()/Pstream::nProcs() + 1)
        );
    }
    else
    {
        // No parallel search: do complete patch
        return masterPatch_.size();
    }
}


template<class MasterPatch, class SlavePatch>
label newGGIInterpolation<MasterPatch, SlavePatch>::parMasterSize() const
{
    return Foam::max
    (
        0,
        this->parMasterEnd() - this->parMasterStart()
    );
}


template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::clearOut()
{
    deleteDemandDrivenData(masterAddrPtr_);
    deleteDemandDrivenData(masterWeightsPtr_);
    deleteDemandDrivenData(slaveAddrPtr_);
    deleteDemandDrivenData(slaveWeightsPtr_);

    if (gapIntegration_)
    {
        deleteDemandDrivenData(masterNeiIntegralGapPtr_);
        deleteDemandDrivenData(slaveNeiIntegralGapPtr_);
        deleteDemandDrivenData(masterNeiContactAreaPtr_);
    }

    deleteDemandDrivenData(uncoveredMasterAddrPtr_);
    deleteDemandDrivenData(uncoveredSlaveAddrPtr_);

    deleteDemandDrivenData(masterPointAddressingPtr_);
    deleteDemandDrivenData(masterPointWeightsPtr_);
    deleteDemandDrivenData(masterPointDistancePtr_);
    deleteDemandDrivenData(masterPointDistanceVectorsPtr_);
    masterEdgeLoopsMap_.clear();

    deleteDemandDrivenData(slavePointAddressingPtr_);
    deleteDemandDrivenData(slavePointWeightsPtr_);
    deleteDemandDrivenData(slavePointDistancePtr_);
    deleteDemandDrivenData(slavePointDistanceVectorsPtr_);
    slaveEdgeLoopsMap_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class MasterPatch, class SlavePatch>
newGGIInterpolation<MasterPatch, SlavePatch>::newGGIInterpolation
(
    const MasterPatch& masterPatch,
    const SlavePatch&  slavePatch,
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep,
    const bool globalData,
    const scalar masterNonOverlapFaceTol,
    const scalar slaveNonOverlapFaceTol,
    const bool rescaleGGIWeightingFactors,
    const quickReject reject,
    const boundBox regionOfInterest
)
:
    masterPatch_(masterPatch),
    slavePatch_(slavePatch),
    forwardT_(forwardT),
    reverseT_(reverseT),
    forwardSep_(forwardSep),
    globalData_(globalData),
    masterNonOverlapFaceTol_(masterNonOverlapFaceTol),
    slaveNonOverlapFaceTol_(slaveNonOverlapFaceTol),
    rescaleGGIWeightingFactors_(rescaleGGIWeightingFactors),
    gapIntegration_(false),
    reject_(reject),
    usePrevCandidateMasterNeighbors_(false),
    prevCandidateMasterNeighbors_(0),
    regionOfInterest_(regionOfInterest),
    masterAddrPtr_(NULL),
    masterWeightsPtr_(NULL),
    masterPointAddressingPtr_(NULL),
    masterPointWeightsPtr_(NULL),
    masterPointDistancePtr_(NULL),
    masterPointDistanceVectorsPtr_(NULL),
    masterEdgeLoopsMap_(0),
    masterNeiIntegralGapPtr_(NULL),
    slaveNeiIntegralGapPtr_(NULL),
    masterNeiContactAreaPtr_(NULL),
    slaveAddrPtr_(NULL),
    slaveWeightsPtr_(NULL),
    slavePointAddressingPtr_(NULL),
    slavePointWeightsPtr_(NULL),
    slavePointDistancePtr_(NULL),
    slavePointDistanceVectorsPtr_(NULL),
    slaveEdgeLoopsMap_(0),
    useNewPointDistanceMethod_(true),
    projectPointsToPatchBoundary_(false),
    checkPointDistanceOrientations_(useNewPointDistanceMethod_ ? false : true),
    uncoveredMasterAddrPtr_(NULL),
    uncoveredSlaveAddrPtr_(NULL)
{
    // Check size of transform.  They should be equal to slave patch size
    // if the transform is not constant
    if (forwardT_.size() > 1 || reverseT_.size() > 1)
    {
        if
        (
            forwardT_.size() != slavePatch_.size()
         || reverseT_.size() != masterPatch_.size()
        )
        {
            FatalErrorIn
            (
                "newGGIInterpolation<MasterPatch, SlavePatch>::"
                "newGGIInterpolation"
            )   << "Incorrectly defined transform: forwardT: "
                << forwardT_.size() << " patch: " << slavePatch_.size()
                << " reverseT: " << reverseT_.size()
                << " patch: " << masterPatch_.size()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
newGGIInterpolation<MasterPatch, SlavePatch>::~newGGIInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const labelListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterAddr() const
{
    if (!masterAddrPtr_)
    {
        calcAddressing();
    }

    return *masterAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterWeights() const
{
    if (!masterWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveAddr() const
{
    if (!slaveAddrPtr_)
    {
        calcAddressing();
    }

    return *slaveAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveWeights() const
{
    if (!slaveWeightsPtr_)
    {
        calcAddressing();
    }

    return *slaveWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterNeiPenVol() const
{
    if (gapIntegration_)
    {
        if (!masterNeiIntegralGapPtr_)
        {
            calcAddressing();
        }
    }
    else
    {
        FatalErrorIn
        (
            "const scalarListList& newGGIInterpolation<MasterPatch, SlavePatch>"
            "::masterNeiPenVol() const"
        )   << "Pointer not set if gapIntegration is false"
            << abort(FatalError);
    }

    return *masterNeiIntegralGapPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveNeiPenVol() const
{
    if (gapIntegration_)
    {
        if (!slaveNeiIntegralGapPtr_)
        {
            calcAddressing();
        }
    }
    else
    {
        FatalErrorIn
        (
            "const scalarListList& newGGIInterpolation<MasterPatch, SlavePatch>"
            "::slaveNeiPenVol() const"
        )   << "Pointer not set if gapIntegration is false"
            << abort(FatalError);
    }

    return *slaveNeiIntegralGapPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterNeiAreaInContact() const
{
    if (gapIntegration_)
    {
        if (!masterNeiContactAreaPtr_)
        {
            calcAddressing();
        }
    }
    else
    {
        FatalErrorIn
        (
            "const scalarListList& newGGIInterpolation<MasterPatch, SlavePatch>"
            "::masterNeiAreaInContact() const"
        )   << "Pointer not set if gapIntegration is false"
            << abort(FatalError);
    }

    return *masterNeiContactAreaPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelList&
newGGIInterpolation<MasterPatch, SlavePatch>::uncoveredMasterFaces() const
{
    if (!uncoveredMasterAddrPtr_)
    {
        calcAddressing();
    }

    return *uncoveredMasterAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelList&
newGGIInterpolation<MasterPatch, SlavePatch>::uncoveredSlaveFaces() const
{
    if (!uncoveredSlaveAddrPtr_)
    {
        calcAddressing();
    }

    return *uncoveredSlaveAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
bool newGGIInterpolation<MasterPatch, SlavePatch>::movePoints
(
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep
)
{
    this->forwardT_ = forwardT;
    this->reverseT_ = reverseT;
    this->forwardSep_ = forwardSep;

    if (prevCandidateMasterNeighbors_.size() > 0)
    {
        if (prevCandidateMasterNeighbors_.size() != parMasterSize())
        {
            Info<< "    " << typeName
                << " : clearing prevCandidateMasterNeighbors" << endl;
            clearPrevCandidateMasterNeighbors();
        }
    }

    clearOut();

    return true;
}

template<class MasterPatch, class SlavePatch>
const Foam::List<labelPair>&
newGGIInterpolation<MasterPatch, SlavePatch>::masterPointAddr() const
{
    if (!masterPointAddressingPtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointAddressingPtr_;
}

template<class MasterPatch, class SlavePatch>
const Foam::FieldField<Field, scalar>&
newGGIInterpolation<MasterPatch, SlavePatch>::masterPointWeights() const
{
    if (!masterPointWeightsPtr_)
    {
        calcMasterPointWeights();
    }

    return *masterPointWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarField&
newGGIInterpolation<MasterPatch, SlavePatch>
::masterPointDistanceToIntersection() const
{
    if (!masterPointDistancePtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointDistancePtr_;
}


template<class MasterPatch, class SlavePatch>
const vectorField&
newGGIInterpolation<MasterPatch, SlavePatch>
::masterPointDistanceVectorsToIntersection() const
{
    if (!masterPointDistanceVectorsPtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointDistanceVectorsPtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::List<labelPair>&
newGGIInterpolation<MasterPatch, SlavePatch>::slavePointAddr() const
{
    if (!slavePointAddressingPtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointAddressingPtr_;
}

template<class MasterPatch, class SlavePatch>
const Foam::FieldField<Field, scalar>&
newGGIInterpolation<MasterPatch, SlavePatch>::slavePointWeights() const
{
    if (!slavePointWeightsPtr_)
    {
        calcSlavePointWeights();
    }

    return *slavePointWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarField&
newGGIInterpolation<MasterPatch, SlavePatch>
::slavePointDistanceToIntersection() const
{
    if (!slavePointDistancePtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointDistancePtr_;
}


template<class MasterPatch, class SlavePatch>
const vectorField&
newGGIInterpolation<MasterPatch, SlavePatch>
::slavePointDistanceVectorsToIntersection() const
{
    if (!slavePointDistanceVectorsPtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointDistanceVectorsPtr_;
}

template<class MasterPatch, class SlavePatch>
const tmp<scalarField>
newGGIInterpolation<MasterPatch, SlavePatch>::masterFacePenVol() const
{
    if (!gapIntegration_)
    {
        FatalErrorIn
        (
            "const tmp<scalarField>  newGGIInterpolation"
            "<MasterPatch, SlavePatch>::masterFacePenVol() const"
        )   << "Not available if gapIntegration is false"
            << abort(FatalError);
    }

    tmp<scalarField> tresult
    (
        new scalarField
        (
            masterPatch_.size(),
            pTraits<scalar>::zero
        )
    );

    scalarField& result = tresult();

    const labelListList& ma = masterAddr();
    const scalarListList& mapv = masterNeiPenVol();

    forAll(result, faceI)
    {
        const labelList& maFaceI = ma[faceI];
        const scalarList& mapvFaceI = mapv[faceI];

        if(maFaceI.size() > 0)
        {
            forAll(mapvFaceI, I)
            {
                result[faceI] += mapvFaceI[I];
            }
        }
    }

    return tresult;
}

template<class MasterPatch, class SlavePatch>
const tmp<scalarField>
newGGIInterpolation<MasterPatch, SlavePatch>::slaveFacePenVol() const
{
    if (!gapIntegration_)
    {
        FatalErrorIn
        (
            "const tmp<scalarField>  newGGIInterpolation"
            "<MasterPatch, SlavePatch>::slaveFacePenVol() const"
        )   << "Not available if gapIntegration is false"
            << abort(FatalError);
    }

    tmp<scalarField> tresult
    (
        new scalarField
        (
            slavePatch_.size(),
            pTraits<scalar>::zero
        )
    );

    scalarField& result = tresult();

    const labelListList& sa = slaveAddr();
    const scalarListList& sapv = slaveNeiPenVol();

    forAll(result, faceI)
    {
        const labelList& saFaceI = sa[faceI];
        const scalarList& sapvFaceI = sapv[faceI];

        if(saFaceI.size() > 0)
        {
            forAll(sapvFaceI, I)
            {
                result[faceI] += sapvFaceI[I];
            }
        }
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
const tmp<scalarField>
newGGIInterpolation<MasterPatch, SlavePatch>::masterAreaInContact() const
{
    if (!gapIntegration_)
    {
        FatalErrorIn
        (
            "const tmp<scalarField>  newGGIInterpolation"
            "<MasterPatch, SlavePatch>::masterAreaInContact() const"
        )   << "Not available if gapIntegration is false"
            << abort(FatalError);
    }

    tmp<scalarField> tresult
    (
        new scalarField
        (
            masterPatch_.size(),
            pTraits<scalar>::zero
        )
    );

    scalarField& result = tresult();

    const scalarListList& mnrca = masterNeiAreaInContact();

    forAll(result, faceI)
    {
        const scalarList& mnrcaFaceI = mnrca[faceI];

        if(mnrcaFaceI.size() > 0)
        {
            forAll(mnrcaFaceI, I)
            {
                result[faceI] += mnrcaFaceI[I];
            }
        }
    }

    return tresult;
}

template<class MasterPatch, class SlavePatch>
const tmp<scalarField>
newGGIInterpolation<MasterPatch, SlavePatch>::slaveAreaInContact() const
{
    if (!gapIntegration_)
    {
        FatalErrorIn
        (
            "const tmp<scalarField>  newGGIInterpolation"
            "<MasterPatch, SlavePatch>::slaveAreaInContact() const"
        )   << "Not available if gapIntegration is false"
            << abort(FatalError);
    }

    tmp<scalarField> tresult
    (
        new scalarField
        (
            slavePatch_.size(),
            pTraits<scalar>::zero
        )
    );

    scalarField& result = tresult();

    const scalarListList& mnrca = masterNeiAreaInContact();
    const labelListList& ma = masterAddr();

    forAll(ma, faceI)
    {
        const labelList& curMa = ma[faceI];
        const scalarList& curMnrca = mnrca[faceI];

        forAll(curMa, mpI)
        {
            const label& slaveFace = curMa[mpI];
            result[slaveFace] += curMnrca[mpI];
        }
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> > newGGIInterpolation<MasterPatch, SlavePatch>::
slaveToMasterPointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != this->slavePatch().nPoints())
    {
        FatalErrorIn
        (
            "newGGIInterpolation::slaveToMasterPointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << this->slavePatch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            this->masterPatch().nPoints(),
            pTraits<Type>::zero
        )
    );

    // Escape the interpolation if there are no faces in the target patch
    if (this->masterPatch().nPoints() == 0)
    {
        return tresult;
    }

    Field<Type>& result = tresult();

    const List<typename SlavePatch::FaceType>& slaveFaces =
        this->slavePatch().localFaces();

    const List<labelPair>& addr = masterPointAddr();

    const FieldField<Field, scalar>& weights = masterPointWeights();

    forAll (result, pointI)
    {
        if (addr[pointI].first() > -1)
        {
            const face& hitFace =
                slaveFaces[addr[pointI].first()];

            label pI = addr[pointI].second();

            Type ctrF = average(Field<Type>(pf, hitFace));

            result[pointI] =
                weights[pointI][0]*pf[hitFace[pI]]
              + weights[pointI][1]*pf[hitFace.nextLabel(pI)]
              + weights[pointI][2]*ctrF;
        }
        else
        {
            FatalErrorIn
            (
                "newGGIInterpolation::masterToSlavePointInterpolate"
                "(const Field<Type> pf)"
            )   << "Master point addressing is not correct"
                << abort(FatalError);
        }
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> > newGGIInterpolation<MasterPatch, SlavePatch>::
masterToSlavePointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != this->masterPatch().nPoints())
    {
        FatalErrorIn
        (
            "newGGIInterpolation::masterToSlavePointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << this->masterPatch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            this->slavePatch().nPoints(),
            pTraits<Type>::zero
        )
    );

    // Escape the interpolation if there are no faces in the target patch
    if (this->slavePatch().nPoints() == 0)
    {
        return tresult;
    }

    Field<Type>& result = tresult();

    const List<typename SlavePatch::FaceType>& masterFaces =
        this->masterPatch().localFaces();

    const List<labelPair>& addr = slavePointAddr();

    const FieldField<Field, scalar>& weights = slavePointWeights();

    forAll (result, pointI)
    {
        if (addr[pointI].first() > -1)
        {
            const face& hitFace =
                masterFaces[addr[pointI].first()];

            label pI = addr[pointI].second();

            Type ctrF = average(Field<Type>(pf, hitFace));

            result[pointI] =
                weights[pointI][0]*pf[hitFace[pI]]
              + weights[pointI][1]*pf[hitFace.nextLabel(pI)]
              + weights[pointI][2]*ctrF;
        }
        else
        {
            FatalErrorIn
            (
                "newGGIInterpolation::masterToSlavePointInterpolate"
                "(const Field<Type> pf)"
            )   << "Slave point addressing is not correct"
                << abort(FatalError);
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newGGIInterpolationPolygonIntersection.C"
#   include "newGGIInterpolationQuickRejectTests.C"
#   include "newGGIInterpolationWeights.C"
#   include "newGGIInterpolate.C"

// ************************************************************************* //
