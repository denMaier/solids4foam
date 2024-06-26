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
    Foam::newMapNearestAMI

Description
    Nearest-mapping Arbitrary Mesh Interface (AMI) method

SourceFiles
    newMapNearestAMI.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_ORG

#ifndef newMapNearestAMI_H
#define newMapNearestAMI_H

#include "newAMIMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class newMapNearestAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class newMapNearestAMI
:
    public newAMIMethod<SourcePatch, TargetPatch>
{
    // Private Member Functions

        // Marching front

            //- Find nearest target face for source face srcFacei
            void findNearestFace
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const label& srcFacei,
                label& tgtFacei
            ) const;

            //- Determine next source-target face pair
            void setNextNearestFaces
            (
                boolList& mapFlag,
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei
            ) const;

            //- Find mapped source face
            label findMappedSrcFace
            (
                const label tgtFacei,
                const List<DynamicList<label>>& tgtToSrc
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("newMapNearestAMI");


    // Constructors

        //- Construct from components
        newMapNearestAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );

        //- Disallow default bitwise copy construction
        newMapNearestAMI(const newMapNearestAMI&) = delete;


    //- Destructor
    virtual ~newMapNearestAMI();


    // Member Functions

        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const newMapNearestAMI&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "newMapNearestAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // end of #ifdef OPENFOAM_ORG

// ************************************************************************* //
