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

\*----------------------------------------------------------------------------*/

#include "solidPointStress.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPointStress, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPointStress,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPointStress::writeData()
{
    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        // Read the stress field
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        // Create a point mesh
        const pointMesh& pMesh = pointMesh::New(mesh);

        // Create a point stress field
        pointSymmTensorField pointSigma
        (
            IOobject
            (
                "pointSigma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("zero", sigma.dimensions(), symmTensor::zero)
        );

        // Lookup the solidModel object
        const solidModel& solMod = lookupSolidModel(mesh);

        // Interpolate vol field to point field
        solMod.mechanical().volToPoint().interpolate(sigma, pointSigma);

        symmTensor pointSigmaValue = symmTensor::zero;
        if (pointID_ > -1)
        {
            pointSigmaValue = pointSigma[pointID_];
        }
        reduce(pointSigmaValue, sumOp<symmTensor>());

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                << " " << pointSigmaValue.xx()
                << " " << pointSigmaValue.xy()
                << " " << pointSigmaValue.xz()
                << " " << pointSigmaValue.yy()
                << " " << pointSigmaValue.yz()
                << " " << pointSigmaValue.zz()
                << endl;
        }
    }
    else
    {
        InfoIn(this->name() + " function object constructor")
            << "volSymmTensorField sigma not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPointStress::solidPointStress
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    pointID_(-1),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    // Lookup the point
    const vector point(dict.lookup("point"));

    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // Find the closest point
        scalar minDist = GREAT;

        forAll(mesh.points(), pI)
        {
            scalar dist = mag(mesh.points()[pI] - point);

            if (dist < minDist)
            {
                minDist = dist;
                pointID_ = pI;
            }
        }

        // Find global closest point
        const scalar globalMinDist = returnReduce(minDist, minOp<scalar>());
        int procNo = -1;
        if (mag(globalMinDist - minDist) < SMALL)
        {
            procNo = Pstream::myProcNo();
        }
        else
        {
            pointID_ = -1;
        }

        // More than one processor can have the point so we will take the proc
        // with the lowest processor number
        const int globalMinProc = returnReduce(procNo, minOp<int>());
        if (mag(globalMinProc - procNo) > SMALL)
        {
            pointID_ = -1;
        }

        if (pointID_ > -1)
        {
            Pout<< this->name()
                << ": distance from specified point is " << minDist
                << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            const word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
            (
                new OFstream(historyDir/"solidPointStress_" + name + ".dat")
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time XX XY XZ YY YZ ZZ" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPointStress::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::solidPointStress::execute(const bool forceWrite)
#else
bool Foam::solidPointStress::execute()
#endif
{
    return writeData();
}


bool Foam::solidPointStress::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidPointStress::write()
{
    return false;
}
#endif

// ************************************************************************* //
