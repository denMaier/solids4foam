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

#include "solidStresses.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidStresses, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidStresses,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidStresses::writeData()
{
    if (patchFound_)
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

        // Calculate the arithmetic average
        symmTensor averageStress = symmTensor::zero;

        // Lookup the stress field
        if (mesh.foundObject<volSymmTensorField>("sigma"))
        {
            // Lookup the stress
            const symmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>
                (
                    "sigma"
                ).boundaryField()[historyPatchID_];

            averageStress = gAverage(sigma);
        }
        else
        {
            FatalErrorIn("bool Foam::solidStresses::writeData()")
                << "No stress tensor found!" << abort(FatalError);
        }

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                << " " << averageStress.xx()
                << " " << averageStress.xy()
                << " " << averageStress.xz()
                << " " << averageStress.yy()
                << " " << averageStress.yz()
                << " " << averageStress.zz()
                << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidStresses::solidStresses
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyPatchID_(-1),
    patchFound_(false),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    word historyPatchName("notSpecified");

    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch")
            >> historyPatchName;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "solidStresses: historyPatch not specified" << endl;
    }

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

    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);

    if (historyPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found"
            << endl;
    }
    else
    {
        patchFound_ = true;
    }

    // Create history file if not already created
    if (historyFilePtr_.empty() && patchFound_)
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
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
                    new OFstream
                    (
                        historyDir/"solidStresses" + historyPatchName + ".dat"
                    )
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

bool Foam::solidStresses::start()
{
    return writeData();
}


#if FOAMEXTEND
bool Foam::solidStresses::execute(const bool forceWrite)
#else
bool Foam::solidStresses::execute()
#endif
{
    return writeData();
}


bool Foam::solidStresses::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidStresses::write()
{
    return false;
}
#endif

// ************************************************************************* //
