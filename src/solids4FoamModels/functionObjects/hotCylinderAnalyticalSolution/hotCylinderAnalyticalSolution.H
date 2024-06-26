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
    hotCylinderAnalyticalSolution

Description
    FunctionObject that generates the analytical temperature (analyticalT) and
    stress fields (analyticalRadialStress, analyticalHoopStress) for the case of
    a thermally-stressed pipe/cylinder.
    The solution is taken from Timoshenko, Theory of Elasticity.

    The analytical solution requires the following inputs:
        - rInner : inner pipe radius
        - rOuter : outer pipe radius
        - TInner : temperature on the inner pipe surface
        - TOuter : temperature on the outer pipe surface
        - E : Young's modulus
        - nu : Poisson's ratio
        - alpha : coefficient of linear thermal expansion

    It is assumed that the inner and outer pipe surfaces are traction-free.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    hotCylinderAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef hotCylinderAnalyticalSolution_H
#define hotCylinderAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class hotCylinderAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class hotCylinderAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Inner pipe radius
        const scalar rInner_;

        //- Outer pipe radius
        const scalar rOuter_;

        //- Inner pipe temperature
        const scalar TInner_;

        //- Outer pipe temperature
        const scalar TOuter_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Coefficient of linear thermal expansion
        const scalar alpha_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        hotCylinderAnalyticalSolution
        (
            const hotCylinderAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const hotCylinderAnalyticalSolution&);


public:

    //- Runtime type information
    TypeName("hotCylinderAnalyticalSolution");


    // Constructors

        //- Construct from components
        hotCylinderAnalyticalSolution
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

#ifdef OPENFOAM_NOT_EXTEND
        //- Write
        virtual bool write();
#endif

#ifndef OPENFOAM_NOT_EXTEND
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
