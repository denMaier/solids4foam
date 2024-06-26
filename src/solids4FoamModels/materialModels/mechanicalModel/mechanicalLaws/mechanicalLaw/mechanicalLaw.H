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
    Foam::mechanicalLaw

Description
    Mechanical constitutive law for the solidModels.

SourceFiles
    mechanicalLaw.C
    newMechanicalLaw.C

Author
    Philip Cardiff, UCD.
    Iago Oliveira

\*---------------------------------------------------------------------------*/

#ifndef mechanicalLaw_H
#define mechanicalLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "nonLinearGeometry.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "scalarMatrices.H"
#else
    #include "scalarSquareMatrix.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class mechanicalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Mechanical law dictionary
        mutable dictionary dict_;

        //- Name of the baseMesh region
        word baseMeshRegionName_;

        //- Nonlinear geometry enumerator
        const nonLinearGeometry::nonLinearType nonLinGeom_;

        //- Linearised shear modulus
        mutable autoPtr<volScalarField> muPtr_;

        //- Linearised shear modulus surface field
        mutable autoPtr<surfaceScalarField> mufPtr_;

        //- Linearised bulk modulus
        mutable autoPtr<volScalarField> KPtr_;

        //- Linearised bulk modulus surface field
        mutable autoPtr<surfaceScalarField> KfPtr_;

        //- Linearised implicit modulus field
        mutable autoPtr<volScalarField> impKPtr_;

        //- Initial (residual) stress field
        //  Note: for nonlinear geometry laws, it is law specific whether this
        //  is the Cauchy stress or 2nd Piola-Kirchhoff or some other stress
        //  measure
        mutable autoPtr<volSymmTensorField> sigma0Ptr_;

        //- Initial (residual) stress surface field
        mutable autoPtr<surfaceSymmTensorField> sigma0fPtr_;

        // Small strain tensor
        mutable autoPtr<volSymmTensorField> epsilonPtr_;

        // Small strain tensor surface field
        mutable autoPtr<surfaceSymmTensorField> epsilonfPtr_;

        // Total deformation gradient volField
        autoPtr<volTensorField> FPtr_;

        // Total deformation gradient surfaceField
        autoPtr<surfaceTensorField> FfPtr_;

        //- Relative deformation gradient volField
        autoPtr<volTensorField> relFPtr_;

        //- Relative deformation gradient surfaceField
        autoPtr<surfaceTensorField> relFfPtr_;

        //- Switch to enable solution of pressure Laplacian equation
        //  This can help quell oscillations in the hydrostatic stress
        const Switch solvePressureEqn_;

        //- Optional: it is possible to scale the amount of smoothing in the
        //  pressure equation with this coefficient
        const scalar pressureSmoothingScaleFactor_;

        // Hydrostatic stress (negative of hydrostatic pressure) volField
        autoPtr<volScalarField> sigmaHydPtr_;

        // Gradient of the hydrostatic stress field
        autoPtr<volVectorField> gradSigmaHydPtr_;

        //- Current time index
        label curTimeIndex_;

        //- Enforce linearity warning flag
        bool warnAboutEnforceLinear_;


    // Private Member Functions

        //- Make the mu field
        void makeMu() const;

        //- Make the muf field
        void makeMuf() const;

        //- Make the K field
        void makeK() const;

        //- Make the Kf field
        void makeKf() const;

        //- Make the sigma0 field
        void makeSigma0() const;

        //- Make the sigma0f field
        void makeSigma0f() const;

        //- Make the epsilon field
        void makeEpsilon() const;

        //- Make the epsilonf field
        void makeEpsilonf() const;

        //- Make the F field
        void makeF();

        //- Make the Ff field
        void makeFf();

        //- Make the relF field
        void makeRelF();

        //- Make the relFf field
        void makeRelFf();

        //- Make the sigmaHyd field
        void makeSigmaHyd();

        //- Disallow copy construct
        mechanicalLaw(const mechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalLaw&);


protected:

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- If the case is plane stress
        bool planeStress() const;

        //- Return const access to solvePressureEqn switch
        const Switch& solvePressureEqn() const
        {
            return solvePressureEqn_;
        }

        //- Return the base mesh region name
        const word& baseMeshRegionName() const
        {
            return baseMeshRegionName_;
        }

        //- Return const access to the base mesh
        const fvMesh& baseMesh() const
        {
            return mesh().time().lookupObject<fvMesh>
            (
                baseMeshRegionName()
            );
        }

        //- Return a const reference to the mu field
        const volScalarField& mu() const;

        //- Return a non-const reference to the mu field
        volScalarField& mu();

        //- Return a const reference to the muf field
        const surfaceScalarField& muf() const;

        //- Return a non-const reference to the muf field
        surfaceScalarField& muf();

        //- Assign a constant mu field and return a reference to it
        const volScalarField& mu(const dimensionedScalar& mu);

        //- Assign a variable mu field and return a reference to it
        const volScalarField& mu(const volScalarField& mu);

        //- Assign a constant muf field and return a reference to it
        const surfaceScalarField& muf(const dimensionedScalar& mu);

        //- Assign a variable muf field and return a reference to it
        const surfaceScalarField& muf(const surfaceScalarField& mu);

        //- Return a const reference to the K field
        const volScalarField& K() const;

        //- Return a non-const reference to the K field
        volScalarField& K();

        //- Return a const reference to the Kf field
        const surfaceScalarField& Kf() const;

        //- Return a non-const reference to the Kf field
        surfaceScalarField& Kf();

        //- Assign a constant Ku field and return a reference to it
        const volScalarField& K(const dimensionedScalar& K);

        //- Assign a variable Ku field and return a reference to it
        const volScalarField& K(const volScalarField& K);

        //- Assign a constant Kuf field and return a reference to it
        const surfaceScalarField& Kf(const dimensionedScalar& K);

        //- Assign a variable Kuf field and return a reference to it
        const surfaceScalarField& Kf(const surfaceScalarField& K);

        //- Assign a constant impK field and return a reference to it
        const volScalarField& impK(const dimensionedScalar& impK);

        //- Return a const reference to the sigma0 field
        const volSymmTensorField& sigma0() const;

        //- Return a non-const reference to the sigma0 field
        volSymmTensorField& sigma0();

        //- Return a const reference to the sigma0 surface field
        const surfaceSymmTensorField& sigma0f() const;

        //- Return a non-const reference to the sigma0 surface field
        surfaceSymmTensorField& sigma0f();

        //- Return a const reference to the epsilon field
        const volSymmTensorField& epsilon() const;

        //- Return a reference to the epsilon field
        volSymmTensorField& epsilon();

        //- Return a const reference to the epsilonf field
        const surfaceSymmTensorField& epsilonf() const;

        //- Return a reference to the Epsilonf field
        surfaceSymmTensorField& epsilonf();

        //- Return a reference to the F field
        volTensorField& F();

        //- Return a reference to the Ff field
        surfaceTensorField& Ff();

        //- Return a reference to the relF field
        volTensorField& relF();

        //- Return a reference to the relFf field
        surfaceTensorField& relFf();

        //- Return a reference to the sigmaHyd field
        volScalarField& sigmaHyd();

        //- Return a reference to the gradSigmaHyd field
        volVectorField& gradSigmaHyd();

        //- Update epsilon
        void updateEpsilon();

        //- Update epsilonf
        void updateEpsilonf();

        //- Update the deformation gradient
        //  The linearised shear and bulk modulii are used to enforce linear
        //  elasticity in the case that the enforceLinear flag is tripped by the
        //  solid model
        //  true is returned if enforceLinear is true (linear elasticity is
        //  enforced)
        bool updateF
        (
            volSymmTensorField& sigma,
            const dimensionedScalar& mu,
            const dimensionedScalar& K
        );

        //- Update the deformation gradient
        bool updateF
        (
            volSymmTensorField& sigma,
            const volScalarField& mu,
            const volScalarField& K
        );

        //- Equivalent to the updateF function, except instead for the Ff
        //  surface field
        bool updateF
        (
            surfaceSymmTensorField& sigma,
            const dimensionedScalar& mu,
            const dimensionedScalar& K
        );

        //- Equivalent to the updateF function, except instead for the Ff
        //  surface field
        bool updateF
        (
            surfaceSymmTensorField& sigma,
            const surfaceScalarField& mu,
            const surfaceScalarField& K
        );

        //- Update the provided hydrostatic stress field
        void updateSigmaHyd
        (
            volScalarField& sigmaHyd,
            const volScalarField& sigmaHydExplicit,
            const dimensionedScalar& impK
        );

        //- Update the "sigmaHyd()" hydrostatic stress field
        void updateSigmaHyd
        (
            const volScalarField& sigmaHydExplicit,
            const dimensionedScalar& impK
        );

        //- Update "sigmaHyd()" and "gradSigmaHyd()"
        void updateSigmaHyd
        (
            const volScalarField& sigmaHydExplicit,
            const volScalarField& impK
        );

        //- Update the provided hydrostatic stress field and its gradient
        void updateSigmaHyd
        (
            volScalarField& sigmaHyd,
            volVectorField& gradSigmaHyd,
            const volScalarField& sigmaHydExplicit,
            const volScalarField& impK
        );

        //- Lookup the enforceLinear Switch in the solidModel
        const Switch& enforceLinear() const;

        //- Does the solver use an incremental approach
        //  i.e. does it solve for DD as opposed to D
        bool incremental() const;

        //- Return nonlinear geometry enumerator
        nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinGeom_;
        }

public:

    //- Runtime type information
    TypeName("mechanicalLaw");


    // Declare run-time constructor selection table

        //- Mechanical law for linear geometry i.e. small strains and rotations
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            linGeomMechLaw,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const nonLinearGeometry::nonLinearType& nonLinGeom
            ),
            (name, mesh, dict, nonLinGeom)
        );

        //- Mechanical law for nonlinear geometry i.e. finite strains and
        //  rotations
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            nonLinGeomMechLaw,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const nonLinearGeometry::nonLinearType& nonLinGeom
            ),
            (name, mesh, dict, nonLinGeom)
        );


    // Selectors

        //- Create a linear geometry mechanical law
        static autoPtr<mechanicalLaw> NewLinGeomMechLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );

        //- Create a nonlinear geometry mechanical law
        static autoPtr<mechanicalLaw> NewNonLinGeomMechLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Constructors

        //- Construct from dictionary
        mechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~mechanicalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return const reference to dict
        dictionary& dict() const
        {
            return dict_;
        }

        //- Return scalar density, if defined
        virtual dimensionedScalar rhoScalar() const;

        //- Return density field
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const = 0;

        //- Return the implicit stiffness surface field
        //  This is the diffusivity for the Laplacian term
        virtual tmp<surfaceScalarField> impKf() const;

#ifdef OPENFOAM_NOT_EXTEND
        //- Material tangent for Newton Raphson loop
        virtual scalarSquareMatrix materialTangent() const
        {
            // We could use numerical differentiation for this if needed
            notImplemented("Not yet implemented for this mechanical law");

            // Keep compiler happy
            return scalarSquareMatrix(0);
        }

        //- Material tangent field for Newton Raphson loop
        //  The size of this field should be mesh.nFaces()
        virtual tmp<Field<scalarSquareMatrix>>
            materialTangentField() const;
#endif

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const
        {
            notImplemented("bulkModulus() const");

            // Keep compiler happy
            return rho();
        }

        //- Calculate the volField stress
        virtual void correct(volSymmTensorField& sigma) = 0;

        //- Calculate the surfaceField stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Calculate the pointField stress given the displacement gradient
        virtual void correct
        (
            pointSymmTensorField& sigma, const pointTensorField& gradD
        );

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update total accumulated fields
        virtual void updateTotalFields()
        {}

        //- Return the desired new time-step
        virtual scalar newDeltaT();

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
