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

#include "oscillateDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillateDisplacementFvPatchVectorField::oscillateDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    amplitude_(vector::zero),
    frequency_(0.0)
{}


oscillateDisplacementFvPatchVectorField::oscillateDisplacementFvPatchVectorField
(
    const oscillateDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_)
{}


oscillateDisplacementFvPatchVectorField::oscillateDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    amplitude_(dict.lookup("amplitude")),
    frequency_(readScalar(dict.lookup("frequency")))
{}

#ifndef OPENFOAM_ORG
oscillateDisplacementFvPatchVectorField::oscillateDisplacementFvPatchVectorField
(
    const oscillateDisplacementFvPatchVectorField& pivpvf
)
:
    fixedDisplacementFvPatchVectorField(pivpvf),
    amplitude_(pivpvf.amplitude_),
    frequency_(pivpvf.frequency_)
{}
#endif

oscillateDisplacementFvPatchVectorField::oscillateDisplacementFvPatchVectorField
(
    const oscillateDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(pivpvf, iF),
    amplitude_(pivpvf.amplitude_),
    frequency_(pivpvf.frequency_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void oscillateDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    totalDisp() =
        amplitude_
       *Foam::sin
        (
            2.0*mathematicalConstant::pi*frequency_*db().time().value()
        );

    fixedDisplacementFvPatchVectorField::updateCoeffs();
}


void oscillateDisplacementFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;

    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    oscillateDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
