/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 474d682a57dba2253bd09eb1b0ba3fe4c4896221
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void oscillatingLinearInlet_474d682a57dba2253bd09eb1b0ba3fe4c4896221(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    oscillatingLinearInletFixedValueFvPatchVectorField
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        printMessage("Construct oscillatingLinearInlet : patch/DimensionedField");
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const oscillatingLinearInletFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct oscillatingLinearInlet : patch/DimensionedField/mapper");
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct oscillatingLinearInlet : patch/dictionary");
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const oscillatingLinearInletFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        printMessage("Copy construct oscillatingLinearInlet");
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const oscillatingLinearInletFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        printMessage("Construct oscillatingLinearInlet : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oscillatingLinearInletFixedValueFvPatchVectorField::
~oscillatingLinearInletFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy oscillatingLinearInlet");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oscillatingLinearInletFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs oscillatingLinearInlet");
    }

//{{{ begin code
    const vector axis(0, 1, 0);
            vectorField v(this->patch().Cf().size());
            scalar t(this->db().time().value());
            forAll(this->patch().Cf(),patchid){
                vector vi((this->patch().Cf()[patchid].y()-0.875)/0.125*(1.0-cos(2.0*3.14159*t/5.0)),0,0);
                v[patchid] = vi;
            }
            operator==(v);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

