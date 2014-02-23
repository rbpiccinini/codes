/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "coreyModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rockFluidModels
{
    defineTypeNameAndDebug(coreyModel, 0);

    addToRunTimeSelectionTable
    (
        rockFluidModel,
        coreyModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::rockFluidModels::coreyModel::calcKrw() const
{
      return krww_*pow(Sw_,nw_);
//    return dimensionedScalar("one", dimless, 1.0)*Sw_;
}

Foam::tmp<Foam::volScalarField>
Foam::rockFluidModels::coreyModel::calcKrn() const
{
      return krnn_*pow(1.0 - Sw_,nn_);
//    return dimensionedScalar("one", dimless, 1.0)*Sw_;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rockFluidModels::coreyModel::coreyModel
(
    const word& name,
    const dictionary& rockFluidProperties,
    const volScalarField& Sw
)
:
    rockFluidModel(name, rockFluidProperties, Sw),
    coreyModelCoeffs_(rockFluidProperties.subDict(typeName + "Coeffs")),
    krww_(coreyModelCoeffs_.lookup("krww")),
    nw_(coreyModelCoeffs_.lookup("nw")),
    krnn_(coreyModelCoeffs_.lookup("krnn")),
    nn_(coreyModelCoeffs_.lookup("nn")),
    krw_(calcKrw()),
    krn_(calcKrn())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::rockFluidModels::coreyModel::read
(
    const dictionary& rockFluidProperties
)
{
    rockFluidModel::read(rockFluidProperties);

    coreyModelCoeffs_ = rockFluidProperties.subDict(typeName + "Coeffs");

    coreyModelCoeffs_.lookup("krww") >> krww_;
    coreyModelCoeffs_.lookup("nw") >> nw_;
    coreyModelCoeffs_.lookup("krnn") >> krnn_;
    coreyModelCoeffs_.lookup("nn") >> nn_;

    return true;
}


// ************************************************************************* //
