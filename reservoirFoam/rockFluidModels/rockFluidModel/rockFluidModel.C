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

#include "rockFluidModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rockFluidModel, 0);
    defineRunTimeSelectionTable(rockFluidModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rockFluidModel::rockFluidModel
(
    const word& name,
    const dictionary& rockFluidProperties,
    const surfaceScalarField& Swf,
    const surfaceScalarField& rhowf,
    const surfaceScalarField& rhonf,
    const surfaceScalarField& muwf,
    const surfaceScalarField& munf
)
:
    name_(name),
    rockFluidProperties_(rockFluidProperties),
    Swf_(Swf),
    rhowf_(rhowf),
    rhonf_(rhonf),
    muwf_(muwf),
    munf_(munf)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::rockFluidModel::read(const dictionary& rockFluidProperties)
{
    rockFluidProperties_ = rockFluidProperties;

    return true;
}


// ************************************************************************* //
