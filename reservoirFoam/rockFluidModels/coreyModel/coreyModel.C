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

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcKrw() const
{
	return krww_*pow
	(
		min(
                max(
                    (Swf_-Swi_) / (1.0 - Swi_ - Snr_),
                    0.0
                ),
                1.0
            ),
		nw_
	);
//    return dimensionedScalar("one", dimless, 1.0)*Swf_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcKrn() const
{
      return krnn_*pow
	(
        min(
                max(
                    (1.0 - Snr_ - Swf_) / (1.0 - Swi_ - Snr_),
                    0.0
                ),
                1.0
            ),
		nn_
	);
//    return dimensionedScalar("one", dimless, 1.0)*Swf_;

}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcDiffKrw() const
{
	scalar dS = 0.001;

	return 0.5/dS * (
		krww_*pow
		(
			(Swf_+dS-Swi_) / (1.0 - Swi_ - Snr_),
			nw_
		) -
		krww_*pow
		(
			(Swf_-dS-Swi_) / (1.0 - Swi_ - Snr_),
			nw_
		)

	);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcDiffKrn() const
{
	scalar dS = 0.001;

	return 0.5/dS * (
		krnn_*pow
		(
			(1.0 - Snr_ - Swf_ - dS) / (1.0 - Swi_ - Snr_),
			nn_
		) -
		krnn_*pow
		(
			(1.0 - Snr_ - Swf_ + dS) / (1.0 - Swi_ - Snr_),
			nn_
		)

	);

}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcFw() const
{
	return  krw_/muwf_ / (krw_/muwf_ + krn_ / munf_ ) ;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcLw() const
{
	return  krw_*rhowf_/muwf_ / ( krw_*rhowf_/muwf_ + krn_*rhonf_/munf_ ) ;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcDFw() const
{
	return ( diffkrw_*krn_ - krw_*diffkrn_) 
		/ pow(krw_*muwf_ + krn_*munf_ , 2)
		* muwf_ * munf_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcM() const
{
	return  krw_/muwf_ + krn_ / munf_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::rockFluidModels::coreyModel::calcL() const
{
	return krw_*rhowf_/muwf_ + krn_*rhonf_/munf_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rockFluidModels::coreyModel::coreyModel
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
    rockFluidModel(name, rockFluidProperties, Swf, rhowf, rhonf, muwf, munf),
    coreyModelCoeffs_(rockFluidProperties.subDict(typeName + "Coeffs")),
    krww_(coreyModelCoeffs_.lookup("krww")),
    nw_(coreyModelCoeffs_.lookup("nw")),
    krnn_(coreyModelCoeffs_.lookup("krnn")),
    nn_(coreyModelCoeffs_.lookup("nn")),
    Swi_(coreyModelCoeffs_.lookup("Swi")),
    Snr_(coreyModelCoeffs_.lookup("Snr")),
    krw_(calcKrw()),
    krn_(calcKrn()),
    diffkrw_(calcDiffKrw()),
    diffkrn_(calcDiffKrn()),
    Fw_(calcFw()),
    Lw_(calcLw()),
    dFw_(calcDFw()),
    M_(calcM()),
    L_(calcL())
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
    coreyModelCoeffs_.lookup("Swi") >> Swi_;
    coreyModelCoeffs_.lookup("Snr") >> Snr_;

    return true;
}


// ************************************************************************* //
