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

Class
    Foam::rockFluidModels::coreyModel

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    coreyModel.C

\*---------------------------------------------------------------------------*/

#ifndef coreyModel_H
#define coreyModel_H

#include "rockFluidModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace rockFluidModels
{

/*---------------------------------------------------------------------------*\
                           Class coreyModel Declaration
\*---------------------------------------------------------------------------*/

class coreyModel
:
    public rockFluidModel
{
    // Private data

        dictionary coreyModelCoeffs_;

        dimensionedScalar krww_;
        dimensionedScalar nw_;
        dimensionedScalar krnn_;
        dimensionedScalar nn_;

        volScalarField krw_;
        volScalarField krn_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcKrw() const;
        
        tmp<volScalarField> calcKrn() const;


public:

    //- Runtime type information
    TypeName("coreyModel");


    // Constructors

        //- Construct from components
        coreyModel
        (
            const word& name,
            const dictionary& rockFluidProperties,
            const volScalarField& Sw
        );


    //- Destructor
    ~coreyModel()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> krw() const
        {
            return krw_;
        }
        
        tmp<volScalarField> krn() const
        {
            return krn_;
        }

        //- Correct the relative permeability
        void correct()
        {
            krw_ = calcKrw();
            krn_ = calcKrn();
        }

        //- Read rockFluidProperties dictionary
        bool read(const dictionary& rockFluidProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace rockFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
