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
     Standard power-law  model of relative permeability.

SourceFiles
    coreyModel.C

\*---------------------------------------------------------------------------*/

#ifndef coreyModel_H
#define coreyModel_H

#include "rockFluidModel.H"
#include "volFields.H"
#include "surfaceFields.H"

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

/*
        dimensionedScalar krww_("krww", dimensionSet(0,0,0,0,0,0,0), scalar(0.0)); 
        dimensionedScalar nw_("nw", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar krnn_("krnn", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar nn_("nn", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Swi_("Swi", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
        dimensionedScalar Snr_("Snr", dimensionSet(0,0,0,0,0,0,0), scalar(0.0));
*/

        dimensionedScalar krww_;
        dimensionedScalar nw_;
        dimensionedScalar krnn_;
        dimensionedScalar nn_;
        dimensionedScalar Swi_;
        dimensionedScalar Snr_;

        surfaceScalarField krw_;
        surfaceScalarField krn_;

        surfaceScalarField diffkrw_;
        surfaceScalarField diffkrn_;

        surfaceScalarField Fw_;
        surfaceScalarField Lw_;
        surfaceScalarField dFw_;

        surfaceScalarField M_;
        surfaceScalarField L_;

    // Private Member Functions

        //- Calculate and return relative permeabilities
        tmp<surfaceScalarField> calcKrw() const;
        
        tmp<surfaceScalarField> calcKrn() const;
        
        tmp<surfaceScalarField> calcdKrw() const;
        
        tmp<surfaceScalarField> calcdKrn() const;

        //- Calculate and return derivatives of relative permeability
        tmp<surfaceScalarField> calcDiffKrw() const;
        
        tmp<surfaceScalarField> calcDiffKrn() const;
        
        //- Calculate and return derivatives of relative permeability
        tmp<surfaceScalarField> calcFw() const;

        //- Calculate and return derivatives of relative permeability
        tmp<surfaceScalarField> calcLw() const;
        
        tmp<surfaceScalarField> calcDFw() const;

        //- Calculate M and L
        tmp<surfaceScalarField> calcM() const;
        
        tmp<surfaceScalarField> calcL() const;
public:

    //- Runtime type information
    TypeName("coreyModel");


    // Constructors

        //- Construct from components
        coreyModel
        (
            const word& name,
            const dictionary& rockFluidProperties,
            const surfaceScalarField& Swf,
            const surfaceScalarField& rhowf,
            const surfaceScalarField& rhonf,
            const surfaceScalarField& muwf,
            const surfaceScalarField& munf
        );


    //- Destructor
    ~coreyModel()
    {}


    // Member Functions

        //- Return relative permeabilities
        tmp<surfaceScalarField> krw() const
        {
            return krw_;
        }
        
        tmp<surfaceScalarField> krn() const
        {
            return krn_;
        }

        tmp<surfaceScalarField> diffkrw() const
        {
            return diffkrw_;
        }
        
        tmp<surfaceScalarField> diffkrn() const
        {
            return diffkrn_;
        }

        tmp<surfaceScalarField> Fw() const
        {
            return Fw_;
        }

        tmp<surfaceScalarField> Lw() const
        {
            return Lw_;
        }
        
        tmp<surfaceScalarField> dFw() const
        {
            return dFw_;
        }

        tmp<surfaceScalarField> M() const
        {
            return M_;
        }
        
        tmp<surfaceScalarField> L() const
        {
            return L_;
        }

	//- Return relative permeability end points Swi and Snr
/*
	scalar Swi() const
	{
	    return Swi_;
	}

	scalar Snr() const
	{
	    return Snr_;
	}
*/

        //- Correct relative permeabilities
        void correct()
        {
            krw_ = calcKrw();
            krn_ = calcKrn();
            diffkrw_ = calcDiffKrw();
            diffkrn_ = calcDiffKrn();
            Fw_  = calcFw();
			Lw_ = calcLw();
            dFw_ = calcDFw();
            M_ = calcM();
            L_ = calcL();
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
