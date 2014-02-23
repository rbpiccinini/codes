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

Namespace
    Foam::rockFluidModels

Description
    A namespace for various incompressible rockFluidModel implementations.

Class
    Foam::rockFluidModel

Description
    An abstract base class for rockFluidModels.


SourceFiles
    rockFluidModel.C
    rockFluidModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef rockFluidModel_H
#define rockFluidModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class rockFluidModel Declaration
\*---------------------------------------------------------------------------*/

class rockFluidModel
{

protected:

    // Protected data

        word name_;
        dictionary rockFluidProperties_;

        const surfaceScalarField& Swf_;
        const surfaceScalarField& rhowf_;
        const surfaceScalarField& rhonf_;
        const surfaceScalarField& muwf_;
        const surfaceScalarField& munf_;


    // Private Member Functions

        //- Disallow copy construct
        rockFluidModel(const rockFluidModel&);

        //- Disallow default bitwise assignment
        void operator=(const rockFluidModel&);


public:

    //- Runtime type information
    TypeName("rockFluidModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rockFluidModel,
            dictionary,
            (
                const word& name,
                const dictionary& rockFluidProperties,
                const surfaceScalarField& Swf,
                const surfaceScalarField& rhowf,
                const surfaceScalarField& rhonf,
                const surfaceScalarField& muwf,
                const surfaceScalarField& munf
            ),
            (name, rockFluidProperties, Swf, rhowf, rhonf, muwf, munf)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<rockFluidModel> New
        (
            const word& name,
            const dictionary& rockFluidProperties,
            const surfaceScalarField& Swf,
            const surfaceScalarField& rhowf,
            const surfaceScalarField& rhonf,
            const surfaceScalarField& muwf,
            const surfaceScalarField& munf
        );


    // Constructors

        //- Construct from components
        rockFluidModel
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
    virtual ~rockFluidModel()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& rockFluidProperties() const
        {
            return rockFluidProperties_;
        }

        //- Return the relative permeability of wetting phase
        virtual tmp<surfaceScalarField> krw() const = 0;

        //- Return the relative permeability of non-wetting phase
        virtual tmp<surfaceScalarField> krn() const = 0;

	 //- Return the relative permeability derivative of wetting phase
        virtual tmp<surfaceScalarField> diffkrw() const = 0;
        
        //- Return the relative permeability derivative of non-wetting phase
        virtual tmp<surfaceScalarField> diffkrn() const = 0;

        //- Return the relative permeability derivative of non-wetting phase
        virtual tmp<surfaceScalarField> Fw() const = 0;

        //- Return the relative permeability derivative of non-wetting phase
        virtual tmp<surfaceScalarField> Lw() const = 0;

        //- Return the relative permeability derivative of non-wetting phase
        virtual tmp<surfaceScalarField> dFw() const = 0;
        
        //- Return the total mobility 
        virtual tmp<surfaceScalarField> M() const = 0;
        
        //- Return the gravity source term
        virtual tmp<surfaceScalarField> L() const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& rockFluidProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
