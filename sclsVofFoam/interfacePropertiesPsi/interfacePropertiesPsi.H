/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::interfacePropertiesPsi

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfacePropertiesPsi.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePropertiesPsi_H
#define interfacePropertiesPsi_H

#include "IOdictionary.H"
#include "surfaceTensionModelPsi.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfacePropertiesPsi Declaration
\*---------------------------------------------------------------------------*/

class interfacePropertiesPsi
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        autoPtr<surfaceTensionModelPsi> sigmaPtr_;

        //- Parameter for S-CLSVOF method 
        dimensionedScalar deltaX_;
        dimensionedScalar gammaCoeff_;
        dimensionedScalar gamma_;
        dimensionedScalar epsilonCoeff_;
        dimensionedScalar epsilon_;
        dimensionedScalar deltaTauCoeff_;
        dimensionedScalar deltaTau_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& psi_;
        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nHatf_;
        volScalarField psi0_;
        volScalarField delta_;
        volScalarField K_;
        volScalarField H_;
        volScalarField Hscale_;
        volScalarField deltaScale_;


    // Private Member Functions

        //- No copy construct
        interfacePropertiesPsi(const interfacePropertiesPsi&) = delete;

        //- No copy assignment
        void operator=(const interfacePropertiesPsi&) = delete;

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();

        //- Re-calculate the delta and Heviside function
        void calculateDelta();
        void calculateH();
        void calculateHscale();
        void calculateDeltaScale();

        //- Re-calculate the psi0
        void calculatePsi0();


public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfacePropertiesPsi
        (
            const volScalarField& psi,
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const dimensionedScalar& epsilon() const
        {
            return epsilon_;
        }

        const dimensionedScalar& deltaTau() const
        {
            return deltaTau_;
        }

        const volScalarField& psi0() const
        {
            return psi0_;
        }

//        const volScalarField& delta() const
//        {
//            return delta_;
//        }

        tmp<volScalarField> sigmaK() const;

        tmp<volScalarField> delta() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        tmp<surfaceScalarField> surfaceTensionDensityScaledForce() const;

        tmp<surfaceScalarField> surfaceTensionDensityScaledBalancedForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        void correct();

        void correctPsi0();

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
