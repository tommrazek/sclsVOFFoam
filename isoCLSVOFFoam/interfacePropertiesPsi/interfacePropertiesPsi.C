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

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesPsi.H"
#include "psiContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfacePropertiesPsi::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradPsif
) const
{
    const fvMesh& mesh = psi_.mesh();
    const volScalarField::Boundary& abf = psi_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<psiContactAngleFvPatchScalarField>(abf[patchi]))
        {
            psiContactAngleFvPatchScalarField& acap =
                const_cast<psiContactAngleFvPatchScalarField&>
                (
                    refCast<const psiContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradPsif[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfacePropertiesPsi::calculateK()
{
    const fvMesh& mesh = psi_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of psi
    const volVectorField gradPsi(fvc::grad(psi_, "nHat"));

    // Interpolated face-gradient of psi
    surfaceVectorField gradPsif(fvc::interpolate(gradPsi));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradPsif/(mag(gradPsif) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradPsif.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */

}

void Foam::interfacePropertiesPsi::calculatePsi0()
{   
    psi0_ == (double(2.0)*alpha1_ - double(1.0))*gamma_;
}

void Foam::interfacePropertiesPsi::calculateDelta()
{   
    forAll(psi_.mesh().cells(),celli)
    {  
       if(mag(psi_[celli]) > epsilon_.value())
          delta_[celli] = double(0.0);
       else
          delta_[celli] = double(1.0)/(double(2.0)*epsilon_.value())*(double(1.0)+cos(M_PI*psi_[celli]/epsilon_.value()));
    }
}

void Foam::interfacePropertiesPsi::calculateH()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          H_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          H_[celli] = double(1.0);
       else
          H_[celli] = double(1.0)/double(2.0)*(double(1.0)+psi_[celli]/epsilon_.value()+sin(M_PI*psi_[celli]/epsilon_.value())/M_PI);
    }
}

void Foam::interfacePropertiesPsi::calculateHscale()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          Hscale_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          Hscale_[celli] = double(1.0);
       else
          Hscale_[celli] = double(1.0)/double(2.0)*(double(1.0)/double(2.0)+psi_[celli]/epsilon_.value()+psi_[celli]*psi_[celli]/(double(2.0)*epsilon_.value()*epsilon_.value())-(cos(double(2.0)*M_PI*psi_[celli]/epsilon_.value())-double(1.0))/(double(4.0)*M_PI*M_PI)+sin(M_PI*psi_[celli]/epsilon_.value())*(epsilon_.value()+psi_[celli])/(M_PI*epsilon_.value()));
    }
}

void Foam::interfacePropertiesPsi::calculateDeltaScale()
{
    deltaScale_ == double(2.0)*H_*delta_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacePropertiesPsi::interfacePropertiesPsi
(
    const volScalarField& psi,
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        psi.mesh().solverDict(psi.name()).get<scalar>("cAlpha")
    ),

    deltaX_
    (
        psi.mesh().solverDict(psi.name()).get<scalar>("deltaX")
    ),

    gammaCoeff_
    (
        psi.mesh().solverDict(psi.name()).get<scalar>("gammaCoeff")
    ),

    gamma_
    (
        "gamma",
        deltaX_*gammaCoeff_
    ),

    epsilonCoeff_
    (
        psi.mesh().solverDict(psi.name()).get<scalar>("epsilonCoeff")
    ),

    epsilon_
    (
        "epsilon",
        deltaX_*epsilonCoeff_
    ),

    deltaTauCoeff_
    (
        psi.mesh().solverDict(psi.name()).get<scalar>("deltaTauCoeff")
    ),

    deltaTau_
    (
        "deltaTau",
        deltaX_*deltaTauCoeff_
    ),

    sigmaPtr_(surfaceTensionModelPsi::New(dict, psi.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(psi.mesh().V()), 1.0/3.0)
    ),

    psi_(psi),
    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),

    K_
    (
        IOobject
        (
            "interfacePropertiesPsi:K",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    psi0_
    (   
        IOobject
        (   
            "psi0",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("psi0", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    delta_
    (
        IOobject
        (
            "delta",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("delta", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    H_
    (
        IOobject
        (
            "H",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("H", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    Hscale_
    (   
        IOobject
        (   
            "Hscale",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("Hscale", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    deltaScale_
    (
        IOobject
        (
            "deltaScale",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("deltaScale", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    )
{
    calculateK();
    calculatePsi0();
    calculateDelta();
    calculateH();
    calculateHscale();
    calculateDeltaScale();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesPsi::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}

Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesPsi::delta() const
{
    return delta_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesPsi::surfaceTensionForce() const
{
//    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
    return fvc::interpolate(sigmaK()*delta())*fvc::snGrad(psi_);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesPsi::surfaceTensionDensityScaledForce() const
{
    return fvc::interpolate(sigmaK()*deltaScale_)*fvc::snGrad(psi_);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesPsi::surfaceTensionDensityScaledBalancedForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(Hscale_);
}

Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesPsi::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::interfacePropertiesPsi::correct()
{
    calculateK();
    calculateDelta();
    delta_.correctBoundaryConditions();
    calculateH();
    H_.correctBoundaryConditions();
    calculateHscale();
    Hscale_.correctBoundaryConditions();
    calculateDeltaScale();
    deltaScale_.correctBoundaryConditions();
}

void Foam::interfacePropertiesPsi::correctPsi0()
{
    calculatePsi0();
    psi0_.correctBoundaryConditions();
}

bool Foam::interfacePropertiesPsi::read()
{
    psi_.mesh().solverDict(psi_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //
