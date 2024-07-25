/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "surfaceTensionModelPsi.H"
#include "constantSurfaceTensionPsi.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceTensionModelPsi> Foam::surfaceTensionModelPsi::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    if (dict.isDict("sigma"))
    {
        const dictionary& sigmaDict = surfaceTensionModelPsi::sigmaDict(dict);

        const word modelType(sigmaDict.get<word>("type"));

        Info<< "Selecting surfaceTensionModelPsi " << modelType << endl;

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown surfaceTensionModelPsi type "
                << modelType << nl << nl
                << "Valid surfaceTensionModelPsi types :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return cstrIter()(sigmaDict, mesh);
    }

    return autoPtr<surfaceTensionModelPsi>
    (
        new surfaceTensionModelsPsi::constant(dict, mesh)
    );
}


// ************************************************************************* //
