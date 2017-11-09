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

Application
    jouleHeatingFoam

Description
    Solves Joule heating with temperature dependand conductivity.
    
    Conductivity is solved with linear interpolation as
        1/(rho0*(1+alpha*(T-T0)))
    see creteFields.H
    
    Modified from laplacianFoam by Antti Mikkonen, a.mikkonen@iki.fi, 2017

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting calculations\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Calculating potential
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix uEqn
            (
                fvm::laplacian(sigma, u)
            );
            uEqn.solve();
        }

        // Calculating electric field
        E = -fvc::grad(u);

        // Calculating current
        J = sigma*E;    
        
        // Calculating Joule heating (W/m3)
        p = magSqr(J)/sigma;

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "Info:" << endl;
    
    // Calculation total Joule Heating
    const DimensionedField<scalar, volMesh>& V = mesh.V();
    scalar totalJouleHeating = 0;
    forAll(p, celli)
    {
        totalJouleHeating += p[celli]*V[celli];
    }
    reduce(totalJouleHeating, sumOp<scalar>());
    
    // Calculation total current
    const fvPatchList& patches = mesh.boundary(); 
    scalar totalCurrent = 0;
    scalar ntotalCurrent = 0;
    forAll(patches, patchI)
    {
        const fvPatch& cPatch = patches[patchI]; 
        const vectorField& inletAreaVectors = cPatch.Sf();
        forAll(cPatch, faceI)
        {
            //calculate boundary flux
            scalar cFaceFlux
                = J.boundaryField()[patchI][faceI] & inletAreaVectors[faceI];
            if (cFaceFlux>0)
            {
                totalCurrent  += cFaceFlux;
            } 
            else 
            {
                ntotalCurrent += cFaceFlux;
            }
        } 
    }
    reduce(totalCurrent, sumOp<scalar>());
    reduce(ntotalCurrent, sumOp<scalar>());
    
    Info<< "Current out - in    = " << totalCurrent + ntotalCurrent << " A" << endl;
    Info<< "Total current       = " << totalCurrent << " A" << endl;
    Info<< "Total Joule heating = " << totalJouleHeating << " W" << endl;
    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
