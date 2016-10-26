/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "wallHeatFluxHTC.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFluxHTC, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFluxHTC, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFluxHTC::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall heat-flux/htc");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "q min");
    writeTabbed(file(), "q max");
    writeTabbed(file(), "q mean");
    writeTabbed(file(), "htc min");
    writeTabbed(file(), "htc max");
    writeTabbed(file(), "htc mean");
    
    file() << endl;
}


void Foam::functionObjects::wallHeatFluxHTC::calcHeatFlux
(
    const compressible::turbulenceModel& model,
    volScalarField& wallHeatFlux,
    volScalarField& htc    
)
{
    surfaceScalarField heatFlux
    (
        fvc::interpolate(model.alphaEff())*fvc::snGrad(model.transport().he())
    );

    volScalarField::Boundary& wallHeatFluxBf =
        wallHeatFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& heatFluxBf =
        heatFlux.boundaryField();

    forAll(wallHeatFluxBf, patchi)
    {
        wallHeatFluxBf[patchi] = heatFluxBf[patchi];
    }

    if (foundObject<volScalarField>("Qr"))
    {
        calcRadiationHeatFlux(wallHeatFluxBf);
    }

    if (calcHtc_)
    {
        calcHtc(htc, wallHeatFluxBf);
    }
}


void Foam::functionObjects::wallHeatFluxHTC::calcHeatFlux
(
    const incompressible::turbulenceModel& model,
    volScalarField& wallHeatFlux,
    volScalarField& htc    
)
{
    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
       "transportProperties"
    );
        
    // Prandtl number, Turbulent Prandtl number, Heat capacity,
    // Fluid density, Thermal conductivity, and Temperature.
    const dimensionedScalar Pr(transportProperties.lookup("Pr"));
    const dimensionedScalar Prt(transportProperties.lookup("Prt"));
    const dimensionedScalar Cp(transportProperties.lookup("Cp"));
    const dimensionedScalar rho(transportProperties.lookup("rho"));
    const volScalarField  alphaEff = (model.nu()/Pr+model.nut()/Prt)*Cp*rho;
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");   

    // Calculate the heat flux
    surfaceScalarField heatFlux
    (
        fvc::interpolate(alphaEff)*fvc::snGrad(T)
    );

    volScalarField::Boundary& wallHeatFluxBf =
        wallHeatFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& heatFluxBf =
        heatFlux.boundaryField();

    forAll(wallHeatFluxBf, patchi)
    {
        wallHeatFluxBf[patchi] = heatFluxBf[patchi];
    }

    if (foundObject<volScalarField>("Qr"))
    {
        calcRadiationHeatFlux(wallHeatFluxBf);
    }

    if (calcHtc_)
    {
        calcHtc(htc, wallHeatFluxBf);
    }
}


void Foam::functionObjects::wallHeatFluxHTC::calcRadiationHeatFlux
(
    volScalarField::Boundary& wallHeatFluxBf
)
{
        const volScalarField& Qr = lookupObject<volScalarField>("Qr");

        const volScalarField::Boundary& radHeatFluxBf =
            Qr.boundaryField();

        forAll(wallHeatFluxBf, patchi)
        {
            wallHeatFluxBf[patchi] += radHeatFluxBf[patchi];
        }
}


void Foam::functionObjects::wallHeatFluxHTC::calcHtc
(
    volScalarField& htc,
    volScalarField::Boundary& wallHeatFluxBf
)
{
        const volScalarField& T =
                mesh_.lookupObject<volScalarField>("T");   
        const volScalarField::Boundary& TBf = T.boundaryField();                
        volScalarField::Boundary& htcBf = htc.boundaryFieldRef();

        forAll(htcBf, patchi)
        {
            htcBf[patchi] = wallHeatFluxBf[patchi] 
                                / (TBf[patchi] - Tref_.value() + SMALL);
        }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxHTC::wallHeatFluxHTC
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    patchSet_(),
    
    // htc
    calcHtc_(false),
    Tref_(dimensionedScalar("Tref", dimTemperature, 0.0))
    
{
    volScalarField* wallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                "wallHeatFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimMass/pow3(dimTime), 0)
        )
    );

    
    volScalarField* htcPtr
    (
        new volScalarField
        (
            IOobject
            (
                "htc",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimMass/pow3(dimTime)/dimTemperature, 0)
        )
    );

    mesh_.objectRegistry::store(wallHeatFluxPtr);
    mesh_.objectRegistry::store(htcPtr);    

    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxHTC::~wallHeatFluxHTC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFluxHTC::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< "wallHeatFlux" << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }
    
    // Heat Transfer Coefficient, htc
    calcHtc_ = dict.lookupOrDefault<Switch>("htc", true);
    if (calcHtc_) 
    { 
        Tref_ = dimensionedScalar("Tref", dimTemperature, 
                                  readScalar(dict.lookup("Tref")));
    }

    return true;
}


bool Foam::functionObjects::wallHeatFluxHTC::execute()
{
    volScalarField& wallHeatFlux = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>("wallHeatFlux")
    );

    volScalarField& htc = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>("htc")
    );

    // Choose compressible or incompressible calculation.
    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
        calcHeatFlux(turbModel, wallHeatFlux, htc);
    }
    else if
    (
        foundObject<incompressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const incompressible::turbulenceModel& turbModel =
            lookupObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
        calcHeatFlux(turbModel, wallHeatFlux, htc);    
    }  
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallHeatFluxHTC::write()
{
    logFiles::write();

    const volScalarField& wallHeatFlux =
        obr_.lookupObject<volScalarField>("wallHeatFlux");
    
    const volScalarField& htc =
        obr_.lookupObject<volScalarField>("htc");    

    Log << "wallHeatFluxHTC" << " " << name() << " write:" << nl
        << "    writing fields: " << wallHeatFlux.name() << " q'' (W/m2)";

    if (calcHtc_) 
    { 
        Log << ", heat transfer coefficient " << htc.name() << " (W/m2K)";
    }
    Log << ", (min/max/mean)" << endl;

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp =
            wallHeatFlux.boundaryField()[patchi];

        // heat flux.
        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        //const scalar integralHfp = gSum(magSf[patchi]*hfp);
        const scalar meanHfp = gSum(magSf[patchi]*hfp) / gSum(magSf[patchi]);

        // htc.
        const scalarField& htcp =
            htc.boundaryField()[patchi];

        const scalar minHtcp = gMin(htcp);
        const scalar maxHtcp = gMax(htcp);
        const scalar meanHtcp = gSum(magSf[patchi]*htcp) / gSum(magSf[patchi]);


        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << meanHfp;
            if (calcHtc_) 
            {     
                file()
                    << token::TAB << minHtcp
                    << token::TAB << maxHtcp
                    << token::TAB << meanHtcp;
            }
            file() << endl;
        }

        Log << "    q'' (" << pp.name() << ")  = "
            << minHfp << ", " << maxHfp << ", " << meanHfp << endl;
        if (calcHtc_) 
        {     
            Log << "    htc (" << pp.name() << ") = "
                << minHtcp << ", " << maxHtcp  << ", " << meanHtcp 
                << endl;
        }
        Log << endl;
    }

    return true;
}


// ************************************************************************* //
