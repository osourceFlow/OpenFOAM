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

Class
    Foam::functionObjects::wallHeatFluxHTC

Group
    grpForcesFunctionObjects

Description
    Calculates and write the heat-flux and heat transfer coefficient at wall 
    patches as the volScalarField field 'wallHeatFlux' and 'htc'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    wallHeatFluxHTC1
    {
        type        wallHeatFluxHTC;
        libs        ("mylibfieldFunctionObjects.so");
        ...
        patches     (".*Wall");
        htc         true;
        Tref        300.15;
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallHeatFluxHTC | yes        |
        patches  | list of patches to process | no         | all wall patches
        htc      | calculate and write htc    | no         | true   
        Tref     | reference htc temperature  | no         | 
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::pressureTools
    Foam::functionObjects::timeControl

SourceFiles
    wallHeatFluxHTC.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatFluxHTC_H
#define functionObjects_wallHeatFluxHTC_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "HashSet.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallHeatFluxHTC Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFluxHTC
:
    public fvMeshFunctionObject,
    public logFiles
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- Calculate heat transfer coefficient, htc (defautl=true)
        bool calcHtc_;
        
        // Reference value for heat transfer coefficient (h = q / (Twall - Tref))
        dimensionedScalar Tref_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the heat-flux, compressible flow version
        void calcHeatFlux
        (
            const compressible::turbulenceModel& turbModel,
            volScalarField& wallHeatFlux,
            volScalarField& htc            
        );
        
        //- Calculate the heat-flux, incompressible flow version
        void calcHeatFlux
        (
            const incompressible::turbulenceModel& turbModel,
            volScalarField& wallHeatFlux,
            volScalarField& htc            
        );


        //- Calculate the radiation heat-flux
        void calcRadiationHeatFlux
        (
            volScalarField::Boundary& wallHeatFluxBf
        );

        //- Calculate the heat transfer coefficient, htc
        void calcHtc
        (
            volScalarField& htc,
            volScalarField::Boundary& wallHeatFluxBf
        );

private:

    // Private member functions

        //- Disallow default bitwise copy construct
        wallHeatFluxHTC(const wallHeatFluxHTC&);

        //- Disallow default bitwise assignment
        void operator=(const wallHeatFluxHTC&);


public:

    //- Runtime type information
    TypeName("wallHeatFluxHTC");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatFluxHTC
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wallHeatFluxHTC();


    // Member Functions

        //- Read the wallHeatFluxHTC data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
