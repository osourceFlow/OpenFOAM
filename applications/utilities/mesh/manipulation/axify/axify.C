/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    axify

Description
    Transforms a planar 2d mesh in xy-plane into wedge type mesh. 
    Turns the mesh around y-axis.
    
    Does not remove the faces at the axis or change boundary types. Recommended
    to run with collapseEdges and changeDictionary. 
    
    Example bash script:
    
    foamCleanPolyMesh >& log/foamCleanPolyMesh
    blockMesh >& log/blockMesh
    axify >& log/axify
    collapseEdges -overwrite >& log/collapseEdges 
    changeDictionary -constant >& log/changeDictionary
    checkMesh >& log/checkMesh

Usage
    axify

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"
#include "mathematicalConstants.H"
#include "polyMeshFilter.H"
#include "polyTopoChange.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
//    argList::addOption
//    (
//        "plane",
//        "vector",
//        "mid plane <vector> - eg, '(1 1 0)'"
//    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    fileName meshDir;

    if (args.optionReadIfPresent("region", regionName))
    {
        meshDir = regionName/polyMesh::meshSubDir;
    }
    else
    {
        meshDir = polyMesh::meshSubDir;
    }

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // this is not actually stringent enough:
//    if (args.options().empty())
//    {
//        FatalErrorInFunction
//            << "No plane suplied."

//            << exit(FatalError);
//    }


    boundBox bb(points);

    Info<< "bounding box: min = " << bb.min()
        << " max = " << bb.max() << " metres."
        << endl;


    scalar midPointZ = gAverage(points).component(vector::Z);
    
    scalar tanAlpha = 0.043660942908512058;//tan(2.5/180.0*pi);
    
    forAll(points, pointi)
    {
        scalar newZ = points[pointi].component(vector::X) * tanAlpha;
        if (points[pointi].component(vector::Z) >= midPointZ)
        {
            points[pointi].component(vector::Z) = newZ;
        }
        else if (points[pointi].component(vector::Z) < midPointZ)
        {
            points[pointi].component(vector::Z) = -newZ;
        }
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();


    // TODO!
    // remove extra faces
    // TODO!
    // change type to wedge
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
