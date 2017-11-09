Info<< "Reading field u\n" << endl;

//Electric potential
volScalarField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Electrical resistivity at reference temperature
dimensionedScalar rho0
(
    "rho0",
    dimensionSet(1, 3, -3, 0, 0, -2, 0),
    transportProperties.lookup("rho0")
);

// Temperature coefficient for linear approximation
dimensionedScalar alpha
(
    "alpha",
    dimensionSet(0, 0, 0, -1, 0, 0, 0),
    transportProperties.lookup("alpha")
);

// Reference temperature for linear approximation
dimensionedScalar T0
(
    "T0",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    transportProperties.lookup("T0")
);


Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Electrical conductivity
volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1/(rho0*(1+alpha*(T-T0)))
);

//Electric field
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimensionedVector("E", dimensionSet(1,1,-3,0,0,-1,0), Foam::vector(0,0,0))
);

//Current density
volVectorField J
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), Foam::vector(0,0,0)) 
);

//Joule heating per volume
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("p", dimPower/dimVolume, 0) 
);

#include "createFvOptions.H"
