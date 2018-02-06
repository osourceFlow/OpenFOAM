// Initialise solid field pointer lists
PtrList<coordinateSystem> coordinates(solidRegions.size());
PtrList<solidThermo> thermos(solidRegions.size());
PtrList<radiation::radiationModel> radiations(solidRegions.size());
PtrList<fv::options> solidHeatSources(solidRegions.size());
PtrList<volScalarField> fieldSolidHS(solidRegions.size());
PtrList<volScalarField> corrFieldSolidHS(solidRegions.size());
PtrList<volScalarField> betavSolid(solidRegions.size());
PtrList<volSymmTensorField> aniAlphas(solidRegions.size());
PtrList<scalar> resistivityRho1(solidRegions.size());
PtrList<scalar> resistivityRho0(solidRegions.size());
PtrList<scalar> resistivityTref(solidRegions.size());

List<bool> residualReachedSolid(solidRegions.size(), true);
List<bool> residualControlUsedSolid(solidRegions.size(), false);
List<bool> firstIterationSolid(solidRegions.size(), true);

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermos\n" << endl;
    thermos.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to radiations\n" << endl;
    radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

    Info<< "    Adding fvOptions\n" << endl;
    solidHeatSources.set
    (
        i,
        new fv::options(solidRegions[i])
    );

    fieldSolidHS.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "fieldSolidHS",
                runTime.timeName(),
                solidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    corrFieldSolidHS.set
    (
	i,
	new volScalarField
	(
	    IOobject
	    (
		"corrFieldSolidHS",
		runTime.timeName(),
		solidRegions[i],
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    solidRegions[i],
	    dimensionedScalar("corrFieldSolidHS", dimPower/dimVolume, 0.0)
	)
    );

    if (!thermos[i].isotropic())
    {
        Info<< "    Adding coordinateSystems\n" << endl;
        coordinates.set
        (
            i,
            coordinateSystem::New(solidRegions[i], thermos[i])
        );

        tmp<volVectorField> tkappaByCp =
            thermos[i].Kappa()/thermos[i].Cp();

        aniAlphas.set
        (
            i,
            new volSymmTensorField
            (
                IOobject
                (
                    "Anialpha",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedSymmTensor
                (
                    "zero",
                    tkappaByCp().dimensions(),
                    Zero
                ),
                zeroGradientFvPatchSymmTensorField::typeName
            )
        );

        aniAlphas[i].primitiveFieldRef() =
            coordinates[i].R().transformVector(tkappaByCp());
        aniAlphas[i].correctBoundaryConditions();

    }

    IOobject betavSolidIO
    (
        "betavSolid",
        runTime.timeName(),
        solidRegions[i],
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (betavSolidIO.typeHeaderOk<volScalarField>(true))
    {
        betavSolid.set
        (
            i,
            new volScalarField(betavSolidIO, solidRegions[i])
        );
    }
    else
    {
        betavSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "betavSolid",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar("1", dimless, scalar(1))
            )
        );
    }

    const dictionary& pimpleDict =
        solidRegions[i].solutionDict().subDict("PIMPLE");

    if (pimpleDict.isDict("residualControl"))
    {
        Info<< "    Reading residual controls\n" << endl;

        bool& residualControlUsed = residualControlUsedSolid[i];
        bool& resReachedSolid = residualReachedSolid[i];
        bool& firstIteration = firstIterationSolid[i];
        const fvMesh& mesh = solidRegions[i];

        #include "readSolidMultiRegionResidualControls.H"
    }

    // Read constants from thermophysicalProperties file
    dictionary subDict = thermos[i].subOrEmptyDict("electricResistivity");
    resistivityRho1.set(i, new scalar(subDict.lookupOrDefault("rho1", 0.0)));
    resistivityRho0.set(i, new scalar(subDict.lookupOrDefault("rho0", 1.0)));
    resistivityTref.set(i, new scalar(subDict.lookupOrDefault("Tref", 1.0)));
}
