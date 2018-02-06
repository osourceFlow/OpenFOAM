fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

tmp<volSymmTensorField> taniAlpha;
if (!thermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.R().transformVector(tkappaByCp());
    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}


volScalarField& h = thermo.he();

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];

const volScalarField& fHS = fieldSolidHS[i];

bool& resReachedSolid = residualReachedSolid[i];
bool& residualControlUsed = residualControlUsedSolid[i];
bool& firstIteration = firstIterationSolid[i];

const scalar& resRho1 = resistivityRho1[i];
const scalar& resRho0 = resistivityRho0[i];
const scalar& resTref = resistivityTref[i];

// Temperature is multiplied by 1.0 1/K to make volScalarfield Tm unitless
const volScalarField Tm = thermo.T()* dimensionedScalar("Temp", dimless/dimTemperature, scalar(1.0)) ;

volScalarField& corrFHS = corrFieldSolidHS[i];
corrFHS = (resRho1 * Tm  + resRho0)/(resRho1 * resTref + resRho0) * fHS;
