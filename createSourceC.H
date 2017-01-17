// BEGIN - ADC MANU
//  	C = solute concentration
Info << "Reading field C \n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading Diffusion coefficient DC\n" << endl;
dimensionedScalar DC
    (
        transportProperties.lookup("DC")
    );

// END - ADC MANU

Info<< "Reading Source Term C\n" << endl;

IOdictionary sourceProperties
(
    IOobject
    (
        "sourceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//scalar activateSource(sourceProperties.lookupOrDefault<scalar>("activateSource",0.)); 

/*
Info<< "Reading field Src (source term C)\n" << endl;
volScalarField Src
(
    IOobject
    (
        "Src",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimensionSet(1,-3,-1,0,0),0)
);
*/

//  Read the position of the excluded cell
Info<< "Reading cell position of fixed cells\n" << endl;
scalarField vals_
(
    sourceProperties.lookup("fixedCValues")
);

//  Read the value of the excluded cell
Info<< "Reading cell position of fixed cells\n" << endl;
vectorField position_
(
    sourceProperties.lookup("fixedCPositions")
);

// ADC MANU
// in controlDict
const bool adjustTimeStepCo =
    runTime.controlDict().lookupOrDefault("adjustTimeStepCo", false);

scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
scalar maxDeltaTCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaTCo", GREAT);

scalar tolC = runTime.controlDict().lookupOrDefault<scalar>("tolC", 1e-6);
scalar tolh = runTime.controlDict().lookupOrDefault<scalar>("tolh", 1e-6);
scalar maxC = runTime.controlDict().lookupOrDefault<scalar>("maxC", 1e10);
