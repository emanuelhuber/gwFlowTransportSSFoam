/*
Utheta.correctBoundaryConditions();
forAll(phi.boundaryField(),patchi)
{
    phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
}
*/

Info<< "Reading thetamax\n" << endl;
dimensionedScalar thetamax
    (
        transportProperties.lookup("thetamax")
    );

Info<< "Reading thetamin\n" << endl;
dimensionedScalar thetamin
    (
        transportProperties.lookup("thetamin")
    );

Info<< "Reading Specific storate (Ss)\n" << endl;
dimensionedScalar Ss
    (
        transportProperties.lookup("Ss")
    );

const dictionary& VanGenuchtenCoeffsDict = transportProperties.subDict("VanGenuchtenCoeffs");

const scalar m_ = VanGenuchtenCoeffsDict.lookupOrDefault<scalar>("m", 1.);
const scalar n_ = 1/(1-m_);
//const scalar alpha_ = VanGenuchtenCoeffsDict.lookupOrDefault<scalar>("alpha", 1.);
dimensionedScalar alpha_
  (
    VanGenuchtenCoeffsDict.lookupOrDefault("alpha",
      dimensionedScalar("alpha",dimless,3.5))
  );

Info<< "Computing theta\n" << endl;
//h.dimensions().reset(dimless);
//theta == neg(h)*(thetamin + (thetamax - thetamin)*pow(1.0 + pow(alpha_*mag(h),n_),-m_)) + pos(h)*(thetamax - SMALL);
theta == neg(h)*(thetamin + (thetamax - thetamin)*pow(1.0 + pow(alpha_*mag(h),n_),-m_)) + pos(h)*thetamax;
//h.dimensions().reset(dimLength);
theta.write();

Info<< "Computing Se\n" << endl;
volScalarField Se_ = (theta - thetamin)/(thetamax - thetamin);

Info<< "Computing krtheta\n" << endl;
volScalarField	krtheta =  pow(Se_, 0.5) * pow(1 - pow(1 - pow(Se_, 1/m_), m_),2);

Info<< "Computing krthetaf\n" << endl;
surfaceScalarField krthetaf = fvc::interpolate(krtheta,"krtheta");

Info<< "Computing specific moisture capacity\n" << endl;
volScalarField Ch_ = alpha_*m_*(thetamax - thetamin)/(1.0-m_)*pow(Se_,1.0/m_)*pow(1.0-pow(Se_,1.0/m_),m_) ;
//Ch_.dimensions().reset(dimless/dimLength);



//- relative permeability (kr)
/*
Info<< "relative permeability (kr)\n" << endl;
krModel->correct();
volScalarField krtheta = krModel->krb();
surfaceScalarField krthetaf ("krthetaf",fvc::interpolate(krtheta,"krtheta"));
*/
//- mobility and fractional flow 
Info<< "mobility and fractional flow \n" << endl;
surfaceScalarField Mf ("Mf",rhotheta*mag(g)*Kf*krthetaf/mutheta);
surfaceScalarField Lf ("Lf",rhotheta*Kf*krthetaf/mutheta);

//- fluxes depending on saturation
Info<< "fluxes depending on saturation \n" << endl;
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());

//- Test if gravity is present
if (mag(g).value() == 0)
{
    FatalErrorIn("createthetaFields.H")
        << " Magnitude of gravity mag(g) equal to zero " << abort(FatalError);
}


//Utheta.write();
//- null Pc field for darcy velocity boundary conditions
surfaceScalarField phiPcNull("phiPc",0*phiG);


// additional log files for residuals
OFstream simuLogFile(runTime.path()/"log.residuals");

Utheta.correctBoundaryConditions();
forAll(phi.boundaryField(),patchi)
{
    phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
}
