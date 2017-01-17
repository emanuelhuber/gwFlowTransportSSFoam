//- update saturation
/*
theta == pcModel->correctAndSb(h);
*/


//h.dimensions().reset(dimless);
//theta = neg(h)*(thetamin + (thetamax - thetamin)*pow(1.0 + pow(alpha_*mag(h),n_),-m_)) + pos(h)*(thetamax - SMALL);
theta = neg(h)*(thetamin + (thetamax - thetamin)*pow(1.0 + pow(alpha_*mag(h),n_),-m_)) + pos(h)*thetamax;
//h.dimensions().reset(dimLength);

Se_ = (theta - thetamin)/(thetamax - thetamin);

//Ch_.dimensions().reset(dimless);
Ch_ = alpha_*m_*(thetamax - thetamin)/(1.0-m_)*pow(Se_,1.0/m_)*pow(1.0-pow(Se_,1.0/m_),m_) ;
//Ch_.dimensions().reset(dimless/dimLength);

//- relative permeability computation
/*
krModel->correct();
krtheta = krModel->krb();
krthetaf = fvc::interpolate(krtheta,"krtheta");
*/
//Se_ = (theta - thetamin)/(thetamax- thetamin);
krtheta =  pow(Se_,0.5) * pow(1-pow(1-pow(Se_,1/m_),m_),2);
krthetaf = fvc::interpolate(krtheta,"krtheta");

//- mobility computation 
Mf = rhotheta*mag(g)*Kf*krthetaf/mutheta;
Lf = rhotheta*Kf*krthetaf/mutheta;

//- compute fluxes
phiG = (Lf * g) & mesh.Sf();

phi == phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
Utheta = fvc::reconstruct(phi);

Utheta.correctBoundaryConditions();
forAll(phi.boundaryField(),patchi)
{
    phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
}
    