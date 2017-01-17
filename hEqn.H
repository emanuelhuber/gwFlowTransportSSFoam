{
    /*
    Utheta.correctBoundaryConditions();
    forAll(phi.boundaryField(),patchi)
    {
        phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
    */
    h.storePrevIter();

    fvScalarMatrix hEqn
        (
           (Ss*theta/thetamax + Ch_) * fvm::ddt( h)
           //pcModel->Ch() * fvm::ddt(h)
           + fvc::div(phiG)
           - fvm::laplacian(Mf,h)
          ==
           (-SrcExt*Wext+SrcInj*Winj)*activateWellbores
        );

    hEqn.solve();
         
/*
    theta.storePrevIter();
    dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
    h_ss = rDeltaT*(thetaprevIter - theta.prevIter()) +
            rDeltaT*Ss*theta/thetamax*(hprevIter - h.prevIter());
  
    fvScalarMatrix hEqn
        (
           Ch_ * fvm::ddt( h)
           + fvc::div(phiG)
           - fvm::laplacian(Mf,h)
          - h_ss
           ==
            (-SrcExt*Wext+SrcInj*Winj)*activateWellbores
        );
    hEqn.solve();
 */     
      
    resPicard = gMax((mag(h-h.prevIter()))->internalField());
    Info << "Picard iteration : " << iterPicard+1 << " residual = " << resPicard << endl; 
/*
    phi == phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
    Utheta = fvc::reconstruct(phi);
*/
}
