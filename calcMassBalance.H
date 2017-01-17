/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
     

Description
     calculate mass balance    

\*---------------------------------------------------------------------------*/

  //calcuate the change of water in the domain
  /*scalar volMassChange =  (theta-thetaprevIter)().weightedAverage(mesh.V()).value()
                         *sum(mesh.V().field());

 */
/*
 scalar totalWaterCompress =  (Ss*theta/thetamax*(h-hprevIter))().weightedAverage(mesh.V()).value()
                         *sum(mesh.V().field());
*/
  //calcuate the net flux of water through the boundaries
  /*surfaceScalarField massFlux =
           - (fvc::interpolate((soil->Ks()*kr)&(fvc::grad(h+z))) & mesh.Sf())
            *runTime.deltaT();

  const surfaceScalarField::GeometricBoundaryField& patchMassFlux =
            massFlux.boundaryField();
*/
  scalar fluxBoundary = 0.0;

  forAll(phi.boundaryField(),patchi)
  {
      scalar fluxTmp = gSum(
                            phi.boundaryField()[patchi]
                            );
      fluxBoundary += fluxTmp;
  }

  /*forAll(patchMassFlux, patchi)
  {
     scalar fluxTmp = 
            gSum
            (
              patchMassFlux[patchi]
            );
      fluxBoundary += fluxTmp;
  }

*/

  totalMassChange    = gSum(theta*mesh.V().field()) - totalMass0;
 // totalWaterCompress += volWaterCompress;
  totalFlux          += fluxBoundary*runTime.deltaT().value();
  massBalanceRatio    = totalMassChange/totalFlux;
/*  
  //current step
  Info << "\nMass balance (current step): " << endl;
  Info << "       volume change         = " << volMassChange << endl;
  Info << "       water compress        = " << volWaterCompress << endl;
  Info << "       flux through boundary = " << fluxBoundary  << endl;
  Info << "       absolute error        = " << (volMassChange+volWaterCompress+fluxBoundary) << endl;

  //cumulative
  Info << "\nMass balance (cumulative): " << endl;
  Info << "       volume change         = " << totalMassChange_in_domain << endl;
  Info << "       water compress        = " << totalWaterCompress_in_domain << endl;
  Info << "       flux through boundary = " << totalFlux_through_boundary  << endl;
  Info << "       absolute error        = " 
       << (totalMassChange_in_domain+totalWaterCompress_in_domain+totalFlux_through_boundary) << endl;
  Info << "       relative error (\%)   = " 
       << 100*(1.0-mag(totalFlux_through_boundary/(totalMassChange_in_domain+totalWaterCompress_in_domain+SMALL))) 
       << endl;
  Info << endl;
   
*/
// ************************************************************************* //
