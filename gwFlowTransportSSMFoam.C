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
    groundwaterFoam

Description
    Transient solver for Richards equation. 
    A Picard loop is used for linearization.
    Permeability is isotropic (K == volScalarField)

Developers
    P. Horgue, J. Franc, R. Guibert and G. Debenest
    
Modified by
    E. Huber

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fvOptions.H" // add Manu
#include "harmonic.H"
#include "incompressiblePhase.H"
#include "capillarityModel.H"
#include "relativePermeabilityModel.H"
#include "OFstream.H" // add Manu

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createthetaFields.H"
    #include "readPicardControls.H"
    //#include "createFvOptions.H" // add MANU
	  #include "createWellbores.H" // add Manu
	  #include "createSourceC.H" // add Manu

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //----------- ADD MANU	
    scalar i = 0;
    // List<Type> values(cells_.size(), injectionRate_[fieldI]);
    labelHashSet selectedCells;
    labelList cells_;
    //forAll(points_, i)
    for(i = 0; i < position_.size(); i++)
    {
        label celli = mesh.findCell(position_[i]);
        if (celli >= 0)
        {
            selectedCells.insert(celli);
            Info<< "Location " << position_[i] << " is at the cell #"
          <<  celli << " + val = " << vals_[i] << "." << nl << endl;
        }
        //vals[i] = 1;
        label globalCelli = returnReduce(celli, maxOp<label>());
        if (globalCelli < 0)
        {
            WarningInFunction
                << "Unable to find owner cell for point " << position_[i]
                << endl;
        }
    }

    cells_ = selectedCells.toc();
    
  /*
    volScalarField h_ss
    (
        IOobject
        (
            "h_ss",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
            (
              "0.0", 
              dimensionSet(0,0,-1,0,0,0,0), 
              0.0 
            ),
        zeroGradientFvPatchScalarField::typeName
    );
  */
    simuLogFile << "runTime" << tab;
    simuLogFile << "resh" << tab;
    simuLogFile << "resTheta" << tab;
    simuLogFile << "resC" << tab;
    simuLogFile << "totMassChange" << tab;
   // simuLogFile << "totWatCompres" << tab;
    simuLogFile << "totFlux" << tab;
    simuLogFile << "ratio" << endl;
    //----------- END ADD MANU
	
    Info<< "\nStarting time loop\n" << endl;
    label iterPicard=0;

    scalar totalMass0         = gSum(theta*mesh.V().field());
    scalar totalFlux          = 0.0;
    scalar totalMassChange    = 0.0;
    //scalar totalWaterCompress = 0.0;
    scalar massBalanceRatio   = 1;
  
    scalar resh = tolh + 1;						// ADD MANU
    scalar resC = tolC + 1;           // ADD MANU
    scalar resTheta = 0;           // ADD MANU
    scalar currentMaxC = 0;           // ADD MANU
    scalar currentMinC = 0;           // ADD MANU
    volScalarField hprevIter = h;
    volScalarField thetaprevIter = theta;
   // while (runTime.run())					  // COMMENTED BY MANU
    while ((resC > tolC || resh > tolh || !runTime.outputTime()) && runTime.run())	// ADD MANU > for steady state solution
    {
       if(resh > tolh)
       {
      
          #include "setDeltaT.H"
          runTime++;
      
          Info << "Time = " << runTime.timeName() << nl << endl;
         
          hprevIter = h;
          thetaprevIter = theta;
          scalar resPicard = 1;
          iterPicard = 0;
          while (resPicard > tolPicard)
          {
              #include "hEqn.H"
              #include "updateProperties.H"
              iterPicard++;
              if (iterPicard >= (2*maxIterPicard))
              {
                  Warning() <<  " Max iteration reached in Picard loop" << endl;
                  break;
              }
          }
          
         
          
          // Calculate mass balance
          #include "calcMassBalance.H"
         
          //Info << "Phi" << " Min(phi) = " << gMin(mag(phi)->internalField()) << " Max(phi) = " << gMax(mag(phi)->internalField()) <<  endl; 
          Info << "Saturation theta " << " Min(theta) = " << gMin(theta) << " Max(theta) = " << gMax(theta) <<  endl;
          Info << "Head pressure h  " << " Min(h) = " << gMin(h) << " Max(h) = " << gMax(h) <<  endl;

          resh = gMax((mag(h - hprevIter))->internalField());
          resTheta = gMax((mag(theta - thetaprevIter))->internalField());
          //resh = gMax((mag(h-h0))->internalField());
          Info << ">>> Pressure residual (h) = " << resh <<  " > " << tolh << endl;  
          Info << ">>> Vol. water content residual (theta) = " << resTheta << endl;  
       }else
       {
          #include "CourantNo.H"        	// ADD MANU
		      #include "setDeltaTCo.H"        // ADD MANU  
          
          runTime++;
        
          Info << "Time = " << runTime.timeName() << nl << endl;
         
       }
      
        C.storePrevIter();

        fvScalarMatrix CEqn
        (
          fvm::ddt(theta, C) 
          + fvm::div(phi, C)
          - fvm::laplacian( theta * DC, C, "laplacian(DC,C)")
          ==
           fvm::Sp(-SrcExt * Wext, C)
        );

        //fvOptions.constrain(CEqn);
        CEqn.setValues(cells_, vals_);
        CEqn.relax();
        CEqn.solve();
        //fvOptions.correct(C);
      
        currentMaxC = gMax(C);
        currentMinC = gMin(C);
        Info << "Concentration C " << " Min(C) = " << currentMinC << " Max(C) = " << currentMaxC << endl;

        resC = gMax((mag(C-C.prevIter()))->internalField());
        Info << ">>> Concentration residual (C) = " << resC <<  " > " << tolC << endl; 
      
      
        //output some information to the log file at current time
        simuLogFile << runTime.value() << tab;
        simuLogFile << resh << tab;
        simuLogFile << resTheta << tab;
        simuLogFile << resC << tab;
        simuLogFile << totalMassChange << tab;
        //simuLogFile << totalWaterCompress_in_domain << tab;
        simuLogFile << totalFlux << tab;
        simuLogFile << massBalanceRatio << endl;
      
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        if (currentMaxC > maxC)
        {
          C.write();	
          h.write();
          Utheta.write();
          phi.write();
          K.write();
          FatalError
            << "Too large concentration ( C > maxC = " << maxC << ")" << nl
            << exit(FatalError);
        }
    }
	
    C.write();	
    h.write();
    Utheta.write();
    phi.write();
    K.write();
    //----- END ADD MANU

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
