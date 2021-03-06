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

Global
    readPicardControls

Description
    Read the control parameters used by setDeltaT for solver with Picard's iterations

\*---------------------------------------------------------------------------*/

const dictionary& PicardDict = mesh.solutionDict().subDict("Picard");

const scalar tolPicard = PicardDict.lookupOrDefault<scalar>("tolerance", 1.);

const label maxIterPicard = PicardDict.lookupOrDefault<label>("maxIter",10);

const label minIterPicard = PicardDict.lookupOrDefault<label>("minIter",3);

const label maxIterStability = PicardDict.lookupOrDefault<label>("maxIterStability",5);
label iterStability=0;

//--- start MOD MANU
const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

//--- end MOD MANU
const scalar dTFactIncrease
(
    readScalar(runTime.controlDict().lookup("dTFactIncrease"))
);
const scalar dTFactDecrease
(
    readScalar(runTime.controlDict().lookup("dTFactDecrease"))
);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
