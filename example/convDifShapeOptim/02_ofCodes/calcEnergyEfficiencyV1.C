/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cuttingPlane.H"

int main(int argc, char *argv[])
{
	argList::addNote
    (
        "Input arguments:\n"
        "----------------\n"
        "  xInl - x-coordinate of the inlet\n"
        "  xSuc - x-coordinate of the suction\n"
        "  xOut - x-coordinate of the outlet \n"
        "  rho - fluid mass density\n"
        "  hCut - horizontal cut between inlet and suction\n"
    );

    // prepare argument list
    argList::noParallel();
    argList::validArgs.append("xInl");
    argList::validArgs.append("xSuc");
    argList::validArgs.append("xOut");
    argList::validArgs.append("rho");
    argList::validArgs.append("hCut");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

	// read input arguments
	const scalar xInl = args.argRead<scalar>(1);
	const scalar xSuc = args.argRead<scalar>(2);
	const scalar xOut = args.argRead<scalar>(3);
	const scalar rho = args.argRead<scalar>(4);
	const scalar hCut = args.argRead<scalar>(5);

	// create the inlet cutting plane
	const point inlPnt(xInl,0.0,0.0);
	const vector inlNorm(1.0,0.0,0.0);

	const plane inlPlane(inlPnt,inlNorm);
	const cuttingPlane inlCutPlane(inlPlane, mesh, 1);

	// create the suction cutting plane
	const point sucPnt(xSuc,0.0,0.0);
	const vector sucNorm(1.0,0.0,0.0);

	const plane sucPlane(sucPnt,sucNorm);
	const cuttingPlane sucCutPlane(sucPlane, mesh, 1);

	// create the outlet cutting plane
	const point outPnt(xOut,0.0,0.0);
	const vector outNorm(1.0,0.0,0.0);

	const plane outPlane(outPnt,outNorm);
	const cuttingPlane outCutPlane(outPlane, mesh, 1);

    // load inlet patch
    word patchName = "inlet";
    const label patchIDInl = mesh.boundaryMesh().findPatchID(patchName);

    const labelList& inletCells = mesh.boundary()[patchIDInl].faceCells();

    // load outlet patch
    patchName = "suction";
    const label patchIDSuc = mesh.boundaryMesh().findPatchID(patchName);

    const labelList& suctionCells = mesh.boundary()[patchIDSuc].faceCells();

    // read the U field
    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // read the p field
    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

	// compute average pressure at the inlet
	scalar pInl = 0;
	scalar VInl = 0;

	forAll(inlCutPlane.cutCells(), i)
	{
        scalar yC = mesh.C()[inlCutPlane.cutCells()[i]].y();
        if (yC < hCut)
        {
		    pInl += p[inlCutPlane.cutCells()[i]]*mesh.V()[inlCutPlane.cutCells()[i]];
		    VInl += mesh.V()[inlCutPlane.cutCells()[i]];
        }
	}

	pInl /= VInl/rho;

	// compute average pressure at the suction
	scalar pSuc = 0;
	scalar VSuc = 0;

	forAll(sucCutPlane.cutCells(), i)
	{
        scalar yC = mesh.C()[sucCutPlane.cutCells()[i]].y();
        if (yC > hCut)
        {
		    pSuc += p[sucCutPlane.cutCells()[i]]*mesh.V()[sucCutPlane.cutCells()[i]];
		    VSuc += mesh.V()[sucCutPlane.cutCells()[i]];
        }
	}

	pSuc /= VSuc/rho;

	// compute average pressure on outlet
	scalar pOut = 0;
	scalar VOut = 0;

	forAll(outCutPlane.cutCells(), i)
	{
	    pOut += p[outCutPlane.cutCells()[i]]*mesh.V()[outCutPlane.cutCells()[i]];
	    VOut += mesh.V()[outCutPlane.cutCells()[i]];
	}

	pOut /= VOut/rho;

    // compute flow rate over the inlet
    scalar QInl = 0;
    forAll(inletCells, i)
    {
        QInl += U[inletCells[i]] & mesh.boundary()[patchIDInl].Sf()[i];
    }

    // compute flor rate over the suction
    scalar QSuc = 0;
    forAll(suctionCells, i)
    {
        QSuc += U[suctionCells[i]] & mesh.boundary()[patchIDSuc].Sf()[i];
    }

	// compute the energy efficiency
    scalar Eeff = QSuc/QInl*(pOut - pSuc)/(pInl - pOut);

    Info << "Flow rate at the inlet is QInl = " << QInl << endl;
    Info << "Flow rate at the suction is QSuc = " << QSuc << endl;
    Info << "Pressure at the inlet is pInl = " << pInl << endl;
    Info << "Pressure at the suction is pSuc = " << pSuc << endl;
    Info << "Pressure at the outlet is pOut = " << pOut << endl;

	Info << "Energy effiency of the ejector" << endl
		<< "is Eeff = " << Eeff << endl;

	Info << endl;
	Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
