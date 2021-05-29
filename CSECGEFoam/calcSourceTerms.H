//alpha cutoff value for source term distribution:(no source terms in cells with cutoff <alpha1<1-cutoff)
scalar cutoff = 1.0e-3;

//- Calculation of source terms
Info << "Calculation of source terms" << endl;

volScalarField psi0Tild = mag(fvc::grad(alpha1));

dimensionedScalar intPsi0Tild = fvc::domainIntegrate(psi0Tild);
dimensionedScalar intAlphaPsi0Tild = fvc::domainIntegrate(alpha1*psi0Tild);

dimensionedScalar N ("N", dimensionSet(0,0,0,0,0,0,0), 2.0);
if (intAlphaPsi0Tild.value() > 1e-99)
{
	N = intPsi0Tild/intAlphaPsi0Tild;
}

volScalarField jMassT = Dl*mag(fvc::grad(Ci))*molarmass;

psi0 = N*jMassT*alpha1*psi0Tild; //Calcualting mass transfer through the interface

dimensionedScalar intPsi0 = fvc::domainIntegrate(psi0);

//- Smearing of source term field
fvScalarMatrix psiEqn
(
	fvm::Sp(scalar(1.0),psi) - fvm::laplacian(DPsi,psi) == psi0
);

psiEqn.solve();

//- Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining source term field
dimensionedScalar intPsiLiquid ("intPsiLiquid", dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar intPsiGas ("intPsiLiquid", dimensionSet(1,0,-1,0,0,0,0), 0.0);

forAll(mesh.C(),iCell)
{
	if (alpha1[iCell] < cutoff)
	{
		intPsiGas.value() += (1.0-alpha1[iCell])*psi[iCell]*mesh.V()[iCell];
	}
	else if (alpha1[iCell] > 1.0-cutoff)
	{
		intPsiLiquid.value() += alpha1[iCell]*psi[iCell]*mesh.V()[iCell];
	}
}

//- Calculate Nl and Ng
dimensionedScalar Nl ("Nl", dimensionSet(0,0,0,0,0,0,0), 2.0);
dimensionedScalar Ng ("Ng", dimensionSet(0,0,0,0,0,0,0), 2.0);

reduce(intPsiLiquid.value(),sumOp<scalar>());
reduce(intPsiGas.value(),sumOp<scalar>());

if (intPsiLiquid.value() > 1e-99)
{
    Nl = intPsi0/intPsiLiquid;
}
if (intPsiGas.value() > 1e-99)
{
    Ng = intPsi0/intPsiGas;
}

// Source term contiuity eqn - p.Eqn        
// Source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
forAll(mesh.C(),iCell)
{
	if (alpha1[iCell] < cutoff)
	{
		rhoSource[iCell] = Ng.value()*(1.0-alpha1[iCell])*psi[iCell];
	}
	else if (alpha1[iCell] > 1.0-cutoff)
	{
		rhoSource[iCell] = 0.0;
	}
	else
	{
		rhoSource[iCell] = 0.0;
	}
}


//- Source term in species concentration equation -Y.H
hCSource = jMassT*psi0Tild*(-N*alpha1);


// Volume fraction equation - alpha.Eqn
        alphaSourceSu = alpha1*fvc::div(phi);//can also be written as alphaSourceSu = alpha1*rhoSource/rho (see Appendix B in PhD thesis); 
	// alphaSourceSu regardless of how its calcualted goes into the prediction step of semi implicit MULES
        alphaSourceSp = alphaSourceSu*0.0;
