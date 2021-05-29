{

volScalarField Di("Di", Dl*Dg/((alpha2*Dl)+(Dg*alpha1)) );
surfaceScalarField phiCi = fvc::interpolate(Di*(1-Hei)/(alpha1 + Hei*alpha2))*fvc::snGrad(alpha1)*mesh.magSf();

surfaceScalarField phicCi(mixture.cAlpha()*mag(phi/mesh.magSf()));
surfaceScalarField ur(phicCi*mixture.nHatf());

surfaceScalarField phiCir = fvc::interpolate((1-Hei)*alpha1*alpha2/(alpha1 + Hei*alpha2))*ur;

solve
	(
		fvm::ddt(Ci)
		+ fvm::div(phi,Ci)
		- fvm::laplacian(Di,Ci)
		+ fvm::div(phiCir,Ci)
		+ fvm::div(phiCi, Ci)
                - (hCSource/molarmass),
		mesh.solver("Ci")
	);

Info << " Min(Ci) ="<< min(Ci).value() << " Max(Ci) =" << max(Ci).value() << endl;

}


