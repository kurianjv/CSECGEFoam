Info<<"Volume integral of bubble: " << fvc::domainIntegrate(alpha2) << endl;
Info<<"Volume integral of water: " << fvc::domainIntegrate(alpha1) << endl;

 Info<<"Total moles of Ci in the system: " << fvc::domainIntegrate(alpha1*Ci) + fvc::domainIntegrate(alpha2*Ci) << endl;
 Info<<"Total moles of Ci in bubble: " << fvc::domainIntegrate(alpha2*Ci) << endl;
 Info<<"Total moles of Ci in liquid: " << fvc::domainIntegrate(alpha1*Ci) << endl;

Info  << "   Uavg (m/s) = " << mag(average(U)).value() << endl;
Info  << "   Umax (m/s) = " << max(mag(U)).value() << endl; 

Info  << "   p_rgh in bubble (kg/ms2) = " << fvc::domainIntegrate(alpha2*p_rgh)/fvc::domainIntegrate(alpha2) << endl; 
Info  << "   p in bubble (kg/ms2) = " << fvc::domainIntegrate(alpha2*p)/fvc::domainIntegrate(alpha2) << endl; 

Info << "Total water flow at out of domain"   <<  fvc::domainIntegrate(fvc::div(phi)) << endl;


volScalarField Ux = U.component(0);
volScalarField Uy = U.component(1);
volScalarField Uz = U.component(2);

Info  << "   Bubble avg. Ux (m/s) = " << fvc::domainIntegrate(alpha2*Ux)/fvc::domainIntegrate(alpha2) << endl;
Info  << "   Bubble avg. Uy (m/s) = " << fvc::domainIntegrate(alpha2*Uy)/fvc::domainIntegrate(alpha2) << endl;
Info  << "   Bubble avg. Uz (m/s) = " << fvc::domainIntegrate(alpha2*Uz)/fvc::domainIntegrate(alpha2) << endl;
Info  << "   Bubble avg. U (m/s) = " << fvc::domainIntegrate( alpha2*mag(U) )/fvc::domainIntegrate(alpha2) << endl;

//calcualting fluxes through boundaries
label Patch1 = mesh.boundaryMesh().findPatchID("right_wall");
label Patch2 = mesh.boundaryMesh().findPatchID("top_wall");
label Patch3 = mesh.boundaryMesh().findPatchID("bottom_wall"); 
label Patch4 = mesh.boundaryMesh().findPatchID("left_wall");

Info << "   Fraction of liquid on right_wall "  << gSum(alpha1.boundaryField()[Patch1]*mesh.magSf().boundaryField()[Patch1])/gSum(mesh.magSf().boundaryField()[Patch1])  << endl; 
Info << "   Fraction of liquid on top_wall "  <<  gSum(alpha1.boundaryField()[Patch2]*mesh.magSf().boundaryField()[Patch2])/gSum(mesh.magSf().boundaryField()[Patch2])  << endl;
Info << "   Fraction of liquid on bottom_wall "  <<  gSum(alpha1.boundaryField()[Patch3]*mesh.magSf().boundaryField()[Patch3])/gSum(mesh.magSf().boundaryField()[Patch3])  << endl;
Info << "   Fraction of liquid on left_wall "  <<  gSum(alpha1.boundaryField()[Patch4]*mesh.magSf().boundaryField()[Patch4])/gSum(mesh.magSf().boundaryField()[Patch4])  << endl;

Info << "   Avg potential at right_wall "  << gSum(EP.boundaryField()[Patch1]*mesh.magSf().boundaryField()[Patch1])/gSum(mesh.magSf().boundaryField()[Patch1])  << endl; 
Info << "   Avg potential at top_wall "  <<  gSum(EP.boundaryField()[Patch2]*mesh.magSf().boundaryField()[Patch2])/gSum(mesh.magSf().boundaryField()[Patch2])  << endl;
Info << "   Avg potential at bottom_wall "  <<  gSum(EP.boundaryField()[Patch3]*mesh.magSf().boundaryField()[Patch3])/gSum(mesh.magSf().boundaryField()[Patch3])  << endl;
Info << "   Avg potential at left_wall "  <<  gSum(EP.boundaryField()[Patch4]*mesh.magSf().boundaryField()[Patch4])/gSum(mesh.magSf().boundaryField()[Patch4])  << endl;

Info << "   Avg current on right_wall "  <<  gSum( currentdensity.boundaryField()[Patch1] & mesh.Sf().boundaryField()[Patch1] ) << endl;
Info << "   Avg current on top_wall "    <<  gSum( currentdensity.boundaryField()[Patch2] & mesh.Sf().boundaryField()[Patch2] )  << endl;
Info << "   Avg current on bottom_wall " <<  gSum( currentdensity.boundaryField()[Patch3] & mesh.Sf().boundaryField()[Patch3] )  << endl;
Info << "   Avg current on left_wall "   <<  gSum( currentdensity.boundaryField()[Patch4] & mesh.Sf().boundaryField()[Patch4] ) << endl;

Info << "   Water flow at right_wall "  <<  gSum(phi.boundaryField()[Patch1])  << endl; //m3/s
Info << "   Water flow at top_wall "    <<  gSum(phi.boundaryField()[Patch2])  << endl;
Info << "   Water flow at bottom_wall " <<  gSum(phi.boundaryField()[Patch3])  << endl;
Info << "   Water flow at left_wall "   <<  gSum(phi.boundaryField()[Patch4])  << endl;

Info << "   Dissolved gas flow at right_wall " <<  gSum(Ci.boundaryField()[Patch1] * phi.boundaryField()[Patch1])  << endl;//mol/s
Info << "   Dissolved gas flow at top_wall "   <<  gSum(Ci.boundaryField()[Patch2] * phi.boundaryField()[Patch2])  << endl;
Info << "   Dissolved gas flow at bottom_wall "<<  gSum(Ci.boundaryField()[Patch3] * phi.boundaryField()[Patch3])  << endl;
Info << "   Dissolved gas flow at left_wall "  <<  gSum(Ci.boundaryField()[Patch4] * phi.boundaryField()[Patch4])  << endl;


