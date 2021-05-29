{
Info << " Algebraic averaging of conductivities" << endl;
	volScalarField k_avg("k_avg", (alpha1*k_liquid)+(k_gas*alpha2) );//k_liquid*k_gas/((alpha2*k_liquid)+(k_gas*alpha1)) );

        fvScalarMatrix EPEqn
        (
           fvm::laplacian(fvc::interpolate(k_avg) ,EP)
        );

        EPEqn.solve();


Info << " Min(EP) ="<< min(EP).value() << " Max(EP) =" << max(EP).value() << endl;

        currentdensity = -k_avg*fvc::grad(EP);

				
}
