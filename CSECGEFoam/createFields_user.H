IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

/************************************ Species transport CCST **************************************/
Info<< "Reading field Ci\n" << endl;
volScalarField Ci
(
    IOobject
    (
        "Ci",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

	dictionary phase1 = transportProperties.subDict("water");
	dictionary phase2 = transportProperties.subDict("air");
	dimensionedScalar Dl = phase1.lookup("D");
	dimensionedScalar Dg = phase2.lookup("D");
        dimensionedScalar Hei (transportProperties.lookup("Hei") );
////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************** SSF *************************************************************/ 
    volVectorField fc
    (
        IOobject
        (
            "fc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("fc", dimMass/(dimLength*dimLength*dimTime*dimTime), vector(0,0,0))
    );
///////////////////////////////////////////////////////////////////////////////////////////////////

/******************************** Source terms ***************************************************/
    Info<< "Reading field psi0\n" << endl;
    volScalarField psi0
    (
        IOobject
        (
            "psi0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        psi0
    );

	volScalarField rhoSource
	(
		IOobject
		(
			"rhoSource",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("rhoSource", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
	);

	volScalarField hCSource
	(
		IOobject
		(
			"hCSource",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("hCSource", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
	);

	volScalarField alphaSourceSu
	(
		IOobject
		(
			"alphaSourceSu",
			runTime.timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("alphaSourceSu", dimensionSet(0,0,-1,0,0,0,0), 0.0)
	);

	volScalarField alphaSourceSp
	(
		IOobject
		(
			"alphaSourceSp",
			runTime.timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("alphaSourceSp", dimensionSet(0,0,-1,0,0,0,0), 0.0)
	);

	dimensionedScalar mGrowth ("mGrowth", dimensionSet(1,0,-1,0,0,0,0), 0.0);
	dimensionedScalar DPsi = transportProperties.lookup("DPsi");
        dimensionedScalar molarmass ( transportProperties.lookup("molarmass") );

//////////////////////////////////////////////////////////////////////////////////////////////////

Info<< "Reading field EP\n" << endl;
volScalarField EP
(
    IOobject
    (
        "EP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
        dimensionedScalar k_liquid ( transportProperties.lookup("k_liquid") );
        dimensionedScalar k_gas ( transportProperties.lookup("k_gas") );

	volVectorField currentdensity
	(
		IOobject
		(
			"currentdensity",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedVector("currentdensity", dimensionSet(0,-2,0,0,0,1,0), vector(0.0, 0.0, 0.0) )
	);

