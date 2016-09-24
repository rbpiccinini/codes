
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Output parcel properties to grid cells
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

if (runTime.write())
{
class0_m=   0.0*class0_m;
class0_mU=   0.0*class0_mU;
class0_N=   0.0*class0_N;
class0_T=   0.0*class0_T;
class0_U=   0.0*class0_U;

class1_m=   0.0*class1_m;
class1_mU=   0.0*class1_mU;
class1_N=   0.0*class1_N;
class1_T=   0.0*class1_T;
class1_U=   0.0*class1_U;

class2_m=   0.0*class2_m;
class2_mU=   0.0*class2_mU;
class2_N=   0.0*class2_N;
class2_T=   0.0*class2_T;
class2_U=   0.0*class2_U;

class3_m=   0.0*class3_m;
class3_mU=   0.0*class3_mU;
class3_N=   0.0*class3_N;
class3_T=   0.0*class3_T;
class3_U=   0.0*class3_U;

drop_smd=0.0*drop_smd;

spray::iterator mol(dieselSpray.begin());


for (mol = dieselSpray.begin(); mol != dieselSpray.end(); ++mol)
{
	label cellI=mol().cell();
	
	if (mol().injector()==scalar(0))
	{
		class0_m[cellI]+=mol().m();
		class0_mU[cellI]+=mol().m()*mol().U();
		class0_N[cellI]+=scalar(1);
	}

	else if (mol().injector()==scalar(1))
	{
		class1_m[cellI]+=mol().m();
		class1_mU[cellI]+=mol().m()*mol().U();
		class1_N[cellI]+=scalar(1);
	}
	else if (mol().injector()==scalar(2))
	{
		class2_m[cellI]+=mol().m();
		class2_mU[cellI]+=mol().m()*mol().U();
		class2_N[cellI]+=scalar(1);
	}
	else if (mol().injector()==scalar(3))
	{
		class3_m[cellI]+=mol().m();
		class3_mU[cellI]+=mol().m()*mol().U();
		class3_N[cellI]+=scalar(1);
	}

}


for (mol = dieselSpray.begin(); mol != dieselSpray.end(); ++mol)
{
	label cellI=mol().cell();

	if (mol().injector()==scalar(0) && class0_N[cellI] != 0)
	{
		class0_T[cellI]+=mol().m()*mol().T()/class0_m[cellI];
		class0_U[cellI]+=mol().U()/class0_N[cellI];
	}
	else if (mol().injector()==scalar(1) && class1_N[cellI] != 0)
	{
		class1_T[cellI]+=mol().m()*mol().T()/class1_m[cellI];
		class1_U[cellI]+=mol().U()/class1_N[cellI];
	}
	else if (mol().injector()==scalar(2) && class2_N[cellI] != 0)
	{
		class2_T[cellI]+=mol().m()*mol().T()/class2_m[cellI];
		class2_U[cellI]+=mol().U()/class2_N[cellI];
	}
	else if (mol().injector()==scalar(3) && class3_N[cellI] != 0)
	{
		class3_T[cellI]+=mol().m()*mol().T()/class3_m[cellI];
		class3_U[cellI]+=mol().U()/class3_N[cellI];
	}

}


forAll(mesh.C(), cellI)
{
	drop_smd[cellI]=dieselSpray.smdcell(cellI);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	class0_m.write();
	class0_U.write();
	class0_T.write();

	class1_m.write();
	class1_U.write();
	class1_T.write();

	class2_m.write();
	class2_U.write();
	class2_T.write();

	class3_m.write();
	class3_U.write();
	class3_T.write();
	
	drop_smd.write();
}


	label Nparcels = dieselSpray.size();
    reduce(Nparcels, sumOp<label>());

    Info<< "\nNumber of parcels in system.... | "
        << Nparcels << endl
        << "Injected liquid mass........... | "
        << 1e6*dieselSpray.injectedMass(runTime.value()) << " mg" << endl
        << "Liquid Mass in system.......... | "
        << 1e6*dieselSpray.liquidMass() << " mg" << endl
        << "SMD, Dmax...................... | "
        << dieselSpray.smd()*1e6 << " mu, "
        << dieselSpray.maxD()*1e6 << " mu"
        << endl;
