/** Author	: Ge Zhnag
 *  Email	:
 *	Date	:	*/

#include "StructureFactor.h"

#include <complex>
#include <cmath>
#include <algorithm>

/** \file StructureFactor.cpp
 *	\brief	Function implementations for structure factor. */

/** \struct A class to store data in each bin. */
struct SkBin
{
	double Sum1,	//< count
		SumK1,		//<	sum of the magnitudes of wavevectors (|k|)
		SumK2,		//< sum of squared magnitudes of wavevectors (|k|^2)
		SumS,		//< sum of spectral functions
		SumS2;		//< sum of squared spectral functions
	SkBin() : Sum1(0), SumK1(0), SumK2(0), SumS(0), SumS2(0)
	{}
};


double StructureFactor(const Configuration & Config, const GeometryVector & k)
{
	std::complex<double> rho;
	size_t NumParticle = Config.NumParticle();
	for(size_t i=0; i<NumParticle; i++)
	{
		rho+=std::exp(std::complex<double>(0, k.Dot(Config.GetCartesianCoordinates(i))));
	}
	return (rho.real()*rho.real()+rho.imag()*rho.imag())/NumParticle;
}

double WeightedStructureFactor(const PeriodicCellList<std::complex<double>> & Config, const GeometryVector & k) {
	std::complex<double> rho;
	size_t NumParticle = Config.NumParticle();
	for (size_t i = 0; i<NumParticle; i++)
	{
		rho += Config.GetCharacteristics(i)*std::exp(std::complex<double>(0, k.Dot(Config.GetCartesianCoordinates(i))));
	}
	return (rho.real()*rho.real() + rho.imag()*rho.imag()) / NumParticle;
}




void IsotropicStructureFactor(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision, double SampleProbability, size_t option, double CircularKMin)
{
	Results.clear();
	if(!(CircularKMax > 0))
		return;
	if(!(LinearKMax > 0))
		return;
	if(CircularKMin > CircularKMax)
		return;

	if (KPrecision == 0.0)
	{
		std::cerr << "Warning in IsotropicStructureFactor : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
		Configuration c = GetConfigsFunction(0);
		KPrecision = std::sqrt(c.GetReciprocalBasisVector(0).Modulus2());
	}
	size_t NumBin = std::floor(LinearKMax / KPrecision) + 1;
	std::vector<SkBin> vSkBin; //(NumBin, SkBin());
	DimensionType d = 0;
	double V = 0.0;
	GeometryVector prevBasis [ ::MaxDimension];
	std::vector<GeometryVector> ks;
	std::vector<double> k_amp;
	if(Verbosity>1){
		std::cout<<"Computing S(k)";
		if (option == 0){
			std::cout << ": averaged over reciprocal lattice points\n";
			vSkBin.resize(NumBin, SkBin());
		}
		else if (option == 1){
			std::cout << ": averaged over volume in the reciprocal space\n";
			vSkBin.resize(NumBin, SkBin());
		}
		else if (option == 2){
			std::cout << ": structure factor at Bragg peaks\n";
		}
		else{
			std::cout << ": Wrong option!!\n";
			return ;
		}
	}
	progress_display pd(NumConfigs);
	for(size_t j=0; j<NumConfigs; j++)
	{
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration CurrentConfig = GetConfigsFunction(j);
		if(CurrentConfig.GetDimension()==0)
			break;
		if(j > 0)
		{
			bool SameBasis = true;
			for(DimensionType i=0; i<CurrentConfig.GetDimension(); i++)
				if( !(prevBasis[i]==CurrentConfig.GetBasisVector(i)) )
					SameBasis=false;

			if(SameBasis==false){
				ks=GetKs(CurrentConfig, CircularKMax, LinearKMax, SampleProbability, true, CircularKMin);
				k_amp.resize(ks.size(), 0.0);
				for (size_t i = 0; i < ks.size(); i++){
					k_amp[i] = std::sqrt(ks[i].Modulus2());
				}
			}
		}
		else{
			ks = GetKs(CurrentConfig, CircularKMax, LinearKMax, SampleProbability, false, CircularKMin);
			k_amp.resize(ks.size(), 0.0);
			for (size_t i = 0; i < ks.size(); i++){
				k_amp[i] = std::sqrt(ks[i].Modulus2());
			}
			d = CurrentConfig.GetDimension();
			V = CurrentConfig.PeriodicVolume();
			if (option == 2){
				vSkBin.resize(ks.size(), SkBin());
			}
		}

		const size_t end = ks.size();

		if (option == 0 || option == 1){
	#pragma omp parallel for schedule(guided)
			for(size_t i=0; i<end; i++)
			{
				double s=StructureFactor(CurrentConfig, ks[i]);
				//double k2 = ks[i].Modulus2();
				size_t Bin = std::floor(k_amp[i] / KPrecision);
		#pragma omp atomic
				vSkBin[Bin].Sum1 += 1.0;
				//-----added-----
		#pragma omp atomic
				vSkBin[Bin].SumK1 += k_amp[i]; 
				//---------------
		#pragma omp atomic
				vSkBin[Bin].SumK2 += k_amp[i]*k_amp[i];
		#pragma omp atomic
				vSkBin[Bin].SumS += s;
		#pragma omp atomic
				vSkBin[Bin].SumS2 += s*s;
					
			}
		}
		else 
		{
			/* at each wavevector */
			#pragma omp parallel for schedule(guided)
			for(size_t i=0; i<end; i++)
			{
				double s=StructureFactor(CurrentConfig, ks[i]);
				//double k2 = ks[i].Modulus2();
				vSkBin[i].Sum1 += 1.0;
				vSkBin[i].SumK1 += k_amp[i];
				vSkBin[i].SumK2 += k_amp[i]*k_amp[i];
				vSkBin[i].SumS += s;
				vSkBin[i].SumS2 += s*s;				
			}			
		}

		for(DimensionType i=0; i<CurrentConfig.GetDimension(); i++)
			prevBasis[i]=CurrentConfig.GetBasisVector(i);
		pd++;
	}

	double Omega = std::pow(2.*pi, d)/V/NumConfigs;
	if  (option == 0 || option == 1){
		// when wavenumbers are binned.
		Results.reserve(NumBin);
	}
	else {
		//Results.resize(vSkBin.size(), GeometryVector(int(4)));
		Results.reserve(vSkBin.size());
	}

	//#pragma omp parallel for schedule(guided)
	for (size_t i = 0; i < vSkBin.size(); i++ ){
		auto iter = vSkBin.begin() + i;
		if (iter->Sum1 > 0.0){
			GeometryVector temp(4);
			//temp.x[0] = std::sqrt(iter->SumK2 / iter->Sum1);
			//temp.x[1] = iter->SumS / iter->Sum1;
			//temp.x[2] = KPrecision;
			//temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1)); // I modified it
			temp.x[0] = iter->SumK1 / iter->Sum1;
			double var = (iter->SumK2 / (iter->Sum1) - temp.x[0] * temp.x[0]) / (iter->Sum1);
			temp.x[2] = (var > 0) ? std::sqrt(var) : 0;
			
			//CircularKMax
			double k = temp.x[0];//KPrecision*(double)(iter - vSkBin.begin());// + 0.5);
			if (option == 0 || (k > CircularKMax)){
				/* average over number of wavevectors */
				temp.x[1] = iter->SumS / iter->Sum1;
				temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1)); // I modified it

				Results.push_back(temp);
			}
			else if (option == 1){
				/* average over volume in reciprocal space*/
				double vd = 0.5*(HyperSphere_Volume(d, k+0.5*KPrecision) - HyperSphere_Volume(d, k-0.5*KPrecision)); //HyperSphere_SurfaceArea(d, k)*KPrecision/2.0;
				temp.x[1] = iter->SumS * Omega / vd;
				temp.x[2] = 0.5*KPrecision; //Omega / vd;
				temp.x[3] = std::sqrt((iter->SumS2  - temp.x[1] * temp.x[1]) / (iter->Sum1) ) * Omega / vd * iter->Sum1;

				Results.push_back(temp);
			}
			else if (option == 2){
				/* Bragg peaks */
				temp.x[1] = iter->SumS / iter->Sum1;
				temp.x[3] = iter->Sum1; //std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1));

				Results.push_back(temp);
			}
			
		}
	}

	if (option == 2){
		std::sort(Results.begin(), Results.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.x[0] < right.x[0];} );
	}

	if(Verbosity>2)
		std::cout<<"done!\n";
}

void IsotropicStructureFactor_weighted(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision, double SampleProbability)
{
	Results.clear();
	if (!(CircularKMax>0))
		return;
	if (!(LinearKMax>0))
		return;

	if (KPrecision == 0.0)
	{
		std::cerr << "Warning in IsotropicStructureFactor : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
		PeriodicCellList<std::complex<double>> c = GetConfigsFunction(0);
		KPrecision = std::sqrt(c.GetReciprocalBasisVector(0).Modulus2());
	}
	size_t NumBin = std::floor(LinearKMax / KPrecision) + 1;
	std::vector<SkBin> vSkBin(NumBin, SkBin());


	GeometryVector prevBasis[::MaxDimension];
	std::vector<GeometryVector> ks;
	if (Verbosity>1)
		std::cout << "Computing S(k)";
	progress_display pd(NumConfigs);



	for (size_t j = 0; j<NumConfigs; j++)
	{
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		PeriodicCellList<std::complex<double>> CurrentConfig = GetConfigsFunction(j);
		if (CurrentConfig.GetDimension() == 0)
			break;
		if (j != 0)
		{
			bool SameBasis = true;
			for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
				if (!(prevBasis[i] == CurrentConfig.GetBasisVector(i)))
					SameBasis = false;

			if (SameBasis == false)
				ks = GetKs(CurrentConfig, CircularKMax, LinearKMax, SampleProbability);
		}
		else
			ks = GetKs(CurrentConfig, CircularKMax, LinearKMax, SampleProbability);

		signed long end = ks.size();

#pragma omp parallel for schedule(guided)
		for (signed long i = 0; i<end; i++)
		{
			double s = WeightedStructureFactor(CurrentConfig, ks[i]);
			double k2 = ks[i].Modulus2();
			size_t Bin = std::floor(std::sqrt(k2) / KPrecision);
#pragma omp atomic
			vSkBin[Bin].Sum1 += 1.0;
			//-----added-----
#pragma omp atomic
			vSkBin[Bin].SumK1 += sqrt(k2);
			//---------------
#pragma omp atomic
			vSkBin[Bin].SumK2 += k2;
#pragma omp atomic
			vSkBin[Bin].SumS += s;
#pragma omp atomic
			vSkBin[Bin].SumS2 += s * s;

		}
		for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
			prevBasis[i] = CurrentConfig.GetBasisVector(i);
		pd++;
	}


	for (auto iter = vSkBin.begin(); iter != vSkBin.end(); iter++)
	{
		if (iter->Sum1 != 0.0)
		{
			GeometryVector temp(4);
			//temp.x[0] = std::sqrt(iter->SumK2 / iter->Sum1);
			//temp.x[1] = iter->SumS / iter->Sum1;
			//temp.x[2] = KPrecision;
			//temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1)); // I modified it
			temp.x[0] = iter->SumK1 / iter->Sum1;
			temp.x[1] = iter->SumS / iter->Sum1;
			double var = (iter->SumK2 / (iter->Sum1) - temp.x[0] * temp.x[0]) / (iter->Sum1);
			temp.x[2] = (var > 0) ? std::sqrt(var) : 0;
			temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1)); // I modified it
			Results.push_back(temp);
		}
	}

	if (Verbosity>2)
		std::cout << "done!\n";
}

