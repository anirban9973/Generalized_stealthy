#include "PairCorrelation.h"
#include <omp.h>

void WeightedIsotropicTwoPairCorrelation(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MaxDistance, std::vector<GeometryVector> & Result, double dR)
{
	if(!(MaxDistance>0.0))
		return;
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"computing radial two-point function";
	progress_display pd(NumConfigs);
	DimensionType dim = GetConfigsFunction(0).GetDimension();

	//prepare for the iteration
	Result.clear();
	for (double R = dR; R <= MaxDistance; R+=dR)
	{
		GeometryVector temp(4);
		temp.x[0]=R;
//		temp.x[2]=dR;
		Result.push_back(temp);
	}

	//discard empty bins at the end. Such bins are useful when we don't have an upper bound of data (e.g. nearest neighbor distance binning)
	//here, they are useless
	//the bin counts
	double TotalAtomCount = 0;
	double AverageDensity = 0.0;
	for(size_t j=0; j<NumConfigs; j++)
	{
		PeriodicCellList<std::complex<double>> Config = GetConfigsFunction(j);
		Config.SetCellSize(0.35*MaxDistance);
		Config.RefineCellList();
		size_t NumParticle=Config.NumParticle();
		//pre-allocate
		Config.PrepareIterateThroughNeighbors(MaxDistance);
		std::vector<double> real(Result.size(),0.0);
		std::vector<double> imag(Result.size(),0.0);

#pragma omp parallel
		{
#pragma omp for schedule(guided)
			for(signed long i=0; i<NumParticle; i++)
			{
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance, 
					[&real, &imag, &i, &dR, &Config](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
				{
					if(SourceAtom>i || (SourceAtom==i && LatticeShift.Modulus2()!=0) )
					{
						double dist=std::sqrt(shift.Modulus2());
						size_t idx = (size_t)floor(dist/dR);
						if(idx<real.size())
						{
							std::complex<double> val = std::conj(Config.GetCharacteristics(i)) * Config.GetCharacteristics(SourceAtom);
						#pragma omp atomic
							real[idx]+= val.real();
						#pragma omp atomic
							imag[idx]+= val.imag();
						}
					}
				});
			}
		}
		// StartParticle=0;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
		pd++;

		for(size_t i=0; i<real.size(); i++)
		{
			double R = dR*(i+1);
			double factor = (HyperSphere_Volume(dim, R)-HyperSphere_Volume(dim, R-dR))*0.5 * Config.NumParticle();
			double val_real = real[i] / factor;
			double val_imag = imag[i] / factor;
			Result[i].x[1]+= val_real;  //values[i];
			Result[i].x[2]+= val_imag;  //values[i];
			Result[i].x[3]+= val_real*val_real + val_imag*val_imag; //values[i]*values[i];
		}
	}
	AverageDensity/=NumConfigs;

	for(int i = 0; i < Result.size(); i ++ ){
		double R = Result[i].x[0];
		Result[i].x[1] /= NumConfigs;
		Result[i].x[2] /= NumConfigs;
		Result[i].x[3] = std::sqrt( std::abs(Result[i].x[3] / NumConfigs - Result[i].x[1] * Result[i].x[1]));
		Result[i].x[0] -= 0.5*dR;
	}

	if(Verbosity>2)
		std::cout<<"done!\n";

}

void WeightedIsotropicSingleSidedTwoPairCorrelation(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, double MaxDistance, std::vector<GeometryVector> &Result, double dR)
{
	if(!(MaxDistance>0.0))
		return;
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"computing radial two-point function";
	progress_display pd(NumConfigs);
	DimensionType dim = GetConfigsFunction(0).GetDimension();

	//prepare for the iteration
	Result.clear();
	for (double R = dR; R <= MaxDistance; R+=dR)
	{
		GeometryVector temp(4);
		temp.x[0]=R;
//		temp.x[2]=dR;
		Result.push_back(temp);
	}

	//discard empty bins at the end. Such bins are useful when we don't have an upper bound of data (e.g. nearest neighbor distance binning)
	//here, they are useless
	//the bin counts
	double TotalAtomCount = 0;
	double AverageDensity = 0.0;
	for(size_t j=0; j<NumConfigs; j++)
	{
		PeriodicCellList<std::complex<double>> Config = GetConfigsFunction(j);
		Config.SetCellSize(0.35*MaxDistance);
		Config.RefineCellList();
		size_t NumParticle=Config.NumParticle();
		//pre-allocate
		Config.PrepareIterateThroughNeighbors(MaxDistance);
		std::vector<double> real(Result.size(),0.0);
		std::vector<double> imag(Result.size(),0.0);

#pragma omp parallel
		{
#pragma omp for schedule(guided)
			for(signed long i=0; i<NumParticle; i++)
			{
				Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(i), MaxDistance, 
					[&real, &imag, &i, &dR, &Config](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
				{
					if(SourceAtom != i || (SourceAtom==i && LatticeShift.Modulus2()!=0) )
					{
						double dist=std::sqrt(shift.Modulus2());
						size_t idx = (size_t)floor(dist/dR);
						if(idx<real.size())
						{
							std::complex<double> val = Config.GetCharacteristics(SourceAtom); //std::conj(Config.GetCharacteristics(i)) * 
						#pragma omp atomic
							real[idx]+= val.real();
						#pragma omp atomic
							imag[idx]+= val.imag();
						}
					}
				});
			}
		}
		// StartParticle=0;
		TotalAtomCount+=Config.NumParticle();
		AverageDensity+=Config.NumParticle()/Config.PeriodicVolume();
		pd++;

		for(size_t i=0; i<real.size(); i++)
		{
			double R = dR*(i+1);
			double factor = (HyperSphere_Volume(dim, R)-HyperSphere_Volume(dim, R-dR)) * Config.NumParticle();
			double val_real = real[i] / factor;
			double val_imag = imag[i] / factor;
			Result[i].x[1]+= val_real;  //values[i];
			Result[i].x[2]+= val_imag;  //values[i];
			Result[i].x[3]+= val_real*val_real + val_imag*val_imag; //values[i]*values[i];
		}
	}
	AverageDensity/=NumConfigs;

	for(int i = 0; i < Result.size(); i ++ ){
		double R = Result[i].x[0];
		Result[i].x[1] /= NumConfigs;
		Result[i].x[2] /= NumConfigs;
		Result[i].x[3] = std::sqrt( std::abs(Result[i].x[3] / NumConfigs - Result[i].x[1] * Result[i].x[1]));
		Result[i].x[0] -= 0.5*dR;
	}

	if(Verbosity>2)
		std::cout<<"done!\n";


}
