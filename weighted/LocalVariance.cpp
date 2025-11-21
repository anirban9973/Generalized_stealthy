/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2025 */

/** \file LocalVariance.cpp
 *	\brief Function implementations of computing local number variances of point patterns under the periodic boundary conditions. */

#include "LocalVariance.h"
//#include <time.h>

double Rmin_local_variance = -1.0;

void GetSamplingPoints(DimensionType dim, size_t NumPoints, std::vector<GeometryVector> &centers, RandomGenerator & rng)
{
	centers.clear();
	centers.resize(NumPoints, GeometryVector(static_cast<int> (dim)));

	for (int i = 0; i < NumPoints; i ++){
		for (int j = 0; j < dim; j++){
			centers[i].x[j] = rng.RandomDouble();
		}
	}
}

void GetF(const MarkedPointPattern &config, std::vector<std::complex<double>> &N, double Rmax, double dR, const GeometryVector &x0_RelativeCoordinates, const std::function<double(const GeometryVector &)> GetRadius)
{
	assert(x0_RelativeCoordinates.Dimension == config.GetDimension()); 
	N.clear();
	N.resize((int)ceil(Rmax / dR), std::complex<double> (0.0,0.0));
	// std::vector<double> RE(N.size(), 0.0);
	// std::vector<double> IM(N.size(), 0.0);
	//Count the number of particles within a shell of radius (R-dR, R)
	double Rc = sqrt(config.GetDimension())*Rmax;
	if (Rmin_local_variance < 0.0){
		config.IterateThroughNeighbors(x0_RelativeCoordinates, Rc, [dR, &N, &GetRadius, &config](const GeometryVector &x0, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void {
			int index = (int)floor(GetRadius(x0) / dR);
			if (index < N.size()){
				// RE[index] +=config.GetCharacteristics(prt_index).real();
				// IM[index] +=config.GetCharacteristics(prt_index).imag();
				N[index] +=config.GetCharacteristics(prt_index);
			}
		}
		);
	}
	else{
		/* Exclude the particles within the Rmin_local_variance */
		config.IterateThroughNeighbors(x0_RelativeCoordinates, Rc, [dR, &N, &GetRadius, &config, &Rmin_local_variance](const GeometryVector &x0, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void {
			double R = GetRadius(x0);
			if (R > Rmin_local_variance)
			{
				int index = (int)floor(GetRadius(x0) / dR);
				if (index < N.size()){
					// RE[index] +=config.GetCharacteristics(prt_index).real();
					// IM[index] +=config.GetCharacteristics(prt_index).imag();
					N[index] +=config.GetCharacteristics(prt_index);
				}
			}
		}
		);
	}


	/* cumulative sum */
	for (int i = 1; i < N.size(); i++) {
		//N[i] += std::complex<double> (RE[i-1], IM[i-1]);  //N[i-1];
		N[i] += N[i-1]; 
	}
}


void MC_LV_single(const MarkedPointPattern & config, std::vector<GeometryVector> & results, double Rmax, double dR, const std::vector<GeometryVector> & window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct)
{
	results.clear();
	results.resize((int)ceil(Rmax / dR), GeometryVector(0,0,0,0));
	for (int i = 0; i < results.size(); i++){
		results[i].x[0] = dR*(i+1.0);
	}
	std::vector<double> sum1_RE(results.size(), 0.0);
	std::vector<double> sum1_IM(results.size(), 0.0);
	std::vector<double> sum2(results.size(), 0.0);
//	config.PrepareIterateThroughNeighbors(Rmax);

	/* second and other centers */
#pragma omp parallel for
	for (int i=0; i< window_centers.size(); i++){
		std::vector<std::complex<double>> N;
		MarkedPointPattern c(config);
		c.PrepareIterateThroughNeighbors(Rmax);

		GetF(c, N, Rmax, dR, window_centers[i], WindowSizeFunct);

		for (int j=0; j<N.size(); j++){
			double x = N[j].real();
			double y = N[j].imag();
			double z = N[j].real()*N[j].real() + N[j].imag()*N[j].imag();
#pragma omp atomic
			sum1_RE[j] += x;
#pragma omp atomic
			sum1_IM[j] += y;
#pragma omp atomic
			sum2[j] += z;
		}
	}

	/* summarize data */
	double num_window_centers = (double)window_centers.size();
	for (int i = 0; i < results.size(); i++){
		std::complex<double> mean = std::complex<double>(sum1_RE[i], sum1_IM[i]) / num_window_centers; 
		double variance = sum2[i]/num_window_centers - (mean.real()*mean.real() + mean.imag()*mean.imag());
		
		// double variance = (results[i].x[2]/num_window_centers - mean*mean) ;
		// if (num_window_centers > 1.0){
		// 	variance *= num_window_centers/(num_window_centers-1.0) ;
		// }

		results[i].x[1] = variance;
		results[i].x[2] = mean.real(); //0.0L;
		results[i].x[3] = mean.imag(); //0.0L;
	}
}

void MC_LV_single(const MarkedPointPattern &config, std::vector<GeometryVector> &results, double Rmax, double dR, size_t num_window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct, int seed)
{
	std::vector<GeometryVector> centers;
	RandomGenerator rng (seed);
	GetSamplingPoints(config.GetDimension(), num_window_centers, centers, rng);
	MC_LV_single(config, results, Rmax, dR, centers, WindowSizeFunct);
}

void MC_LV_Ensemble(const std::function<MarkedPointPattern(size_t i)> &GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> &Result, double Rmax, double dR, const std::vector<GeometryVector> &WindowCenters, std::vector<GeometryVector> & means, const std::function<double(const GeometryVector &)> WindowSizeFunct)
{
	Result.clear();
	std::cout<< "\nComputing Local Number Variance (given window centers):";
	progress_display pd(NumConfigs);

	std::vector<GeometryVector> mean_NR;
	/* repeat over configurations */
	for (int i = 0; i < NumConfigs; i++){
		MarkedPointPattern c = GetConfigsFunction(i);
		std::vector<GeometryVector> Result_;
		MC_LV_single(c, Result_, Rmax, dR, WindowCenters, WindowSizeFunct);
		pd ++ ;
		if (i == 0){
			Result = std::vector<GeometryVector> (Result_);
			mean_NR.reserve(Result_.size());
			for (int j = 0; j < Result.size(); j++){
				Result[j].x[3] = Result[j].x[1];
				Result[j].x[1] = 0.0;

				mean_NR.push_back(GeometryVector(Result_[j].x[2],Result_[j].x[3],0.0,0.0));
			}
		}
		else{
			for (int j = 0; j < Result.size(); j++){
				double val = Result_[j].x[1] - Result[j].x[3];
				Result[j].x[1] += val;
				Result[j].x[2] += val*val;

				mean_NR[j].x[0] += Result_[j].x[2];
				mean_NR[j].x[1] += Result_[j].x[3];
				mean_NR[j].x[2] += Result_[j].x[2]*Result_[j].x[2];
				mean_NR[j].x[3] += Result_[j].x[3]*Result_[j].x[3];
			}
		}
	}
	/* Summarize data */
	for(int i = 0; i < Result.size(); i++){
		Result[i].x[1] /= (double)NumConfigs;
		/*SE*/
		if(NumConfigs > 1){
			Result[i].x[2] = std::sqrt(abs(Result[i].x[2]/NumConfigs - Result[i].x[1]*Result[i].x[1])/(NumConfigs - 1));	
		}
		else{
			Result[i].x[2] = 0.0;
		}
		double mean = Result[i].x[1];
		Result[i].x[1] += Result[i].x[3];
		Result[i].x[3] = 0.0;

		/* Mean of N_f(R) */
		mean_NR[i].x[0] /= (double)NumConfigs;
		mean_NR[i].x[1] /= (double)NumConfigs;
		mean_NR[i].x[2] = (mean_NR[i].x[2]/(double)NumConfigs - mean_NR[i].x[0]*mean_NR[i].x[0]) / ((double)NumConfigs -1.0);
		mean_NR[i].x[3] = (mean_NR[i].x[3]/(double)NumConfigs - mean_NR[i].x[1]*mean_NR[i].x[1]) / ((double)NumConfigs -1.0);
		for (int j = 2; j < 4; j++){
			if (mean_NR[i].x[j] < 0.0){
				mean_NR[i].x[j] = 0.0; // avoid negative variance
			}
			else{
				mean_NR[i].x[j] = std::sqrt(mean_NR[i].x[j]);
			}
		} 
	}

	if(true){
		means.swap(mean_NR);
		std::cout << "Mean values are computed.\n";
	}

	std::cout << "Done!\n";
}

void MC_LV_Ensemble(const std::function<MarkedPointPattern(size_t i)> &GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> &Result, double Rmax, double dR, size_t num_centers, std::vector<GeometryVector> & means, bool use_particle_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct)
{
	Result.clear();
	std::cout<< "\nComputing Local Number Variance (given window centers):";
	progress_display pd(NumConfigs);
	RandomGenerator rng(12345); // fixed seed for reproducibility

	std::vector<GeometryVector> mean_NR;
	/* repeat over configurations */
	for (int i = 0; i < NumConfigs; i++){
		MarkedPointPattern c = GetConfigsFunction(i);
		std::vector<GeometryVector> Result_;
		/* Create window centers */
		std::vector<GeometryVector> WindowCenters;
		if (use_particle_centers){
			WindowCenters.reserve(c.NumParticle());
			for (int j = 0; j < c.NumParticle(); j++){
				GeometryVector v = c.GetRelativeCoordinates(j);
				WindowCenters.push_back(v);
			}
		}
		else{
			GetSamplingPoints(c.GetDimension(), num_centers, WindowCenters, rng);
		}

		MC_LV_single(c, Result_, Rmax, dR, WindowCenters, WindowSizeFunct);
		pd ++ ;
		if (i == 0){
			Result = std::vector<GeometryVector> (Result_);
			mean_NR.reserve(Result_.size());
			for (int j = 0; j < Result.size(); j++){
				Result[j].x[3] = Result[j].x[1];
				Result[j].x[1] = 0.0;

				mean_NR.push_back(GeometryVector(Result_[j].x[2],Result_[j].x[3],0.0,0.0));
			}
		}
		else{
			for (int j = 0; j < Result.size(); j++){
				double val = Result_[j].x[1] - Result[j].x[3];
				Result[j].x[1] += val;
				Result[j].x[2] += val*val;

				mean_NR[j].x[0] += Result_[j].x[2];
				mean_NR[j].x[1] += Result_[j].x[3];
				mean_NR[j].x[2] += Result_[j].x[2]*Result_[j].x[2];
				mean_NR[j].x[3] += Result_[j].x[3]*Result_[j].x[3];
			}
		}
	}
	/* Summarize data */
	for(int i = 0; i < Result.size(); i++){
		Result[i].x[1] /= (double)NumConfigs;
		/*SE*/
		if(NumConfigs > 1){
			Result[i].x[2] = std::sqrt(abs(Result[i].x[2]/NumConfigs - Result[i].x[1]*Result[i].x[1])/(NumConfigs - 1));	
		}
		else{
			Result[i].x[2] = 0.0;
		}
		double mean = Result[i].x[1];
		Result[i].x[1] += Result[i].x[3];
		Result[i].x[3] = 0.0;

		/* Mean of N_f(R) */
		mean_NR[i].x[0] /= (double)NumConfigs;
		mean_NR[i].x[1] /= (double)NumConfigs;
		mean_NR[i].x[2] = (mean_NR[i].x[2]/(double)NumConfigs - mean_NR[i].x[0]*mean_NR[i].x[0]) / ((double)NumConfigs -1.0);
		mean_NR[i].x[3] = (mean_NR[i].x[3]/(double)NumConfigs - mean_NR[i].x[1]*mean_NR[i].x[1]) / ((double)NumConfigs -1.0);
		for (int j = 2; j < 4; j++){
			if (mean_NR[i].x[j] < 0.0){
				mean_NR[i].x[j] = 0.0; // avoid negative variance
			}
			else{
				mean_NR[i].x[j] = std::sqrt(mean_NR[i].x[j]);
			}
		} 
	}

	if(true){
		means.swap(mean_NR);
		std::cout << "Mean values are computed.\n";
		std::cout << mean_NR.size() << ", " << means.size() << "\n";
	}

	std::cout << "Done!\n";


}
