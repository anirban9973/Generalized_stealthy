/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: March 2022 */

/** \file Computation.h 
 * \brief Header file for a Computation class.
 * It simplified the input and output for some calculations involving an ensemble of point configurations.
 * Computation PairStatisticsCLI code written by Ge.*/


#ifndef COMPUTATION_H__
#define COMPUTATION_H__

#include <set>
#include <cmath>
#include <list>
#include <vector>
#include <omp.h>

/* header files in $(cores) */
#include <etc.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>

#include "WeightedStructureFactor.h"
#include "PairCorrelation.h"
#include "LocalVariance.h"

typedef PeriodicCellList<std::complex<double>> MarkedPointPattern;



class Computation
{
public:
	int num_threads = 1;

	virtual void Compute(std::function<const MarkedPointPattern (size_t i)> GetConfigsFunction, size_t NumConfig) =0;
	virtual void Write(const std::string OutputPrefix) =0;
	virtual void Plot(const std::string OutputPrefix, const std::string & Title) =0;

	//some computations allow additional options
	//any command not understood by PairStatisticsCLI() is passed to here
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}

	virtual ~Computation()
	{
	}

	void SetNumThreads(int n){
		this->num_threads = n;
	}
};

/* class to compute pair correlation functions as a function of radius $r$ */
class g2Computation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double g2rmax, dr;
	double resolution;
	bool single_sided = false;
	
	g2Computation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"g2 R_max=";
		ifile>>g2rmax;
		ofile<<"dr = ";
		ifile>>dr;
//		resolution=1.0;
	}

	virtual void Compute(std::function<const MarkedPointPattern(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	virtual ~g2Computation()
	{
	}
};


/* a class to compute the structure factor as a function of a wavenumber */
class SkComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double CircularKMax, LinearKMax, KPrecision, CircularKMin = 0.0;
	double SampleProbability;
	size_t average_option = 0;
	
	SkComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"circular K_max=";
		ifile>>CircularKMax;
		ofile<<"linear K_max=";
		ifile>>LinearKMax;
		ofile<<"K precision (binning width)=";
		ifile>>KPrecision;

		SampleProbability = 1.0;
	}

	virtual void Compute(std::function<const MarkedPointPattern(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


/* a class to compute the local number variance as a function of a wavenumber */
class LocalVariance : public Computation
{
public:
	std::vector<GeometryVector> result;
	std::vector<GeometryVector> mean_NR; 
	double Rmax, dR, Rmin = -1.0;
	size_t num_samp_centers, seed;
	RandomGenerator rng;
	bool use_prt_centers = false;
	LocalVariance(std::istream & ifile, std::ostream & ofile)
	{
		ofile << "Rmax = ";
		ifile >> Rmax;
		ofile << "dR = ";
		ifile >> dR;
		ofile << "number of random sampling centers = ";
		ifile >> num_samp_centers;
		ofile << "random seed = ";
		ifile >> seed;
		rng.seed(seed);
	}

	virtual void Compute(std::function<const MarkedPointPattern(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


#endif // COMPUTATION_H__