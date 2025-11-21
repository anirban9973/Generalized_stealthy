/** Author	: Ge Zhang
 *  Email	:
 *	Date	:	*/

/** \file PairCorrelation.h
 *	\brief Header file to compute g2 function for point configurations, and S2 function for sphere packings. */

#ifndef PAIRCORRELATION_INCLUDED
#define PAIRCORRELATION_INCLUDED

#include "PeriodicCellList.h"
#include <vector>
#include <complex>

//void Weighted IsotropicTwoPairCorrelation_FullAdaptive(std::vector<const Configuration *> Configs, double MaxDistance, std::vector<GeometryVector> & Result, double ResolutionPreference=1.0);

//clear Result, then fill it with the pair correlation data calculated from Config
//Result is filled with 4-dimensional GeometryVectors, the elements are:
//(r, g_2(r), \delta r, \delta g_2(r) )
//gets SampleDistanceSize number of sample pair distances, use it to generate bins, then count all pair distances
//ResolutionPreference : >1 to get better r resolution, <1 to get better g_2 resolution
/** \brief Compute an isotropic g2 function from a number of point configurations.
 * @param[in] (GetConfigsFunction, NumConfigs)	A Lambda function to generate "NumConfigs" point configurations. 
 * @param[in] MaxDistance	Up to this distance, the g2 function is computed.
 * @param[out] Result		A resulting g2 function.
							First clear this vector and store data as follows:
							Result[i] = (radius, g2(radius), error in radius, error in g2).
 * @param[in] SampleDistanceSize	The number of pairs to determine a histogram.
									The larger this argument, the finer bins the histogram has. This sampling is performed in a single thread.
 * @param[in] ResolutionPreference	>1 to get better r resolution, <1 to get better g_2 resolution. */
void WeightedIsotropicTwoPairCorrelation(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MaxDistance, std::vector<GeometryVector> & Result, double dR = 0.1);


void WeightedIsotropicSingleSidedTwoPairCorrelation(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MaxDistance, std::vector<GeometryVector> & Result, double dR = 0.1);



#endif