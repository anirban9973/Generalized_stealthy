/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Date	:  */



#ifndef WeightedSTRUCTUREFACTOR_INCLUDED
#define WeightedSTRUCTUREFACTOR_INCLUDED

#include <functional>
#include <vector>

/* header fileds in $(cores)*/
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>
/** \file WeightedStructureFactor.h
 *	\brief Header file for computing the structure factor of point configurations with weights */


//calculate the Structure factof with weights; here config.getCharacteristics(i) = the weight on the i-th particle
double WeightedStructureFactor(const PeriodicCellList<std::complex<double>> & Config, const GeometryVector & k);

void IsotropicStructureFactor_weighted(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision = 0.01, double SampleProbability = 1.0, size_t option = 0, double CircularKMin = 0.0);


#endif