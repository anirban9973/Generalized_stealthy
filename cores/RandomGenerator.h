#ifndef RANDOMGENERATOR_INCLUDED
#define RANDOMGENERATOR_INCLUDED
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <random>
#include <functional>

/** \brief Random-number utilities for uniform, Gaussian, and Poisson draws. */
class RandomGenerator
{
private:
	boost::mt19937 gen;
	boost::mt19937::result_type min, max;
	void initialize(void);
public:
	/** \brief Seed the underlying generator with a specific integer. */
	void seed(int s);

	/** \brief Construct a generator seeded from current time. */
	RandomGenerator();

	/** \brief Construct a generator with an explicit seed value. */
	RandomGenerator(int Seed);

	/** \brief Draw a uniform deviate in [0,1). */
	double RandomDouble(void);

	/** \brief Draw a deviate using an inverse cumulative distribution. 
	 *  \param inverseCp Inverse CDF mapping from uniform [0,1) to the target distribution. */
	double RandomDouble(const std::function<double(double)> &inverseCp);

	/** \brief Draw a Gaussian deviate with mean 0 and supplied standard deviation. */
	double RandomDouble_normal(double std) {
		boost::normal_distribution<> nd(0.0, std);
		boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(this->gen, nd);
		return var_nor();
	}

	/** \brief Draw a Poisson deviate with given mean. */
	int PoissonRNG(double mean); 

	/** \brief Draw an integer from the full underlying engine range. */
	inline int RandomInt(void) {return this->gen();}
};

/** \brief Generate a random unit vector in a given dimension. 
 *  \param Dimension Dimension of the target space. 
 *  \param gen Random number generator to sample from. 
 *  \param result Output buffer to receive the unit-vector components. */
void GetRandomVector(int Dimension, RandomGenerator & gen, double * result);


#endif
