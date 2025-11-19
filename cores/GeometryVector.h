/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Created	:	2015? 
 *	Comments: Jaeuk Kim (November 2025).
 */

/** \file GeometryVector.h
 *	\brief Header file to define vectors in Euclidean space and their operatios. 
 * 			(November 2025)
 *	Some type-related warnings are resolved. */


#ifndef GEOMETRYVECTOR_INCLUDED
#define GEOMETRYVECTOR_INCLUDED

#include "etc.h"

#include <cstring>
#include <cassert>
#include <fstream>
#include <cmath>
#include "RandomGenerator.h"

#ifndef NDEBUG
#define GEOMETRYVECTOR_RECORDDIMENSION
#endif


//In debug mode, this class also records the dimension of the GeometryVector.
/** \brief Fixed-size vector in Euclidean space supporting basic arithmetic. */
class GeometryVector
{
public:
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
	DimensionType Dimension;
#endif

	/** \brief Record the current dimension (debug builds only). */
	void SetDimension(DimensionType d)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=d;
#endif
	}

	double x[ ::MaxDimension];

	/** \brief Default-construct a zero vector of unspecified dimension. */
	GeometryVector()
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=0;
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=0;
	}
	//do not change this "int" to "DimensionType"!
	//otherwize GeometryVector(2) would be ambiguous.

	/** \brief Construct a zero vector with an explicit dimension. */
	GeometryVector(int dimension)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=dimension;
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=0.0L;
	}

	/** \brief Construct a 1D vector. */
	GeometryVector(double xx)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=1;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		for(int i=1; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Construct a 2D vector. */
	GeometryVector(double xx, double yy)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=2;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		this->x[1]=yy;
		for(int i=2; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Construct a 3D vector. */
	GeometryVector(double xx, double yy, double zz)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=3;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		this->x[1]=yy;
		this->x[2]=zz;
		for(int i=3; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Construct a 4D vector. */
	GeometryVector(double xx, double yy, double zz, double aa)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		this->Dimension=4;
		assert(this->Dimension<=::MaxDimension);
#endif
		this->x[0]=xx;
		this->x[1]=yy;
		this->x[2]=zz;
		this->x[3]=aa;
		for(int i=4; i< ::MaxDimension; i++)
			this->x[i]=0;
	}

	/** \brief Copy-construct from another vector. */
	GeometryVector(const GeometryVector & src)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(src.Dimension<= ::MaxDimension);
		this->Dimension=src.Dimension;
#endif
		//std::memcpy(this->x, src.x, ::MaxDimension*sizeof(double));
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=src.x[i];
	}

	/** \brief Add components from another vector in-place. */
	void AddFrom(const GeometryVector & right)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]+=right.x[i];
	}

	/** \brief Subtract components from another vector in-place. */
	void MinusFrom(const GeometryVector & right)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]-=right.x[i];
	}

	/** \brief Scale all components by a constant. */
	void MultiplyFrom(const double & right)
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(this->Dimension<= ::MaxDimension);
#endif
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]*=right;
	}

	/** \brief Inner product with another vector. */
	double Dot(const GeometryVector & right) const		// equivalent to the inner product
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
#endif
		double result=0;
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		for(int i=0; i<this->Dimension; i++)
#else
		for(int i=0; i< ::MaxDimension; i++)
#endif
				result+=this->x[i]*right.x[i];
		return result;
	}
	//cross product

	/** \brief 3D cross product. */
	GeometryVector Cross(const GeometryVector & right) const
	{
		GeometryVector result((int)3);
	//#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension == 3);
		assert(this->Dimension == 3);
	//#endif

		result.x[0] = this->x[1] * right.x[2] - this->x[2] * right.x[1];
		result.x[1] = this->x[2] * right.x[0] - this->x[0] * right.x[2];
		result.x[2] = this->x[0] * right.x[1] - this->x[1] * right.x[0];

		return result;
	}


	/** \brief Squared Euclidean norm. */
	double Modulus2(void) const
	{
		return this->Dot(*this);
	}

	/** \brief Exact component-wise equality test. */
	bool IsEqual(const GeometryVector & right) const
	{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
		assert(right.Dimension==this->Dimension);
		assert(right.Dimension<= ::MaxDimension);
		for(int i=0; i<right.Dimension; i++)
			if(this->x[i]!=right.x[i])
				return false;
#else
		for(int i=0; i< ::MaxDimension; i++)
			if(this->x[i]!=right.x[i])
				return false;
#endif
		return true;
	}

	/** \brief Stream out the first dim coordinates separated by tabs. */
	void OutputCoordinate(std::ostream & os, DimensionType dim) const
	{
		for(int i=0; i< dim; i++)
			os<<this->x[i]<<" \t";
	}

	/** \brief Read dim coordinates from a stream. */
	void InputCoordinate(std::istream & is, DimensionType dim)
	{
		for(int i=0; i< dim; i++)
			is>>this->x[i];
	}

	/** \brief Write coordinates in binary form. */
	void WriteBinary(std::ostream & ofile, DimensionType dimension) const;

	/** \brief Read coordinates from binary form. */
	void ReadBinary(std::istream & ifile, DimensionType dimension);

	friend inline GeometryVector SameDimensionZeroVector(const GeometryVector & src);
};
//create another GeometryVector with same dimension as src, but the value is zero
inline GeometryVector SameDimensionZeroVector(const GeometryVector & src)
{
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
	return GeometryVector(src.Dimension);
#else
	return GeometryVector(0);
#endif
}



inline std::ostream & operator << (std::ostream & os, const GeometryVector & a)
{
	//os<<a.Dimension<<" \t";
	a.OutputCoordinate(os, ::MaxDimension);
	os<<'\n';
	return os;
}
inline std::istream & operator >> (std::istream & is, GeometryVector & a)
{
	//is>>a.Dimension;
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
	assert(a.Dimension<= ::MaxDimension);
#endif
	a.InputCoordinate(is, ::MaxDimension);
	return is;
}
inline GeometryVector operator + (const GeometryVector & left, const GeometryVector & right)
{
	GeometryVector result(left);
	result.AddFrom(right);
	return result;
}
inline GeometryVector operator - (const GeometryVector & left, const GeometryVector & right)
{
	GeometryVector result(left);
	result.MinusFrom(right);
	return result;
}
inline GeometryVector operator * (const GeometryVector & left, const double & right)
{
	GeometryVector result(left);
	result.MultiplyFrom(right);
	return result;
}
inline GeometryVector operator * (const double & left, const GeometryVector & right)
{
	GeometryVector result(right);
	result.MultiplyFrom(left);
	return result;
}
inline bool operator == (const GeometryVector & left, const GeometryVector & right)
{
	return left.IsEqual(right);
}
inline bool operator != (const GeometryVector & left, const GeometryVector & right)
{
	return !left.IsEqual(right);
}

//random unitvector in a unit cube
inline void RandomVector(DimensionType d, RandomGenerator & gen, GeometryVector & result) {
	result.SetDimension(d);
	for (DimensionType i = 0; i < d; i++) {
		result.x[i] = 2.0*(gen.RandomDouble()-0.5);
	}
}

inline void RandomVector(DimensionType d, GeometryVector & result) {
	result.SetDimension(d);
	RandomGenerator gen;
	for (DimensionType i = 0; i < d; i++) {
		result.x[i] = 2.0*(gen.RandomDouble() - 0.5);
	}
}

GeometryVector RandomUnitVector(DimensionType d);
GeometryVector RandomUnitVector(DimensionType d, RandomGenerator & gen);
double Volume (const GeometryVector * vecs, DimensionType dimension);
double SimplexVolume (const GeometryVector * vecs, DimensionType d);
void StandardizeOrientation(GeometryVector * bas, DimensionType d);

void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields = std::vector<std::string>());
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields = std::vector<std::string>());
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns);
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns);


//another implementation of vector
//the difference between GeometryVector and BigVector is that BigVectors can have dimensions exceeding MaxDimension
//BigVector can be used to represent vectors in configuration space
class BigVector
{
public:
	double * x;
	size_t dim;
	BigVector()
	{
		this->dim=0;
		this->x=nullptr;
	}
	BigVector(size_t dimension)
	{
		this->dim=dimension;
		this->x=new double[this->dim];
		for(size_t i=0; i<dimension; i++)
			this->x[i]=0.0;
	}
	BigVector(const BigVector & src)
	{
		this->dim=src.dim;
		this->x=new double[this->dim];
		std::memcpy(this->x, src.x, sizeof(double)*this->dim);
	}
	BigVector(BigVector && src)
	{
		this->dim=src.dim;
		this->x=src.x;
		src.dim=0;
		src.x=nullptr;
	}
	~BigVector()
	{
		if(this->x!=nullptr)
			delete [] this->x;
	}
	void AddFrom(const BigVector & right)
	{
		assert(right.dim==this->dim);
		for(size_t i=0; i<this->dim; i++)
			this->x[i]+=right.x[i];
	}
	void MinusFrom(const BigVector & right)
	{
		assert(right.dim==this->dim);
		for(size_t i=0; i<this->dim; i++)
			this->x[i]-=right.x[i];
	}
	void MultiplyFrom(const double & right)
	{
		for(size_t i=0; i<this->dim; i++)
			this->x[i]*=right;
	}
	double Dot(const BigVector & src) const
	{
		assert(src.dim==this->dim);
		double result=0.0;
		for(size_t i=0; i<this->dim; i++)
			result+=this->x[i]*src.x[i];
		return result;
	}
	double Modulus2(void) const
	{
		return this->Dot(*this);
	}
	void OutputCoordinate(std::ostream & os) const
	{
		for(size_t i=0; i< dim; i++)
			os<<this->x[i]<<" \t";
	}
	void InputCoordinate(std::istream & is)
	{
		for(size_t i=0; i< dim; i++)
			is>>this->x[i];
	}
};
inline BigVector operator + (const BigVector & left, const BigVector & right)
{
	BigVector result(left);
	result.AddFrom(right);
	return result;
}

inline BigVector operator - (const BigVector & left, const BigVector & right)
{
	BigVector result(left);
	result.MinusFrom(right);
	return result;
}

inline BigVector operator * (const BigVector & left, const double & right)
{
	BigVector result(left);
	result.MultiplyFrom(right);
	return result;
}

inline BigVector operator * (const double & left, const BigVector & right)
{
	return  right*left;
}


#endif
