/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Created	:	2015? 
 *	Comments: Jaeuk Kim (November 2025).
 */

/** \file PeriodicCellList.h
 *	\brief Header file to define general point patterns under periodic boundary conditions. It and its derived classes (e.g., Configuration and SpherePacking) implement cell-list acceleration techniques. 
 * 			(November 2025)
 *	Some type-related and unintialized warnings are resolved. 
 * */


#ifndef PeriodicCellList_INCLUDED
#define PeriodicCellList_INCLUDED

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <list>

#include <omp.h>

#include "GeometryVector.h"
#include "RandomGenerator.h"
#include "etc.h"
//#include "etc1.h"

const size_t ParticlePerCell = 1;

//class for point configurations under periodic boundary conditions;
//unit cells can be any Bravais lattice structure.
//Characteristics store additional information on each particle...

//Configuration & SpehrePacking are derived from this.

/** \brief Cell-list accelerator for periodic point configurations parameterized by particle traits. */
template<typename OtherCharacteristics>
class PeriodicCellList
{
	template<typename a> friend class PeriodicCellList;

	//some concepts:
	//there are 3 types of coordinates:
	//1. Cartesian Coordinates
	//2. Relative Coordinates, this is the coordinates relative to the periodic box base vectors.
	//		The periodic box base vectors are PeriodicCellList::BasisVector
	//3. Lattice Coordinates, this is the coordinates relative to the cell lattice(i.e. the base vectors of a cell list cell).
public:
	typedef unsigned long Index;

protected:
	mutable double MaxRadius;//=0 indicates this variable is not valid. call UpdateMaxRadius to update it //moved from SpherePacking
	DimensionType Dimension;
	GeometryVector BasisVector[::MaxDimension];//the base vector of a periodic cell
	signed long CellRank[::MaxDimension];//number of cells along each basis vector
	std::vector<std::vector<size_t > > Cells;//content Particles in the cells
	std::vector<GeometryVector> ParticleRelatives;
	std::vector<GeometryVector> ParticleCartesians;
	std::vector<OtherCharacteristics> ParticleCharacteristics;
	double ApproximateCellSize;// the best choice of cell size obtained in constructor when user input cell size and cell rank. We will adjust CellRank so that their size is roughly equal to this quantity.
	double vol; // the volume of a fundamental cell: add recently

	//reciprocal basis vector, calculated when needed
	mutable bool ReciprocalBasisVectorValid[::MaxDimension];
	mutable GeometryVector ReciprocalBasisVector[::MaxDimension];

	//member data for vicinity lattice
	struct VicinityLatticesMember
	{
		GeometryVector CartesianCoordinates;
		signed long LatticeShift[MaxDimension];
	};
	mutable std::vector<VicinityLatticesMember> VicinityLattices;
	mutable std::map<double, size_t> VicinityLatticePartitions;
	mutable double VicinityLatticesRc2;
	bool SortedVicinityList;//use sorted vicinity lattice list, which takes some time to sort, but makes iteration more efficient

	void UpdateCartesianCoordinates(size_t a)
	{
		if (!(this->ParticleRelatives[a].x[0] < 1 + ::LengthPrecision) || !(this->ParticleRelatives[a].x[0] >= 0 - ::LengthPrecision))
		{
			std::cerr << "Error in PeriodicCellList::UpdateCartesianCoordinates : invalid ParticleRelatives. ParticleRelatives[" << a << "]=" << this->ParticleRelatives[a];
			exit(0);
		}
		this->ParticleCartesians[a] = GeometryVector(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->ParticleCartesians[a].AddFrom((this->BasisVector[i])*(this->ParticleRelatives[a].x[i]));
	}
	void PutToPeriodicBox(size_t a)
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->ParticleRelatives[a].x[i] -= std::floor(this->ParticleRelatives[a].x[i]);
	}
	Index ConvertIndex(const signed long * index) const//convert dimensional index to the 1-dimension cell index
	{
		Index result = index[0];
		for (DimensionType i = 1; i<this->Dimension; i++)
		{
			assert(index[i] >= 0);
			assert(index[i]<this->CellRank[i]);
			result *= this->CellRank[i];
			result += index[i];
		}
		return result;
	}
	void GetFullIndex(signed long * index, const GeometryVector & RelativeCoordinates) const
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			index[i] = (signed long)(std::floor(RelativeCoordinates.x[i] * this->CellRank[i]));
			if (index[i]<0 || index[i] >= this->CellRank[i])
			{
				std::cerr << "Error in Configuration::GetFullIndex : index invalid. index=" << index[i] << ", Relative Coordinate=" << RelativeCoordinates.x[i] << ", Cell Rank=" << this->CellRank[i] << '\n';
				assert(false);
			}
		}
	}


	Index GetIndex(const size_t a) const //get the index of the cell that a should be placed in, the LatticeCoordinates of particle must be in box
	{
		signed long index[::MaxDimension] = {0};
		this->GetFullIndex(index, this->ParticleRelatives[a]);
		return this->ConvertIndex(index);
	}


	bool UpdateCellList(std::vector<size_t>::iterator a, std::vector<size_t> & from)//move a to another cell if necessary, from is the cell that a comes from, return true if moved
	{
		//calculate new index
		Index newindex = this->GetIndex((*a));
		if (&(from) == &(this->Cells[newindex]))//don't need to move
			return false;

		//add to new cell
		this->Cells[newindex].push_back(*a);

		//remove from original cell
		*a = from.back();
		from.pop_back();
		return true;
	}


	bool _VicinityLatticeIteration(GeometryVector now, DimensionType DimIDo, GeometryVector nowRelative, GeometryVector * ortho, double Rmax2, double Rmin2) const
	{
		const double VicinityLatticeRmin_Tolerance = 1e-10;
		if (DimIDo == 0)
		{
			double temp = now.Modulus2();
			if (temp<Rmax2)
			{
				bool ToAdd;
				if (temp >(1.0 + VicinityLatticeRmin_Tolerance)*Rmin2)
					ToAdd = true;
				else if (temp < (1.0 - VicinityLatticeRmin_Tolerance)*Rmin2)
					ToAdd = false;
				else
				{
					//cannot tell if the lattice point was already included just from the distance
					//instead, go through the entire VicinityLattices vector to see if it is already added
					VicinityLatticesMember temp;
					for (DimensionType i = 0; i<this->Dimension; i++)
					{
						double & x = nowRelative.x[i];
						temp.LatticeShift[i] = x >= 0 ? static_cast<signed long>(x + 0.5) : static_cast<signed long>(x - 0.5);
					}
					ToAdd = true;
					for (auto iter = this->VicinityLattices.begin(); iter != this->VicinityLattices.end(); iter++)
					{
						bool LatticeShiftIdentical = true;
						for (DimensionType d = 0; d<this->Dimension; d++)
							if (iter->LatticeShift[d] != temp.LatticeShift[d])
							{
								LatticeShiftIdentical = false;
								break;
							}
						if (LatticeShiftIdentical)
						{
							ToAdd = false;
							break;
						}
					}
				}

				if (ToAdd)
				{
					VicinityLatticesMember temp;
					temp.CartesianCoordinates = now;
					for (DimensionType i = 0; i<this->Dimension; i++)
					{
						double & x = nowRelative.x[i];
						temp.LatticeShift[i] = x >= 0 ? static_cast<signed long>(x + 0.5) : static_cast<signed long>(x - 0.5);
					}

					this->VicinityLattices.push_back(temp);
				}
				return true;
			}
			else
				return false;
		}


		GeometryVector k = now;
		for (DimensionType i = 0; i<DimIDo; i++)
		{
			GeometryVector & e = ortho[i];
			k.MinusFrom(e*(k.Dot(e) / e.Dot(e)));
		}
		if (k.Dot(k)>Rmax2)
			return false;
		//right now k is the normal vector of plane (e1, e2,...) in subspace (e1, e2, ... , now)
		GeometryVector a = now - k;//the part of now in that plane
		GeometryVector ethis = this->BasisVector[DimIDo - 1] * (1.0 / this->CellRank[DimIDo - 1]);//base vector that this function will loop through
		for (DimensionType i = 0; i<DimIDo - 1; i++)
		{
			GeometryVector & e = ortho[i];
			a.MinusFrom(e*(a.Dot(e) / e.Dot(e)));
			ethis.MinusFrom(e*(ethis.Dot(e) / e.Dot(e)));
		}
		double coeff = 0;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (std::abs(ethis.x[i])> ::LengthPrecision)
			{
				coeff = a.x[i] / ethis.x[i];
				break;
			}
		}//a equals to ethis*coeff
		assert((a - coeff * ethis).Modulus2()< ::LengthPrecision);

		double CoeffMaxInPlane = std::sqrt((Rmax2 - k.Modulus2()) / ethis.Modulus2());//sqrt(R^2-k^2)/|ethis|
		signed long imin = (signed long)(std::ceil((-1)*CoeffMaxInPlane - coeff));
		signed long imax = (signed long)(std::floor(CoeffMaxInPlane - coeff));

		//int nearest=std::floor(coeff);
		GeometryVector Unit(this->Dimension);
		Unit.x[DimIDo - 1] = 1;
		for (signed long i = imin; i <= imax; i++)
		{
			this->_VicinityLatticeIteration(now + i * this->BasisVector[DimIDo - 1] * (1.0 / this->CellRank[DimIDo - 1]), DimIDo - 1, nowRelative + i * Unit, ortho, Rmax2, Rmin2);
		}
		return true;
	}

public:
	//This class can work in one of the two modes:
	//UseSortedList=true : efficient when the basis vectors almost never change.
	//UseSortedList=false : efficient when the basis vectors change.

	/** \brief Construct an empty list with no dimension set. */
	PeriodicCellList()
	{
		//this->VicinityLatticesRc2=0;
		this->Dimension = 0;
		Index CellCount = 1;//need to calculate overall number of cells
		this->Cells.assign(CellCount, std::vector<size_t >());
		this->MaxRadius = 0.0;
		this->ApproximateCellSize = ::MaxDimension;
		this->SortedVicinityList = false;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->vol = -1.0;
	}

	/** \brief Construct a periodic cell list for a given box geometry.
	 *  \param Dimension Spatial dimension.
	 *  \param BasisVectors Basis vectors defining the periodic box.
	 *  \param CellSize Target linear cell size.
	 *  \param UseSortedList Choose sorted vicinity list for faster neighbor checks when geometry is fixed. */
	PeriodicCellList(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList = false)
	{
		//this->VicinityLatticesRc2=0;
		this->Dimension = Dimension;
		Index CellCount = 1;//need to calculate overall number of cells
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			GeometryVector & vec = BasisVectors[i];
			this->CellRank[i] = (signed long)(std::floor(std::sqrt(vec.Modulus2()) / CellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;

			//error check
			if (this->CellRank[i] <= 0)
			{
				std::cerr << "error in Configuration constructor: invalid cell rank. basis vector length squared=" << vec.Modulus2() << ", CellSize=" << CellSize << '\n';
				assert(false);
			}
			this->BasisVector[i] = vec;
			CellCount *= this->CellRank[i];
		}
		this->MaxRadius = 0.0;
		for (DimensionType i = this->Dimension; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = 0;
			this->BasisVector[i] = GeometryVector();
		}
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = CellSize;
		this->SortedVicinityList = UseSortedList;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->vol = ::Volume(BasisVectors, Dimension);
	}

	/** \brief Copy-construct, duplicating geometry and particle data. */
	PeriodicCellList(const PeriodicCellList<OtherCharacteristics> & source)
		: MaxRadius(source.MaxRadius), Dimension(source.Dimension),
		Cells(source.Cells), ApproximateCellSize(source.ApproximateCellSize),
		SortedVicinityList(source.SortedVicinityList), VicinityLatticesRc2(0.0),
		ParticleCartesians(source.ParticleCartesians), ParticleRelatives(source.ParticleRelatives),
		ParticleCharacteristics(source.ParticleCharacteristics)
	{
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->BasisVector[i] = source.BasisVector[i];
		}
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = source.CellRank[i];
		}
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
	}

	/** \brief Construct from another list while assigning a uniform characteristic.
	 *  \param add Source characteristic type.
	 *  \param cha Characteristic value assigned to every particle. */
	template<typename add> PeriodicCellList(const PeriodicCellList<add> & source, const OtherCharacteristics & cha)
		: Dimension(source.Dimension),
		Cells(source.Cells), ApproximateCellSize(source.ApproximateCellSize),
		SortedVicinityList(source.SortedVicinityList), VicinityLatticesRc2(0.0),
		ParticleCartesians(source.ParticleCartesians), ParticleRelatives(source.ParticleRelatives),
		ParticleCharacteristics(source.NumParticle(), cha)
	{
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->BasisVector[i] = source.BasisVector[i];
		}
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = source.CellRank[i];
		}
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
	}


	/** \brief Trivial destructor. */
	virtual ~PeriodicCellList() = default;

	/** \brief Exact comparison of dimension, basis, and particle layout. */
	bool operator == (const PeriodicCellList<OtherCharacteristics> & right) const
	{
		if (right.Dimension != this->Dimension)
			return false;
		if (right.NumParticle() != this->NumParticle())
			return false;
		for (DimensionType d = 0; d<this->Dimension; d++){
			if (right.BasisVector[d] != this->BasisVector[d])
				return false;
		}
		for (size_t i = 0; i<this->NumParticle(); i++){
			if (right.ParticleRelatives[i] != this->ParticleRelatives[i])
				return false;
		}
		return true;
	}

	/** \brief Reset cached vicinity lattice data. */
	void ClearVicinityLatticeList(void)
	{
		this->VicinityLattices.clear();
		this->VicinityLatticesRc2 = 0.0;
		this->VicinityLatticePartitions.clear();
		this->VicinityLatticePartitions.insert(std::make_pair((double)(0.0), (size_t)(0)));
	}


	void TryRefineBasisVectors_inner(int ilongest)
	{
		//find the other basis vector
		int iother;
		signed int movepart = 0;
		for (iother = 0; iother<this->Dimension; iother++)
		{
			if (iother == ilongest) continue;
			movepart = (signed int)(std::floor(0.5 + (this->BasisVector[ilongest].Dot(this->BasisVector[iother])) / (this->BasisVector[iother].Dot(this->BasisVector[iother]))));
			if (movepart != 0)
				break;
		}
		if (movepart == 0)
			return;

		this->BasisVector[ilongest] = this->BasisVector[ilongest] - static_cast<double>(movepart)*this->BasisVector[iother];
		for (size_t i = 0; i<this->ParticleRelatives.size(); i++)
		{
			this->ParticleRelatives[i].x[iother] += this->ParticleRelatives[i].x[ilongest] * movepart;
			this->PutToPeriodicBox(i);
			this->UpdateCartesianCoordinates(i);
		}
		this->ClearVicinityLatticeList();
		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
		size_t NumCells = 1;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double length = std::sqrt(this->BasisVector[i].Modulus2());
			this->CellRank[i] = (signed long)(std::floor(length / this->ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			NumCells *= this->CellRank[i];
		}
		this->Cells.resize(NumCells);
		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->Cells[this->GetIndex(i)].push_back(i);
		}
	}


	void TryRefineBasisVectors(void)
	{
		if (this->Dimension == 1)
			return;
		for (int i = 0; i<this->Dimension; i++)
			this->TryRefineBasisVectors_inner(i);
	}


	double PeriodicVolume(void) const
	{
		//double result = ::Volume(this->BasisVector, this->Dimension);
		//return result;
		double result = 0;
		if (this->vol <= 0.0) {
			result = ::Volume(this->BasisVector, this->Dimension);
			//this->vol = result;
		}
		else {
			result = this->vol;
		}
		return result;
	}

	/** \brief Return the dimensionality of the configuration box. */
	DimensionType GetDimension(void) const
	{
		return this->Dimension;
	}

	/** \brief Insert a particle with a characteristics at the given relative coordinate (wrapped to the box). */
	virtual void Insert(OtherCharacteristics cha, const GeometryVector & RelativeCoordinate)
	{
		//this->Particles.push_back(particle(this->Dimension, cha));
		this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
		this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
		this->ParticleCharacteristics.push_back(cha);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			assert(RelativeCoordinate.x[i] == RelativeCoordinate.x[i]);
			//this->Particles.back().RelativeCoordinates.x[i]=RelativeCoordinate.x[i]-std::floor(RelativeCoordinate.x[i]);
			//this->Particles.back().RelativeCoordinates.x[i]-=std::floor(this->Particles.back().RelativeCoordinates.x[i]);

			//NEED TO subtract the floor TWICE
			this->ParticleRelatives.back().x[i] = RelativeCoordinate.x[i] - std::floor(RelativeCoordinate.x[i]);
			this->ParticleRelatives.back().x[i] -= std::floor(this->ParticleRelatives.back().x[i]);
		}
		this->UpdateCartesianCoordinates(this->NumParticle() - 1);
		//---test code
		int test = this->GetIndex(this->NumParticle() - 1);
		this->Cells[test].push_back(this->NumParticle() - 1);
		//this->Cells[this->GetIndex(this->NumParticle() - 1)].push_back(this->NumParticle() - 1);
	}

	/** \brief Insert a particle at a random position using the supplied RNG. */
	void Insert(OtherCharacteristics cha, RandomGenerator & gen)
	{
		GeometryVector a(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			a.x[i] = gen.RandomDouble();
		this->Insert(cha, a);
	}

	/** \brief Populate vicinity lattice up to distance Rc (squared cached). */
	void GetVicinityLattice(double Rc) const
	{
		size_t PrevVicinityLatticeNum = this->VicinityLattices.size();
		//prepare orthogonalized basis
		GeometryVector * ortho = new GeometryVector[this->Dimension];
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			ortho[i] = this->BasisVector[i];
			for (DimensionType j = 0; j<i; j++)
			{
				GeometryVector & e = ortho[j];
				ortho[i].MinusFrom(e*(ortho[i].Dot(e) / e.Dot(e)));
			}
		}

		this->_VicinityLatticeIteration(GeometryVector(this->Dimension), this->Dimension, GeometryVector(this->Dimension), ortho, Rc*Rc, this->VicinityLatticesRc2);

		delete[] ortho;

		if (this->SortedVicinityList)
			std::sort(this->VicinityLattices.begin() + PrevVicinityLatticeNum, this->VicinityLattices.end(), [](const VicinityLatticesMember & left, const VicinityLatticesMember & right) -> bool {return left.CartesianCoordinates.Modulus2()<right.CartesianCoordinates.Modulus2(); });

		this->VicinityLatticesRc2 = Rc*Rc;
	}

	/** \brief Debug-print vicinity lattice vectors. */
	void _PrintVicinityLattices(std::ostream & out) const
	{
		for (auto iter = this->VicinityLattices.begin(); iter<this->VicinityLattices.end(); iter++)
			out << (*iter).CartesianCoordinates << '\n';
	}

	/** \brief Pick a random particle index; optionally return its cell id. */
	size_t GetRandomParticle(RandomGenerator & gen, Index * SourceCell = nullptr) const//pick up a random particle, return it's pointer, if SourceCell is not null, also return the cell that the particle come from
	{
		size_t result;
		result = std::floor(gen.RandomDouble()*this->NumParticle());
		if (SourceCell != nullptr)
			*SourceCell = this->GetIndex(result);
		return result;
	}

	/** \brief Total number of particles stored. */
	size_t NumParticle(void) const
	{
		return this->ParticleRelatives.size();
	}

	/** \brief Access mutable per-particle characteristic. */
	virtual OtherCharacteristics & GetCharacteristics(size_t Num)
	{
		return this->ParticleCharacteristics[Num];
	}

	/** \brief Access const per-particle characteristic. */
	virtual const OtherCharacteristics & GetCharacteristics(size_t Num) const
	{
		return this->ParticleCharacteristics[Num];
	}

	/** \brief Relative (box) coordinates of a particle. */
	const GeometryVector & GetRelativeCoordinates(size_t Num) const
	{
		return this->ParticleRelatives[Num];
	}

	/** \brief Cartesian coordinates of a particle. */
	const GeometryVector & GetCartesianCoordinates(size_t Num) const
	{
		return this->ParticleCartesians[Num];
	}

	//return the length of the longest diagonal of the simulation box
	double MaxDiagonalLength(void) const
	{
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i]);
			else
				maxdist.MinusFrom(this->BasisVector[i]);
		}
		return std::sqrt(maxdist.Modulus2());
	}
	//allocate and generate internal cell lists, so that calling IterateThroughNeighbors with the same Rc is thread-safe
	void PrepareIterateThroughNeighbors(double Rc) const
	{
		//debug temp
		//std::cerr<<"IterateThroughNeighbors for Rc="<<Rc<<'\n';

		Rc = std::abs(Rc);
		//find the true Rc, which consideres the distance from a particle to a lattice point
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
			else
				maxdist.MinusFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
		}
		double newRc = Rc + std::sqrt(maxdist.Modulus2());
		double newRc2 = newRc*newRc;

		auto IterationEnd = this->VicinityLattices.end();
		//Get new VicinityLattice if necessary
		if (this->SortedVicinityList)
		{
			if ((newRc2)>this->VicinityLatticesRc2)
				this->GetVicinityLattice(newRc*1.1);
			IterationEnd = this->VicinityLattices.end();
		}
		else
		{
			auto ub = this->VicinityLatticePartitions.lower_bound(newRc2);
			auto lb = ub;
			if (this->VicinityLatticePartitions.size()>1)
				lb--;
			if ((newRc2)>this->VicinityLatticesRc2)
			{
				this->GetVicinityLattice(newRc);
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, this->VicinityLattices.size()));
				IterationEnd = this->VicinityLattices.end();
			}
			else if (ub->first == newRc2)
			{
				IterationEnd = ub->second + this->VicinityLattices.begin();
			}
			else if (ub->second - lb->second>10)
			{
				auto titer = std::partition(this->VicinityLattices.begin() + lb->second, this->VicinityLattices.begin() + ub->second,
					[&newRc2](const VicinityLatticesMember & a) ->bool
				{
					return a.CartesianCoordinates.Modulus2()<newRc2;
				});
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, titer - this->VicinityLattices.begin()));
				IterationEnd = titer;
			}
			else
				IterationEnd = ub->second + this->VicinityLattices.begin();
		}
	}

	/** \brief Iterate over neighbors within Rc of a point in relative coordinates.
	 *  \param RelativeCoordinate Center point in relative (box) coordinates.
	 *  \param Rc	Neighbor cutoff radius.
	 *  \param callback	A lambda function repeated for each neighbor with shift, lattice shift, periodic index shift, and source index.
	 *  \param pAbort	An optional flag to abort iteration early when set to true. */
	void IterateThroughNeighbors(const GeometryVector & RelativeCoordinate, double Rc, const std::function<void(const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)> callback, bool * pAbort = nullptr) const
	{
		//debug temp
		//std::cerr<<"IterateThroughNeighbors for Rc="<<Rc<<'\n';

		Rc = std::abs(Rc);
		double Rc2 = Rc*Rc;
		GeometryVector StartRelativeCoordinate = RelativeCoordinate;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			StartRelativeCoordinate.x[i] -= std::floor(StartRelativeCoordinate.x[i]);
			StartRelativeCoordinate.x[i] -= std::floor(StartRelativeCoordinate.x[i]);
		}
		//find the true Rc, which consideres the distance from a particle to a lattice point
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
			else
				maxdist.MinusFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
		}
		double newRc = Rc + std::sqrt(maxdist.Modulus2());
		double newRc2 = newRc*newRc;

		auto IterationEnd = this->VicinityLattices.end();
		//Get new VicinityLattice if necessary
		if (this->SortedVicinityList)
		{
			if ((newRc2)>this->VicinityLatticesRc2)
				this->GetVicinityLattice(newRc*1.1);
			IterationEnd = this->VicinityLattices.end();
		}
		else
		{
			auto ub = this->VicinityLatticePartitions.lower_bound(newRc2);
			auto lb = ub;
			if (this->VicinityLatticePartitions.size()>1)
				lb--;
			if ((newRc2)>this->VicinityLatticesRc2)
			{
				this->GetVicinityLattice(newRc);
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, this->VicinityLattices.size()));
				IterationEnd = this->VicinityLattices.end();
			}
			else if (ub->first == newRc2)
			{
				IterationEnd = ub->second + this->VicinityLattices.begin();
			}
			else if (ub->second - lb->second>10)
			{
				auto titer = std::partition(this->VicinityLattices.begin() + lb->second, this->VicinityLattices.begin() + ub->second,
					[&newRc2](const VicinityLatticesMember & a) ->bool
				{
					return a.CartesianCoordinates.Modulus2()<newRc2;
				});
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, titer - this->VicinityLattices.begin()));
				IterationEnd = titer;
			}
			else
				IterationEnd = ub->second + this->VicinityLattices.begin();
		}

		//Get cartesian coordinate of start point
		GeometryVector StartCartesianCoordinate(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			StartCartesianCoordinate.AddFrom(this->BasisVector[i] * StartRelativeCoordinate.x[i]);

		//Get the (dimensional) index of the cell that the start point is in
		signed long index[::MaxDimension] = {0};
		this->GetFullIndex(index, StartRelativeCoordinate);
		//iterate through cells
		for (auto iter = this->VicinityLattices.begin(); iter != IterationEnd; iter++)
		{
			//the vicinity lattice list can be larger than what's needed
			if (this->SortedVicinityList)
			{
				if (iter->CartesianCoordinates.Modulus2()>newRc2)
					break;
			}
			else
			{
				if (iter->CartesianCoordinates.Modulus2()>newRc2)
					continue;
			}


			signed long newindex[::MaxDimension]={0}; //the dimensional index for the cell
			signed long newshift[::MaxDimension]={0}; //the difference of periodic cell of new particle ant current particle, in terms of relative coordinate
			for (DimensionType i = 0; i<this->Dimension; i++)
			{
				signed long n = index[i] + iter->LatticeShift[i];

				if (n >= 0)
				{
					newindex[i] = (n) % this->CellRank[i];
					newshift[i] = (n) / this->CellRank[i];
				}
				else
				{
					newshift[i] = (-1) - (((-1) - n) / this->CellRank[i]);
					newindex[i] = this->CellRank[i] - 1 - (((-1) - n) % this->CellRank[i]);
				}
			}

			const std::vector<size_t > & thiscell = this->Cells[this->ConvertIndex(newindex)];

			if (thiscell.size() != 0)
			{
				GeometryVector LatticeShift(this->Dimension);
				for (int i = 0; i<this->Dimension; i++)
					LatticeShift.AddFrom(static_cast<double>(newshift[i])*this->BasisVector[i]);
				//iterate through Particles
				for (auto iter2 = thiscell.begin(); iter2 != thiscell.end(); iter2++)
				{
					if (pAbort != nullptr)
						if (*pAbort == true)
							return;

					GeometryVector shift = this->ParticleCartesians[*iter2] - StartCartesianCoordinate;

					shift.AddFrom(LatticeShift);

					if (shift.Modulus2()<Rc2)
						callback(shift, LatticeShift, newshift, *iter2);
				}
			}
		}
	}
	void MoveParticle(size_t Num, const GeometryVector & newRelativeCoordinate)
	{
		//particle * a=&this->Particles[Num];
		Index fromIndex;
		fromIndex = this->GetIndex(Num);
		std::vector<size_t > & currentCell = this->Cells[fromIndex];
		std::vector<size_t >::iterator itera = currentCell.begin();
		for (;; itera++)
		{
			assert(itera != currentCell.end());
			if (*itera == Num)
				break;
		}

		this->ParticleRelatives[Num] = newRelativeCoordinate;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			//while(this->ParticleRelatives[Num].x[i]>=1)
			//	this->ParticleRelatives[Num].x[i]-=1;
			//while(this->ParticleRelatives[Num].x[i]<0)
			//	this->ParticleRelatives[Num].x[i]+=1;
			this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
			this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
		}
		this->UpdateCartesianCoordinates(Num);
		this->UpdateCellList(itera, currentCell);
	}
	void DeleteParticle(size_t Num)
	{
		//particle * a=&this->Particles[Num];
		Index afromIndex = this->GetIndex(Num);
		std::vector<size_t > & acurrentCell = this->Cells[afromIndex];
		std::vector<size_t >::iterator itera = acurrentCell.begin();
		for (;; itera++)
		{
			assert(itera != acurrentCell.end());
			if (*itera == Num)
				break;
		}
		if (Num == this->NumParticle() - 1)
		{
			this->ParticleCartesians.pop_back();
			this->ParticleRelatives.pop_back();
			this->ParticleCharacteristics.pop_back();
			if (itera == acurrentCell.end() - 1)
				acurrentCell.pop_back();
			else
			{
				(*itera) = acurrentCell.back();
				acurrentCell.pop_back();
			}
			return;
		}

		//particle * back=&this->Particles.back();
		size_t backi = this->NumParticle() - 1;
		Index bfromIndex = this->GetIndex(backi);
		std::vector<size_t > & bcurrentCell = this->Cells[bfromIndex];
		std::vector<size_t >::iterator iterb = bcurrentCell.begin();
		for (;; iterb++)
		{
			assert(iterb != bcurrentCell.end());
			if (*iterb == (backi))
				break;
		}

		//(*a)=this->Particles.back();
		//this->Particles.pop_back();
		this->ParticleCartesians[Num] = this->ParticleCartesians.back();
		this->ParticleRelatives[Num] = this->ParticleRelatives.back();
		this->ParticleCharacteristics[Num] = this->ParticleCharacteristics.back();
		this->ParticleCartesians.pop_back();
		this->ParticleRelatives.pop_back();
		this->ParticleCharacteristics.pop_back();

		(*iterb) = Num;
		if (itera == acurrentCell.end() - 1)
			acurrentCell.pop_back();
		else
		{
			(*itera) = acurrentCell.back();
			acurrentCell.pop_back();
		}
	}
	GeometryVector GetBasisVector(::DimensionType i) const
	{
		return this->BasisVector[i];
	}
	
	GeometryVector GetReciprocalBasisVector(::DimensionType i) const
	{
		if (this->ReciprocalBasisVectorValid[i])
			return this->ReciprocalBasisVector[i];
		else
		{

			GeometryVector result = this->BasisVector[i];

			/////////////////////////////////////////////
			std::vector<GeometryVector> ortho;
			ortho.reserve(this->Dimension);
			for (DimensionType j = 0; j<this->Dimension; j++)
			{
				if (j == i) continue;
				GeometryVector temp = this->BasisVector[j];
				for (auto iter = ortho.begin(); iter != ortho.end(); iter++)
					temp.MinusFrom((*iter)*((temp.Dot(*iter)) / (iter->Modulus2())));
				ortho.push_back(temp);
			}
			////////////////////////////////////////////

			for (auto iter = ortho.begin(); iter != ortho.end(); iter++)
				result.MinusFrom((*iter)*((iter->Dot(result)) / (iter->Modulus2())));

			double product = result.Dot(this->BasisVector[i]);
			result.MultiplyFrom(2 * ::pi / product);
			this->ReciprocalBasisVector[i] = result;
			this->ReciprocalBasisVectorValid[i] = true;
			return result;
		}
	}

	//obtain indices of cells neighboring the RelativeCoordinates around radius Rc
	void GetNeighborList(const GeometryVector & RelativeCoordinates, double Rc, std::vector<int> & prt_index_list) const
	{
		prt_index_list.resize(0);
		DimensionType d = this->Dimension;
		std::vector<GeometryVector> reciprocalBasis;
		for (int i = 0; i < d; i++) {
			if (this->ReciprocalBasisVectorValid[i]) {//reciprocal vectors are already computed
				reciprocalBasis.push_back(this->ReciprocalBasisVector[i]);
			}
			else {//reciprocal vectors are yet to computed
				this->GetReciprocalBasisVector(i);
				this->ReciprocalBasisVector[i] = ((double)this->CellRank[i] / (2.0*::pi)) * ReciprocalBasisVector[i];
				reciprocalBasis.push_back(this->ReciprocalBasisVector[i]);
			}
		}
		
		//the cell within which the given point is located
		//the range of cells to span the sphere of radius Rc
		//Cell[i].x[0] = min. of the cell index in i-th relative coordinates
		//Cell[i].x[1] = max. of ...
		std::vector<GeometryVector> Cell;
		double center = 0;
		double change = 0;
		for (DimensionType i = 0; i < d; i++) {
			center = RelativeCoordinates.x[i] * (double)CellRank[i];
			change = Rc*std::sqrt(reciprocalBasis[i].Modulus2());
			double min = floor(center - change);
			double max = floor(center + change) + 1;
			Cell.push_back(GeometryVector(min, max));
		}

		if (this->NumParticle() > 0) {
			signed long *Indices = new signed long[d];
			for (DimensionType i = 0; i < d; i++) {
				Indices[i] = (long)Cell[i].x[0];
			}

			do {
				signed long *Indices_temp = new signed long[d];
				//Apply the periodic boundary condition
				for (int i = 0; i < d; i++) {
					if (Indices[i] >= CellRank[i]) {
						Indices_temp[i] = Indices[i] - CellRank[i];
					}
					else if (Indices[i] < 0) {
						Indices_temp[i] = Indices[i] + CellRank[i];
					}
					else {
						Indices_temp[i] = Indices[i];
					}
				}

				//convert index in the form of a vector into a long index
				signed long long_index = this->ConvertIndex(Indices_temp);
				for (int i = 0; i < this->Cells[long_index].size(); i++) {
					prt_index_list.push_back(this->Cells[long_index].at(i));
				}

				//to the next cell
				Indices[0]++;
				for (int i = 0; i < d - 1; i++) {
					if (Indices[i] == (long)Cell[i].x[1]) {
						Indices[i] = (long)Cell[i].x[0];
						Indices[i + 1] ++;
					}
				}
			} while (Indices[d - 1] < (long)Cell[d - 1].x[1]);
		}
		

				
	}

	void RefineCellList(void)
	{
		this->ClearVicinityLatticeList();
		size_t NumCells = 1;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double length = std::sqrt(this->BasisVector[i].Modulus2());
			this->CellRank[i] = (signed long)(std::floor(length / this->ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			NumCells *= this->CellRank[i];
		}
		this->Cells.resize(NumCells);
		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->Cells[this->GetIndex(i)].push_back(i);
		}
	}
	double GetCellSize(void)
	{
		return this->ApproximateCellSize;
	}
	void SetSortedCellList(bool SortedCellList)
	{
		if (SortedCellList != this->SortedVicinityList)
			this->ClearVicinityLatticeList();
		this->SortedVicinityList = SortedCellList;
	}
	void SetCellSize(double CellSize)
	{
		this->ApproximateCellSize = CellSize;
		bool NeedRefine = false;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double celllength = std::sqrt(this->BasisVector[i].Modulus2()) / this->CellRank[i];
			if (celllength > 2 * this->ApproximateCellSize)
				NeedRefine = true;
			else if (celllength < 0.5 * this->ApproximateCellSize && this->CellRank[i] > 1)
				NeedRefine = true;
		}
		if (NeedRefine)
			this->RefineCellList();
	}
	void ChangeBasisVector(const GeometryVector NewBasisVectors[])
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->BasisVector[i] = NewBasisVectors[i];
		this->ClearVicinityLatticeList();
		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->UpdateCartesianCoordinates(i);
		}

		bool NeedRefine = false;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double celllength = std::sqrt(this->BasisVector[i].Modulus2()) / this->CellRank[i];
			if (celllength > 2 * this->ApproximateCellSize)
				NeedRefine = true;
			else if (celllength < 0.5 * this->ApproximateCellSize && this->CellRank[i] > 1)
				NeedRefine = true;
		}
		if (NeedRefine)
			this->RefineCellList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
		this->vol = -1;
		this->vol = this->PeriodicVolume();
	}
	void RemoveParticles(void)//remove all Particles
	{
		//this->Particles.clear();
		this->ParticleCartesians.clear();
		this->ParticleCharacteristics.clear();
		this->ParticleRelatives.clear();

		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();
	}
	virtual void Rescale(double Factor)
	{
		this->ApproximateCellSize *= Factor;
		std::vector<GeometryVector> basis;
		for (DimensionType i = 0; i<this->Dimension; i++)
			basis.push_back(this->BasisVector[i] * Factor);
		this->ChangeBasisVector(&basis[0]);
		/*reset volume*/
		this->vol = -1.0;
	}
	virtual void Resize(double NewVolume)
	{
		double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
		this->Rescale(factor);
	}
	GeometryVector CartesianCoord2RelativeCoord(const GeometryVector & CartesianCoord) const
	{
		GeometryVector result(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			result.x[i] = CartesianCoord.Dot(this->GetReciprocalBasisVector(i));

		return result*(0.5 / pi);
	}
	GeometryVector RelativeCoord2CartesianCoord(const GeometryVector & RelativeCoord) const
	{
		GeometryVector result(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			result.AddFrom((this->BasisVector[i])*(RelativeCoord.x[i]));

		return result;
	}
	void RelativeCoordToMinimumImage(GeometryVector & a) const
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			a.x[i] -= std::floor(a.x[i]+0.5);
	}
	void WriteBinary(std::ostream & ofile) const
	{
		ofile.write((char *)(&this->Dimension), sizeof(this->Dimension));
		ofile.write((char *)(&ApproximateCellSize), sizeof(ApproximateCellSize));
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->BasisVector[i].WriteBinary(ofile, this->Dimension);
		int nbrp = this->NumParticle();
		ofile.write((char *)(&nbrp), sizeof(nbrp));
		for (int i = 0; i<nbrp; i++)
		{
			this->ParticleRelatives[i].WriteBinary(ofile, this->Dimension);
			ofile.write((char *)(&this->ParticleCharacteristics[i]), sizeof(OtherCharacteristics));
		}
	}
	// there is no such thing as void ReadBinary(std::istream & ifile)! Use the following constructor instead!
	PeriodicCellList(std::istream & ifile, bool UseSortedList = false)
	{
		ifile.read((char *)(&this->Dimension), sizeof(this->Dimension));
		ifile.read((char *)(&ApproximateCellSize), sizeof(ApproximateCellSize));
		Index CellCount = 1;//need to calculate overall number of cells
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			this->BasisVector[i].ReadBinary(ifile, this->Dimension);
			GeometryVector & vec = this->BasisVector[i];
			this->CellRank[i] = (signed long)(std::floor(std::sqrt(vec.Modulus2()) / ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			//error check
			if (this->CellRank[i] <= 0)
			{
				std::cerr << "error in Configuration constructor: invalid cell rank. basis vector length squared=" << vec.Modulus2() << ", CellSize=" << ApproximateCellSize << '\n';
				assert(false);
			}
			CellCount *= this->CellRank[i];
		}
		for (DimensionType i = this->Dimension; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = 0;
			this->BasisVector[i] = GeometryVector();
		}
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = ApproximateCellSize;
		this->SortedVicinityList = UseSortedList;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		int nbrp;
		ifile.read((char *)(&nbrp), sizeof(nbrp));
		for (int i = 0; i<nbrp; i++)
		{
			GeometryVector temp;
			temp.ReadBinary(ifile, this->Dimension);
			OtherCharacteristics otemp;
			ifile.read((char *)(&otemp), sizeof(OtherCharacteristics));
			this->Insert(otemp, temp);
		}
		this->vol = ::Volume(this->BasisVector, this->Dimension);
	}

	//rotate the current structure so that
	//the 1st basis vector is (c11, 0, 0, ...)
	//the 2nd basis vector is (c21, c22, 0, ...)
	// ......
	//AND the basis vector length l1>=l2>=l3...
	void StandardizeOrientation(void)
	{
		::StandardizeOrientation(this->BasisVector, this->Dimension);
		this->ChangeBasisVector(this->BasisVector);
	}
	//PeriodicCellList< of anything > can be converted to PeriocicCellList<Empty>
	operator PeriodicCellList<Empty> const () const
	{
		return PeriodicCellList<Empty>(*this, Empty());
	}

	//return the distance from fromRelative to the closest particle
	double NearestParticleDistance(const GeometryVector & fromRelative)
	{
		double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension())*1.5;
		std::vector<GeometryVector> neighbors;
		double l = TypicalLength;
		while (neighbors.size() < 1)
		{
			neighbors.clear();
			(*this).IterateThroughNeighbors(fromRelative, l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 1.5;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin() + 1, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2() < right.Modulus2(); });
		return std::sqrt(neighbors[0].Modulus2());
	}
	//return the distance from particle No. i to the closest particle
	double NearestParticleDistance(size_t i)
	{
		double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension())*1.5;
		std::vector<GeometryVector> neighbors;
		double l = TypicalLength;
		while (neighbors.size() < 2)
		{
			neighbors.clear();
			(*this).IterateThroughNeighbors((*this).GetRelativeCoordinates(i), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 1.5;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2() < right.Modulus2(); });
		return std::sqrt(neighbors[1].Modulus2());
	}

	//Least image displacement : input vector is in relative coordinates
	GeometryVector displacement_periodic_relative (int i, const GeometryVector & r0_Relative) const {
		if (i < this->NumParticle()) {
			GeometryVector Displacement = ParticleRelatives[i] - r0_Relative;
			for (int k = 0; k < this->Dimension; k++) {
				if (Displacement.x[k] < -0.5) { Displacement.x[k] += 1.0; }
				else if (Displacement.x[k] >0.5) { Displacement.x[k] -= 1.0; }
			}
			Displacement = this->RelativeCoord2CartesianCoord(Displacement);
			return Displacement;
		}
		else {
			std::cout << "PeriodicCellList::displacement_periodic - i is too large \n";
			return GeometryVector(r0_Relative.Dimension);
		}
	}
	//input vector is in Cartesian coordinates
	GeometryVector displacement_periodic(int i, const GeometryVector & r0) const {
		if (i < this->NumParticle()) {
			GeometryVector Displacement = ParticleCartesians[i] - r0;
			this->CartesianCoord2RelativeCoord(Displacement);
			for (int k = 0; k < this->Dimension; k++) {
				if (Displacement.x[k] < -0.5) { Displacement.x[k] += 1.0; }
				else if (Displacement.x[k] >0.5) { Displacement.x[k] -= 1.0; }
			}
			Displacement = this->RelativeCoord2CartesianCoord(Displacement);
			return Displacement;
		}
		else {
			std::cout << "PeriodicCellList::displacement_periodic - i is too large \n";
			return GeometryVector(r0.Dimension);
		}
	}


	GeometryVector GetRandomVector() {
		RandomGenerator gen;
		GeometryVector result;
		for (int i = 0; i < MaxDimension; i++) {
			result = result + gen.RandomDouble()*this->BasisVector[i];
		}
		return result;
	}

	//return squared distance btw particle i and j
	double distance2(int i, int j) {
		if (i < this->NumParticle() && j < this->NumParticle()) {
			GeometryVector dr = this->ParticleRelatives[i] - this->ParticleRelatives[j];
			for (int k = 0; k < this->Dimension; k++) {
				if (dr.x[k] < -0.5) { dr.x[k] += 1.0; }
				else if (dr.x[k] >0.5) { dr.x[k] -= 1.0; }
			}
			GeometryVector temp = this->RelativeCoord2CartesianCoord(dr);
			return temp.Dot(temp);
		}
		else {
			std::cout << "wrong index\n";
			return 0;
		}
	}
	//the squared distance between i-th particle and test particle		
	double distance2(int i, GeometryVector test) {
		if (i < this->NumParticle()) {
			GeometryVector dr = ParticleCartesians[i] - test;
			this->CartesianCoord2RelativeCoord(dr);
			for (int k = 0; k < this->Dimension; k++) {
				if (dr.x[k] < -0.5) { dr.x[k] += 1.0; }
				else if (dr.x[k] >0.5) { dr.x[k] -= 1.0; }
			}
			GeometryVector temp = this->RelativeCoord2CartesianCoord(dr);
			return temp.Modulus2();
		}
		else {
			std::cout << "PeriodicCellList::distance2 - i is too large \n";
			return 0;
		}
	}

	double distance(int i, int j) {
		return sqrt(this->distance2(i, j));
	}

};
//overlap all the configurations in that std::vector
//require that all the configurations in the std::vector have the same basis vectors
template<typename a>PeriodicCellList<a> CombinedList(const std::vector<PeriodicCellList<a> > & components)
{
	if(components.size()==0)
		return PeriodicCellList<a>();
	PeriodicCellList<a> result = components[0];
	for(size_t i=1; i<components.size(); i++)
	{
		const PeriodicCellList<a> & list = components[i];
		for(size_t j=0; j<list.NumParticle(); j++)
			result.Insert(list.GetCharacteristics(j), list.GetRelativeCoordinates(j));
	}
	return result;
}
//config is a CombinedList of m PeriodicCellLists, extract the nth list
template<typename a>PeriodicCellList<a> SeparateList(const PeriodicCellList<a> & config, size_t m, size_t n)
{
	size_t nt=config.NumParticle();
	if(nt%m!=0)
	{
		std::cerr<<"Error in SeparateList : the total number of particles is not multiple of the number of configurations!\n";
		assert(false);
	}
	size_t np=nt/m;
	size_t start=n*np;
	size_t end=start+np;
	PeriodicCellList<a> result = config;
	result.RemoveParticles();
	for(size_t i=start; i<end; i++)
	{
		result.Insert(config.GetCharacteristics(i), config.GetRelativeCoordinates(i));
	}
	return result;
}

struct AtomName
{
	char name[4];
	AtomName()
	{
		this->name[0]='\0';
	}
	AtomName(const char * n)
	{
		std::strncpy(this->name, n, 3);
		this->name[3]='\0';
	}
};
typedef PeriodicCellList<AtomName> Configuration;

void Output(std::string prefix, const Configuration & List);
void Output(std::ostream & out, const Configuration & List);
inline void WriteBinary(std::ostream & out, const Configuration & List)
{
	List.WriteBinary(out);
}
inline void WriteBinary(std::string prefix, const Configuration & List)
{
	prefix+=std::string(".Config");
	std::fstream ofile(prefix.c_str(), std::fstream::out | std::fstream::binary);
	WriteBinary(ofile, List);
}
inline Configuration ReadBinary(std::istream & ifile)
{
	return Configuration(ifile);
}
inline Configuration ReadBinary(std::string prefix)
{
	prefix+=std::string(".Config");
	std::fstream ifile(prefix.c_str(), std::fstream::in | std::fstream::binary);
	return ReadBinary(ifile);
}

Configuration ReadPos(std::string prefix);
Configuration ReadPos(std::istream & ifile);
Configuration ReadHoomdXml(std::istream & ifile);
Configuration ReadHoomdXml(std::string prefix);
Configuration ReadStealthOutput(std::istream & file, double SideLength);
Configuration ReadStealthOutput(std::istream & file);
Configuration MultiplicateStructure(const Configuration & src, size_t NumberInEachSide);
Configuration MultiplicateStructure(const Configuration & src, const std::vector<size_t> & NumberInEachSide);
Configuration ReadGro(std::istream & file);
Configuration GetUnitCubicBox(DimensionType d, double CellSize=0.1);

Configuration ReadCoordinate(std::istream & file);
Configuration ReadCoordinate(std::istream & file, double SideLength);
Configuration ReadCartesianCoordinate(std::istream & file);
Configuration ReadCartesianCoordinate(std::istream & file, double SideLength);



// a pack of multiple configurations, recorded in files Prefix.ConfigPack and Prefix.ConfigPackIndex
// .ConfigPack file contains configurations in binary format, one by one
// .ConfigPackIndex has the format:
// long long Number of Config, long long location of config0, long long location of config1, ...
class ConfigurationPack
{
private:
	std::string IndexName, PackName;
	long long NConfig;
public:
	ConfigurationPack()
	{
		this->NConfig=0;
	}
	ConfigurationPack(const std::string & prefix)
	{
		this->Open(prefix);
	}
	void Open(const std::string & prefix);
	long long NumConfig(void);
	void AddConfig(const Configuration & c);
	Configuration GetConfig(long long i);

	//discard all configurations that are already in the file
	void Clear(void);
};


//polydisperse sphere packing
class SpherePacking : public PeriodicCellList<double>
{
public:
	SpherePacking() : PeriodicCellList()//, MaxRadius(0.0)
	{
	}
	SpherePacking(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList=false);
	template<typename add> SpherePacking(const PeriodicCellList<add> & source, const double & radius) : PeriodicCellList(source, radius)
	{
		this->MaxRadius=radius;
	}
	virtual double & GetCharacteristics(size_t Num)
	{
		//this->MaxRadius = 0.0;
		return this->ParticleCharacteristics[Num];
	}
	virtual const double & GetCharacteristics(size_t Num) const
	{
		return this->ParticleCharacteristics[Num];
	}

	SpherePacking(std::istream & ifile, bool UseSortedList=false) : PeriodicCellList<double>(ifile, UseSortedList)
	{
		this->MaxRadius=0.0;
		for(size_t i=0; i<this->NumParticle(); i++)
			this->MaxRadius=std::max(this->MaxRadius, this->GetCharacteristics(i));
	}
	void UpdateMaxRadius(void) const
	{
		if (this->MaxRadius == 0.0)
		{
			for (size_t i = 0; i < this->NumParticle(); i++)
				this->MaxRadius = std::max(this->MaxRadius, this->ParticleCharacteristics[i]);
		}
	}
	double GetMaxRadius(void) const
	{
		UpdateMaxRadius();
		return this->MaxRadius;
	}
	virtual void Insert(double radius, const GeometryVector & RelativeCoordinate);

	//check if the voxel is COMPLETELY occupied
	//no need to distinguish partial availability and full availability
	//relative coord and halfsize 
	//Radius : cartesian radius of the upcoming sphere
	bool CheckVoxelOverlap(const GeometryVector & cr, double halfsize, double Radius) const;

	//returns 0 if voxel is completely empty
	//2 if voxel is fully occupied
	//1 otherwise
	//relative coord and halfsize 
	//Radius : cartesian radius of the upcoming sphere
	int CheckVoxelOverlap2(const GeometryVector & cr, double halfsize, double Radius) const;


	//relative coord
	//Radius : cartesian radius of the upcoming sphere
	bool CheckOverlap(const GeometryVector & cr, double Radius) const;

	//check if there are overlapping spheres
	/* @return false = no overlapping. */
	bool CheckOverlap(void) const
	{
		UpdateMaxRadius();
		for (size_t i = 0; i < this->NumParticle(); i++)
		{
			bool Overlap = false;
			const std::vector<double> * pSphereRadii = &this->ParticleCharacteristics;
			double Radius = this->GetCharacteristics(i);
			this->IterateThroughNeighbors(this->GetRelativeCoordinates(i), this->GetCharacteristics(i) + MaxRadius, [&Overlap, &Radius, &pSphereRadii, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
			{
				if (Sourceparticle != i || LatticeShift.Modulus2() != 0)
				{
					double OverlapRadius = Radius + pSphereRadii->at(Sourceparticle);
					if (shift.Modulus2() < OverlapRadius*OverlapRadius)
						Overlap = true;
				}
			}, &Overlap);
			if (Overlap)
				return true;
		}
		return false;
	}

	double PackingFraction(void) const
	{
		double v1=0.0;
		DimensionType d=(DimensionType) this->GetDimension();
		for(size_t i=0; i<this->NumParticle(); i++)
			v1+= ::HyperSphere_Volume(d, this->GetCharacteristics(i));
		return v1/this->PeriodicVolume();
	}

	//-----rescale system box and particle sizes....
	virtual void Rescale(double Factor)
	{
		this->ApproximateCellSize *= Factor;
		this->MaxRadius *= Factor;
		std::vector<GeometryVector> basis;
		for (DimensionType i = 0; i<this->Dimension; i++)
			basis.push_back(this->BasisVector[i] * Factor);
		this->ChangeBasisVector(&basis[0]);
		
		for (size_t i = 0; i < this->NumParticle(); i++) {
			this->ParticleCharacteristics[i] *= Factor;
		}
		/* inform that the volume is updated. */
		this->vol = -1.0;
	}
	virtual void Resize(double NewVolume)
	{
		double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
		this->Rescale(factor);
	}
	//Rescale particle sizes only
	void Rescale_prtOnly(double factor) {
		for (size_t i = 0; i < this->NumParticle(); i++) {
			this->ParticleCharacteristics[i] *= factor;
		}
		this->MaxRadius *= factor;
	}

	void UpdateCharacteristics(double newRadius, int index) {
		this->ParticleCharacteristics[index] = newRadius;
		this->MaxRadius = std::max(this->MaxRadius, newRadius);
	}
	
	// Newly defined for First-passage time simulations

		/** Check whether a test point is inside a particle of this Dispersion class. 
	 * @param[in] RelativeCoordinate 	The position of a test point in the relative coordinates.	 
	 * @return Particle ID, if the point is inside of a particle. -1, otherwise. */
	int IsInsideParticle(const GeometryVector & RelativeCoordinate) const;

	/** Find the three nearest neighbors. 
	 * @param[in] RelativeCoordinate	The position of a test point in the relative coordinates.
	 * @param[out] list	A list of some nearest interfaces. list[i] = GeometryVector(distance, side ID, particle ID)
	 * @param[out] x0_Relative	The point on the interface closest to RelativeCoordinate. 	
	 * @param[out] n2host	The unit outward normal vector at x0. 	*/
	void FindNearestNeighbors (const GeometryVector & RelativeCoordinate, std::vector<GeometryVector> & list, GeometryVector & x0_Relative, GeometryVector & n2host) const;

};


SpherePacking ReadTJOutput(std::istream & inFile);
SpherePacking ReadTJOutput(const std::string & prefix);


void WriteTJOutput(std::fstream & ofile, const SpherePacking & pk);
void WriteTJOutput(const std::string & prefix, const SpherePacking & pk);

/*	 a SpherePacking counterpart of ConfigurationPack
*.SpherePack files and .SpherePackIndex
*/
class SpherePackingPack 
{
private:
	std::string IndexName, PackName;
	long long NConfig;
public:
	SpherePackingPack() { this->NConfig = 0; }
	SpherePackingPack(const std::string & prefix) {
		this->Open(prefix);
	}
	void Open(const std::string &prefix);
	long long NumConfig(void);
	void AddConfig(const SpherePacking &c);
	SpherePacking GetConfig(long long i);
	void Clear(void);
};

//VoxelType should be either char or bool. Char is faster, while bool can work with smaller memory.
typedef char VoxelType;
//pixelizing the configuration to MeshSide pixes in each side
void DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> & result);
inline std::vector<VoxelType> DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide)
{
	std::vector<VoxelType> result;
	DigitizeConfiguration(c, radius, MeshSide, result);
	return result;
}
//calculate volume fraction if each particle is replaced with a hypersphere of certain radius
//do this by pixelizing the configuration to MeshSide pixes in each side
double Volume(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> * pbuffer = nullptr);


//return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
double MinConnectDiameter(const Configuration & c, size_t a, size_t b);
//return the minimum D such that if each particle is replaced with a sphere with diameter D, 
//particle n will connect to its periodic image in direction dir.
double PercolationDiameter(const Configuration & c, size_t n, DimensionType dir);


//Replace each particle in c with a hypersphere with radius,
//what's the cluster size InitParticle is in?
//FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, 
//which means the cluster size is actually infinite rather than the returned value.
size_t ClusterSize(const Configuration & c, size_t InitParticle, double radius, bool & FoundPeriodicImages);



#include <complex>

//This class initializes by assigining each particle in a PeriodicCellList to its own cluster
//Then, as the user call member function Proceed(double r), particles within distance r are merged to a single cluster
//The user can then query various member objects of this class for various aspects of clusters:
//member clusters is a list of all clusters. Each cluster is represented as a list of its members' index
//affiliation[i] is the cluster particle i currently belong to
class ClusterMerger
{
protected:
	typedef std::list<size_t> cluster;
private:
	//forbid copy
	ClusterMerger(const ClusterMerger & src);
	struct Bond
	{
		double l2;
		size_t i, j;
	};
	std::vector<Bond> vBond;
	std::vector<Bond>::iterator next;
	size_t _N;
	std::list<cluster> _clusters;
	std::vector<std::list<cluster>::iterator> _affiliation;
	size_t _LargestClusterSize;
public:
	const size_t & N=_N;

//clusters is a list of all clusters. Each cluster is represented as a list of its members' index
	const std::list<cluster> & clusters=_clusters;

//affiliation[i] is the cluster particle i currently belong to
	const std::vector<std::list<cluster>::iterator> & affiliation=_affiliation;
	const size_t & LargestClusterSize=_LargestClusterSize;


public:
	ClusterMerger(const PeriodicCellList<Empty> & c, double Rmax)
	{
		_N = c.NumParticle();
		for (size_t i = 0; i < _N; i++)
		{
			//at the beginning, each particle consist its own cluster
			cluster temp;
			temp.insert(temp.begin(), i);
			_clusters.insert(_clusters.begin(), temp);
			_affiliation.push_back(_clusters.begin());
			std::vector<Bond> & vb = vBond;
			//generate a vector of bonds for future computation
			c.IterateThroughNeighbors(c.GetRelativeCoordinates(i), Rmax, [&vb, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				if (i < Sourceparticle)
				{
					Bond temp;
					temp.i = i;
					temp.j = Sourceparticle;
					temp.l2 = shift.Modulus2();
					vb.push_back(temp);
				}
			});
		}
		std::sort(vBond.begin(), vBond.end(), [](const Bond & left, const Bond & right) ->bool
		{
			return left.l2 < right.l2;
		});
		_LargestClusterSize = 1;
		next = vBond.begin();
	}
	//merge _clusters until distances of any cluster is larger than r
	void Proceed(double r)
	{
		double r2 = r*r;
		while (next != vBond.end() && next->l2 < r2)
		{
			std::list<cluster>::iterator pi = _affiliation[next->i], pj = _affiliation[next->j];
			if (pi != pj)
			{
				//merge the two _clusters
				cluster & cj = *_affiliation[next->j];
				cluster & ci = *_affiliation[next->i];
				//1. change the _affiliation of all particles in the same cluster as j
				for (auto iter = cj.begin(); iter != cj.end(); ++iter)
				{
					_affiliation[*iter] = _affiliation[next->i];
				}
				//2. merge the two _clusters
				ci.splice(ci.begin(), cj);
				_clusters.erase(pj);
				_LargestClusterSize = std::max(_LargestClusterSize, ci.size());
			}
			next++;
		}
	}
};


//uncomfortably added functions
//output Cartesian coordinates of particles with the dimension and basis vectors
void WriteConfiguration(const Configuration & conf, const std::string & OutputName);
void WriteConfiguration(const SpherePacking & conf, const std::string & OutputName);
template<typename T> void WriteConfiguration(const PeriodicCellList<T> & conf, const std::string & OutputName);
void ReadConfiguration(Configuration & conf, const std::string & InputName);
void ReadConfiguration(SpherePacking & conf, const std::string & InputName);


/** Find all reciprocal lattice vectors (wavevectors) of a given periodic simulation box.
 *  It is modified to deal with errors.
 *	@param tempList	A periodic simulation box.
 *	@param (CircularKMax, LinearKMax)	Find all reciprocal lattice vectors (excluding parity pairs) whose magnitude is smaller than CircularKMax but greater than CircularKMin, and then fina all integral multiples of the aforementioned vectors whose magnitudes are smaller than LinearKMax.
 *	@param SampleProbability	When <1, randomly choose wavevectors within CircularKMax by the prescribed probability.
 *	@param use_iteratethrough	True = use IterateThrough() function, which may surge memory use. 
 *	@param CircularKMin			0 by default.  
 */

std::vector<GeometryVector> GetKs(const PeriodicCellList<Empty> & tempList, double CircularKMax, double LinearKMax, double SampleProbability = 1.0, bool use_iteratethrough = true, double CircularKMin = 0.0);



#endif
