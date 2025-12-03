/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com
*	Date	: November 2025 */

/** \file MultiStealthyPotential.h
 * \brief Header file to declare a class to compute multi-stealthy collective coordinate potential + soft-core repulsions */

#ifndef MultiStealthyPotential_H
#define MultiStealthyPotential_H

/* in $(core) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
/* in $(pot) */
#include <Potential.h>
#include <CollectiveCoordinatePotential.h>

#include <algorithm>
#include <memory>

struct Constraint_Fourier
{
	GeometryVector k;//absolute k, not relative to basis vectors
	double V;
	double rho_real, rho_imag;
	Constraint_Fourier(const GeometryVector & k, double V): k(k), V(V)
	{
	}
};

/* Constraint energies in Fourier space for a given subset of particles */
struct EnergyConstraints
{
	std::vector<size_t> idx ;
	bool RhoValid = false;
	std::vector<Constraint_Fourier> constraints;

	/**
	 * @brief Add a Fourier-space constraint at a wavevector for a species
	 * 
	 * @param k 	A wavevector
	 * @param v 	Fourier-space energy, must be non-negative.
	 */
	void AddConstraint(const GeometryVector & k, double v){
		this->constraints.emplace_back(k, v);
	}

	/**
	 * @brief Print the defined constraints with the indices of species
	 * 
	 * @param ofile 
	 */
	void PrintConstraints(std::ostream & ofile = std::cout){
		// particle indices
		ofile << "indices:";
		for (auto i = idx.begin(); i != idx.end(); i++)
			ofile << "\t" << &i ;

		ofile <<"\n";
		// constraints
		for (size_t i =0; i<this->constraints.size(); i++){
			ofile << "@ k= " << this->constraints[i].k << ": v="<<this->constraints[i].V <<"\n";
		}
		ofile << "\n\n";
	}
};

/* collective-coordinate potential in a species. */
struct species_potential{
	std::unique_ptr<ShiftedCCPotential> CCPotential;
	std::vector<size_t> idx ;   // particle i in CCPotential -> pConfig  =  particle idx[i] in the mother configuration.

	species_potential (DimensionType dim): CCPotential(std::make_unique<ShiftedCCPotential>(dim)), idx() {}
	// non-copyable to preserve unique ownership of CCPotential
	species_potential(const species_potential&) = delete;
	species_potential& operator=(const species_potential&) = delete;
	species_potential(species_potential&&) = default;
	species_potential& operator=(species_potential&&) = default;

	/**
	 * @brief Print the defined constraints with the indices of species
	 * 
	 * @param ofile 
	 */
	void PrintConstraints(std::ostream & ofile = std::cout){
		// particle indices
		ofile << "indices:";
		for (auto i = idx.begin(); i != idx.end(); i++)
			ofile << "\t" << &i ;

		ofile <<"\n";
		// constraints
		for (size_t i =0; i<CCPotential->constraints.size(); i++){
			ofile << "@ k= " << CCPotential->constraints[i].k << ": v="<<CCPotential->constraints[i].V <<"\n";
		}
		ofile << "\n\n";
	}

};


struct spec_index{
	int species;
	size_t index;

	spec_index(int s, size_t idx): species(s), index(idx){}
};

class MultiStealthyPotential : public Potential
{
protected:
	std::unique_ptr<Configuration> pConfig;	//< packing of interest
	std::vector<double> radii; 		//< particle radii
	double rc;						//< the critical radius over which the contact potential become zero.
	bool WithRepulsion = true;
	bool SubconfigValid = false;
	//	bool RhoValid;
	//	double Volume;
	std::vector<spec_index> idx2species;
	std::vector<Configuration> subc;


public :
	double e0 = 1.0;	//< energy scale for contact potential
	double u0 = 1.0;	//< energy scale for collective-coordinate potential
	
	// constraints[i]: Fourier-space energy constraints for a species i.
	std::vector<species_potential> species_constraints;

	/// @brief Constructor
	/// @param dimension 
	MultiStealthyPotential(DimensionType dimension, double epsilon): Potential(dimension){
		this->pConfig=nullptr;
		this->rc = -1.0;
		this->ParallelNumThread=1;
		this->e0 = epsilon;
		this->WithRepulsion = (std::abs(this->e0) > 1e-10);
		if(this->WithRepulsion){
			std::cout<< "Soft-core repulsion is on\n";
		}
		else{
			std::cout<< "Soft-core repulsion is off\n";
		}
	}

	~MultiStealthyPotential() override = default;
	MultiStealthyPotential(const MultiStealthyPotential&) = delete;
	MultiStealthyPotential& operator=(const MultiStealthyPotential&) = delete;
	MultiStealthyPotential(MultiStealthyPotential&&) = default;
	MultiStealthyPotential& operator=(MultiStealthyPotential&&) = default;

	MultiStealthyPotential * clone(){}


	/* -------------------------------------
	 * 	Methods to define configuration and potential
	 * -------------------------------------*/

	/** @brief Set the number of species in the system. */
	/**
	 * @brief Set the number of species in this configuration.
	 * 
	 * @param num 
	 */
	void SetNumSpecies(const size_t num){
		for (size_t i = 0; i < num ; i++){
			species_constraints.emplace_back(this->Dimension);
		}
	}

	/// @brief Determine the species of a particle
	/// @param prt_idx 
	/// @param species 
	void SetSpeciesOfParticle(size_t prt_idx, size_t species);

	/** @brief Add a constraint to a species. */
	void AddConstraint(size_t species, const GeometryVector & k, const std::vector<double> V);
  
	/** @brief Update rho values for selected species. */
	void GetRho(size_t species);

	/** @brief Print constraints to an output stream. */
	void PrintConstraints(std::ostream & ofile = std::cout);

	/** @brief Set the working configuration. */
	virtual void SetConfiguration(const Configuration & Position);
	
	/** @brief Set particle radii explicitly. */
	void SetSecureRadii(const std::vector<double> & radii);

	/** @brief Set the configuration from a sphere packing. */
	void SetPacking(const SpherePacking & pck);

	/** @brief Rebuild species sub-configurations. */
	void UpdateSubconfigs();

	/** @brief Get the contact cutoff radius. */
	double GetRc(void) const {return this->rc;};

	/* -------------------------------------
	 * 	Implementations of this class.
	 * -------------------------------------*/

	/** @brief Compute total potential energy. */
	virtual double Energy();
	/** @brief Compute force on a particle. */
	virtual void Force(GeometryVector & result, size_t i);
	/** @brief Compute forces on all particles. */
	virtual void AllForce(std::vector<GeometryVector> & results);
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure);

};


#endif 
