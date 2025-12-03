#include "MultiStealthyPotential.h"

/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: November 2025 */

/** \file MultiStealthyPotential.cpp 
 * \brief c++ implementations of multi-stealthy potential */



/**
 * @brief Assign a particle to a species and refresh bookkeeping indices.
 * 
 * @param prt_idx Global particle index in the full configuration.
 * @param species Target species id.
 */
void MultiStealthyPotential::SetSpeciesOfParticle(size_t prt_idx, size_t species)
{
	this->species_constraints[species].idx.emplace_back(prt_idx);

	if(this->idx2species[prt_idx].species!= species){
		/* species is replaced! */
		if (this->idx2species[prt_idx].species > -1){ // the species of a particle was designated previously.
			std::cout << "change species of particle " << prt_idx << " from " << this->idx2species[prt_idx].species << " to " << species << "\n";
		}
	    this->species_constraints[this->idx2species[prt_idx].species].idx.erase(
			std::remove(
				this->species_constraints[this->idx2species[prt_idx].species].idx.begin(), 
				this->species_constraints[this->idx2species[prt_idx].species].idx.end(), 
				prt_idx
			), 
			this->species_constraints[this->idx2species[prt_idx].species].idx.end()); // remove the previously defined species of particle `prt_idx`.
	}
	this->idx2species[prt_idx].species = species;
	this->idx2species[prt_idx].index = this->species_constraints[species].idx.size()-1;
	
	this->SubconfigValid = false;
}



/**
 * @brief Append a Fourier-space constraint to a given species.
 * 
 * @param species Species id to target.
 * @param k Wavevector for the constraint.
 * @param V List of constraint coefficients.
 */
void MultiStealthyPotential::AddConstraint(size_t species, const GeometryVector &k, std::vector<double> V)
{
	this->species_constraints[species].CCPotential->AddConstraint(k, V);
}

/**
 * @brief Compute the structure factor contributions for a species.
 * 
 * @param s Species id to evaluate.
 */
void MultiStealthyPotential::GetRho(size_t s)
{
	assert(this->pConfig!=nullptr);
	assert(s < this->species_constraints.size());
	this->species_constraints[s].CCPotential->GetRho();

	// double Volume=ptr_pot->pConfig->PeriodicVolume();

	// size_t end = this->species_constraints[s].constraints.size();
	// size_t num_in_species = this->species_constraints[s].idx.size();
	// #pragma omp parallel for num_threads(this->ParallelNumThread)
	// for (size_t i = 0; i < end; i++){
	// 	/* At a wavevector for the constraint*/
	// 	double tr=0, ti=0;//real and imag parts of result
		
	// 	for(size_t j = 0; j < num_in_species; j++)
	// 	{
	// 		/* repeat over particles in the same species */
	// 		size_t p_idx = this->species_constraints[s].idx[j];
	// 		GeometryVector x=pConfig->GetCartesianCoordinates(p_idx);
	// 		tr+=std::cos(this->species_constraints[s].constraints[j].k.Dot(x));
	// 		ti+=std::sin(this->species_constraints[s].constraints[j].k.Dot(x));
	// 	}
	// 	this->species_constraints[s].constraints[i].rho_real=tr;
	// 	this->species_constraints[s].constraints[i].rho_imag=ti;
	// }
	// this->species_constraints[s].RhoValid = true;
}

/**
 * @brief Print all species constraints and indices.
 * 
 * @param ofile Destination output stream.
 */
void MultiStealthyPotential::PrintConstraints(std::ostream &ofile)
{
	for (size_t s = 0; s < this->species_constraints.size(); s ++ ){
		ofile << "speciecs "<< s << "\n";
		ofile << "\t (particle index, radius)";
		for (size_t i = 0; i < this->species_constraints[s].idx.size(); i++){
			size_t idx = this->species_constraints[s].idx[i];
			ofile << "\t(" << idx << ",\t" << this->radii[idx] <<")\n"; //this->pConfig->GetCharacteristics(idx) <<")\n";
		}

		ofile << "\n";
		for (size_t i =0; i< this->species_constraints[s].CCPotential->constraints.size(); i++){
			ofile << "@ k= " << this->species_constraints[s].CCPotential->constraints[i].k << ": v=" << this->species_constraints[s].CCPotential->constraints[i].V <<"\n";
		}
		ofile << "----- \n\n";
	}
}

/**
 * @brief Set the working configuration from particle positions.
 * @param Position Full particle configuration.
 */
void MultiStealthyPotential::SetConfiguration(const Configuration &Position)
{
	if (this->rc < 0.0){
		// an empty configuration
		this->idx2species.resize(Position.NumParticle(), spec_index(-1,0));
		std::cout << "species of all particles are set to be -1 (unspecified)\n";
		std::cout << "the critical radius of the contact potential is set to be "<< this->rc <<std::endl;

		std::cout << "No input for particle radii\n";
		std::cout << "Particle radii are set to be 0.5 by default.\n";	
		this->rc = 1.0;
		if (this->radii.size() != Position.NumParticle()){
			this->radii.resize(Position.NumParticle(), -1.0);
		}
	}
	
	pConfig = std::make_unique<Configuration>(Position, "a");
	pConfig->PrepareIterateThroughNeighbors(this->rc);


	this-> SubconfigValid = false;
	this-> UpdateSubconfigs();
	for (size_t s = 0; s < this->species_constraints.size(); s++){
		this->species_constraints[s].CCPotential->RhoValid = false;
	}
}

/**
 * @brief Set particle radii for a configuration.
 * 
 * @param radii Radii list matching particles.
 */
void MultiStealthyPotential::SetSecureRadii(const std::vector<double> &radii)
{
	assert(radii.size() == this->pConfig->NumParticle());
	this->radii = std::vector<double> (radii);
}

/**
 * @brief Initialize configuration and radii from a packing object.
 * 
 * @param pck Input sphere packing.
 */
void MultiStealthyPotential::SetPacking(const SpherePacking &pck)
{
	if (this->rc < 0.0){
		// an empty configuration
		this->idx2species.resize(pck.NumParticle(),spec_index(-1,0));
		std::cout << "species of all particles are set to be -1 (unspecified)\n";
		std::cout << "the critical radius of the contact potential is set to be "<< this->rc <<std::endl;
	}
	this->radii.resize(pck.NumParticle(), -1.0);
	this->rc = 0.0;
	for(size_t i =0; i < pck.NumParticle(); i++){
		this->rc = std::max(this-> rc , pck.GetCharacteristics(i));
		this->radii[i] = pck.GetCharacteristics(i);
	}
	this->rc *= 2.01;
	pConfig = std::make_unique<Configuration>(pck, "a");
	pConfig->PrepareIterateThroughNeighbors(this->rc);
	
	/* define sub-configurations for species. */
	this-> SubconfigValid = false;
	this-> UpdateSubconfigs();
	for (size_t s = 0; s < this->species_constraints.size(); s++){
		this->species_constraints[s].CCPotential->RhoValid = false;
	}
}

/**
 * @brief Rebuild sub-configurations for each species.
 */
void MultiStealthyPotential::UpdateSubconfigs()
{
	if (! this->SubconfigValid){
		Configuration temp (* this->pConfig);
		temp.RemoveParticles();
		this->subc.clear();
		for (size_t s = 0; s < this->species_constraints.size(); s ++ )
			this->subc.emplace_back(temp);
		
		for (size_t i = 0; i < this->idx2species.size(); i++){
			if (idx2species[i].species > -1){
				this->subc[this->idx2species[i].species].Insert("a", this->pConfig->GetRelativeCoordinates(i));
			}
		}

		for (size_t s = 0; s < this->species_constraints.size(); s ++ ){
			this->species_constraints[s].CCPotential->SetConfiguration(this->subc[s]);
		}
	}
	this->SubconfigValid = true;
}

/**
 * @brief Compute total energy from stealthy and repulsive contributions.
 * 
 * @return Total potential energy.
 */
double MultiStealthyPotential::Energy()
{
	double CC = 0.0, CP = 0.0;
	/* compute stealthy potential contribution (CC)*/
	for (size_t s = 0; s < this->species_constraints.size(); s ++ ){
		double CC_s = this->species_constraints[s].CCPotential->Energy();
		CC += ((double)this->species_constraints[s].idx.size()/(double)this->pConfig->NumParticle())*CC_s;
	}
	// 	double CC_s = 0.0;
	// 	EnergyConstraints * ptr_s = & (this->species_constraints[s]);
		
	// 	if(ptr_s->RhoValid == false){
	// 		this->GetRho(s);
	// 	}
		
	// 	// iterate over wavevectors 
	// 	#pragma omp parallel for schedule (guided, 1) num_threads(this->ParallelNumThread) reduction(+:CC_s)
	// 	for (size_t i = 0; i < ptr_s -> constraints.size(); i++){
	// 		CC_s += ptr_s -> constraints[i].V * (ptr_s -> constraints[i].rho_real*ptr_s -> constraints[i].rho_real + ptr_s -> constraints[i].rho_imag * ptr_s -> constraints[i].rho_imag );
	// 	}

	// 	CC += CC_s ;
	// }	
	// CC /= this->Volume;

	/* compute the contact potential contribution (CP)*/
	if (this->WithRepulsion){
#pragma omp parallel for schedule (guided, 1) num_threads(this->ParallelNumThread) reduction(+:CP)
		for (size_t i = 0; i < this->pConfig->NumParticle(); i++){

			double E = 0.0L;
			pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc, 
				[&i, &E, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
				/*Prevent double-counting */
				if (i < SourceParticle) { 
					double r = sqrt(shift.Modulus2()) / (radii[i] + radii[SourceParticle]);
					if (r < 1.0){
						E += (1.0 - r)*(1.0 - r);
					}
				}
				}
			);
			
			CP += E;
		}
	} 
	return this->u0 * CC + this->e0 * CP;
}

/**
 * @brief Evaluate force on a single particle.
 * 
 * @param force Output force vector.
 * @param i Particle index.
 */
void MultiStealthyPotential::Force(GeometryVector &force, size_t i)
{
	// stealthy potential contribution
	GeometryVector force_CC(static_cast<int> (this->Dimension));
	size_t s = (size_t) this->idx2species[i].species; // species
	size_t idx_in_s = this->idx2species[i].index;
	this->species_constraints[s].CCPotential -> Force(force_CC, idx_in_s);

	force_CC.MultiplyFrom(this->u0 * (double)this->species_constraints[s].idx.size()/(double)this->pConfig->NumParticle());

	// EnergyConstraints * ptr_s = & (this->species_constraints[s]);
	// if(ptr_s->RhoValid == false){
	// 	this->GetRho(s);
	// }
	
	// // iterate over wavevectors 
	// for (int q = 0; q < ptr_s -> constraints.size(); q++){
	// 	GeometryVector k = ptr_s->constraints[q].k;
	// 	GeometryVector r = this->pConfig->GetCartesianCoordinates(i);

	// 	force_CC = force_CC + k* ptr_s->constraints[q].V * (ptr_s->constraints[q].rho_real * sin(k.Dot(r))  +  ptr_s->constraints[q].rho_imag * cos(k.Dot(r)));
	// }
	// force_CC.MultiplyFrom(this->u0/this->Volume);

	// contact potential contribution (CP)
	GeometryVector force_CP(static_cast<int> (this->Dimension));
	if (this->WithRepulsion){
		pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->rc,
			[&i, &force_CP, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
			/*Prevent self-interaction */
			if (i != SourceParticle) {
				double r = sqrt(shift.Modulus2());
				double sigma = this->radii[i] + this->radii[SourceParticle];

				if (r < 1e-15){
					/* two particles are completely overlapped.*/
					GeometryVector x(shift);
					x.x[0] = (i > SourceParticle)? -1.0 : 1.0 ;
					force_CP = force_CP - 2.0 / sigma * x;
				}
				else if ( r < sigma ){
					force_CP = force_CP + 2.0 / sigma * (1.0 / r - 1.0 / sigma) * (-1.*shift);
				}
			}
		});
	}
	force_CP.MultiplyFrom(this->e0);

	force = force_CP + force_CC;
}

/**
 * @brief Evaluate forces for all particles.
 * 
 * @param results Output vector of forces.
 */
void MultiStealthyPotential::AllForce(std::vector<GeometryVector> &results)
{
	long end = this->pConfig->NumParticle();
	results.resize(end);
	unsigned long tempNumThread = 1;
	std::swap(tempNumThread, this->ParallelNumThread);
#pragma omp parallel for num_threads(tempNumThread) default(shared)
	for (long i = 0; i<end; i++)
		this->Force(results[i], i);
	std::swap(tempNumThread, this->ParallelNumThread);
}
