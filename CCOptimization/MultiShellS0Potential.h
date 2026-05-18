/**
*	Author	: Jaeuk Kim / Anirban
*	Date	: May 2026 */

/** \file MultiShellS0Potential.h
 * \brief Collective-coordinate potential for M annular shells each with its own S0 target.
 *        Stealthy shells (S0=0) use ShiftedCCPotential; equiluminous shells (S0>0) use ShiftedCCPotential_S2.
 *        Optional soft-core repulsion is included (same form as RepulsiveCCPotential). */

#ifndef MultiShellS0Potential_H
#define MultiShellS0Potential_H

/* in $(core) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
/* in $(pot) */
#include <Potential.h>
#include <CollectiveCoordinatePotential.h>

class MultiShellS0Potential : public Potential
{
protected:
	const Configuration * pConfig;
	bool WithRepulsion = true;
public:
	double e0,		//< energy scale for soft-core repulsion
		sigma;		//< length scale for soft-core repulsion
	ShiftedCCPotential * CCPotential1;		//< stealthy shells (S0 = 0)
	ShiftedCCPotential * CCPotential2;		//< equiluminous shells (S0 > 0)

	MultiShellS0Potential(DimensionType dim, double epsilon, double sigma) : Potential(dim) {
		CCPotential1 = new ShiftedCCPotential(dim);
		CCPotential2 = new ShiftedCCPotential_S2(dim);
		this->e0 = epsilon;
		this->sigma = sigma;
		this->WithRepulsion = (std::abs(this->e0) > 1e-10);
		if (this->WithRepulsion)
			std::cout << "Soft-core repulsion is on\n";
		else
			std::cout << "Soft-core repulsion is off\n";
	}

	virtual ~MultiShellS0Potential() {}

	virtual double Energy() {
		double LJ = 0.0, CC = 0.0;
		if (this->WithRepulsion) {
#pragma omp parallel for schedule(guided, 1) num_threads(this->ParallelNumThread) reduction(+:LJ)
			for (int i = 0; i < (int)pConfig->NumParticle(); i++) {
				double E = 0.0;
				pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->sigma,
					[&i, &E, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
					if (i < (int)SourceParticle) {
						double r = sqrt(shift.Modulus2()) / this->sigma;
						if (r < 1.0)
							E += this->e0*(1.0 - r)*(1.0 - r);
					}
					});
				LJ += E;
			}
		}
		CC = this->CCPotential1->Energy() + this->CCPotential2->Energy();
		return LJ + CC;
	}

	virtual void Force(GeometryVector & force, size_t i) {
		GeometryVector force_LJ(static_cast<int>(this->Dimension));
		if (this->WithRepulsion) {
			pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->sigma,
				[&i, &force_LJ, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
				if (i != SourceParticle) {
					double r = sqrt(shift.Modulus2());
					if (r < 1e-15) {
						GeometryVector x(shift);
						x.x[0] = (i > SourceParticle) ? -1.0 : 1.0;
						force_LJ = force_LJ - 2.0*this->e0 / this->sigma * x;
					} else {
						force_LJ = force_LJ - 2.0*this->e0 / this->sigma * (1.0 / r - 1.0 / this->sigma) * shift;
					}
				}
				});
		}
		GeometryVector force_cc1, force_cc2;
		this->CCPotential1->Force(force_cc1, i);
		this->CCPotential2->Force(force_cc2, i);
		force = force_LJ + force_cc1 + force_cc2;
	}

	virtual void AllForce(std::vector<GeometryVector> & result) {
		this->CCPotential1->AllForce(result);
		std::vector<GeometryVector> result2;
		this->CCPotential2->AllForce(result2);
		for (size_t i = 0; i < result.size(); i++)
			result[i] = result[i] + result2[i];
		if (this->WithRepulsion) {
#pragma omp parallel for schedule(guided) num_threads(this->ParallelNumThread)
			for (int i = 0; i < (int)pConfig->NumParticle(); i++) {
				GeometryVector force_lj(static_cast<int>(this->Dimension));
				pConfig->IterateThroughNeighbors(pConfig->GetRelativeCoordinates(i), this->sigma,
					[&i, &force_lj, this](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, size_t SourceParticle)->void {
					if (i != (int)SourceParticle) {
						double r = sqrt(shift.Modulus2());
						if (r < 1e-15) {
							GeometryVector x(shift);
							x.x[0] = (i > (int)SourceParticle) ? -1.0 : 1.0;
							force_lj = force_lj - 2.0*this->e0 / this->sigma * x;
						} else {
							force_lj = force_lj - 2.0*this->e0 / this->sigma * (1.0 / r - 1.0 / this->sigma) * shift;
						}
					}
					});
				result[i] = result[i] + force_lj;
			}
		}
	}

	virtual void SetConfiguration(const Configuration & c) {
		this->pConfig = &c;
		if (this->WithRepulsion)
			this->pConfig->PrepareIterateThroughNeighbors(this->sigma);
		this->CCPotential1->SetConfiguration(c);
		this->CCPotential2->SetConfiguration(c);
	}
};

#endif
