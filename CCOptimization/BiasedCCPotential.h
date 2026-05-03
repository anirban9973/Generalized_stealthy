/**
 * BiasedCCPotential: wraps an existing CC potential and adds a radial exclusion
 * field centred at the box centre.  Used by feasibility_scan (hole CLI mode) to
 * find the maximum void radius R_c compatible with stealthiness.
 *
 * Exclusion field (phi_ex):
 *   phi_ex = sum_{i: r_ic < R_f} ( R_f / r_ic - 1 )
 *   r_ic   = |r_i - r_c|,   r_c = 0.5 * (a1 + a2 + ... + ad)
 *
 * Force convention: AllForce() returns  -grad(Phi_total),  matching the rest of
 * the codebase (NLOPT uses  grad = -Force  inside RelaxStructure_NLOPT_inner).
 */

#ifndef BIASEDCCPOTENTIAL_H
#define BIASEDCCPOTENTIAL_H

#include <cmath>
#include <cassert>
#include <vector>
#include <Potential.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>

class BiasedCCPotential : public Potential
{
public:
    Potential *            pWrapped;
    double                 R_f;
    const Configuration *  pConfig;

    BiasedCCPotential(Potential * wrapped, double R_f_init, DimensionType dim)
        : Potential(dim), pWrapped(wrapped), R_f(R_f_init), pConfig(nullptr)
    {}

    void set_R_f(double r) { R_f = r; }

    // Box centre from lattice vectors (not particle centroid).
    GeometryVector BoxCenter() const
    {
        assert(pConfig != nullptr);
        GeometryVector rc(this->Dimension);
        for (DimensionType i = 0; i < this->Dimension; i++)
            rc.AddFrom(0.5 * pConfig->GetBasisVector(i));
        return rc;
    }

    // Delegate to wrapped potential; wrapped sets RhoValid=false.
    // Never call this inside Energy() or AllForce() — that would thrash the cache.
    virtual void SetConfiguration(const Configuration & c) override
    {
        this->pConfig = &c;
        this->pWrapped->SetConfiguration(c);
    }

    // Energy = E_cc + phi_ex.
    // Calling wrapped->Energy() triggers GetRho() (sets RhoValid=true) if needed.
    virtual double Energy() override
    {
        assert(pConfig != nullptr);
        double E_cc = this->pWrapped->Energy();

        GeometryVector rc = BoxCenter();
        double phi_ex = 0.0;
        size_t N = pConfig->NumParticle();
        for (size_t i = 0; i < N; i++)
        {
            GeometryVector ric = pConfig->GetCartesianCoordinates(i) - rc;
            double dist = std::sqrt(ric.Modulus2());
            if (dist < R_f)
                phi_ex += R_f / dist - 1.0;
        }
        return E_cc + phi_ex;
    }

    // Force on particle i: F_cc + F_ex.
    // F_ex = +R_f/r_ic^3 * (r_i - r_c)  for r_ic < R_f  (repulsive, outward).
    virtual void Force(GeometryVector & result, size_t i) override
    {
        assert(pConfig != nullptr);
        this->pWrapped->Force(result, i);

        GeometryVector ric = pConfig->GetCartesianCoordinates(i) - BoxCenter();
        double dist2 = ric.Modulus2();
        double dist  = std::sqrt(dist2);
        if (dist < R_f)
            result.AddFrom((R_f / (dist2 * dist)) * ric);
    }

    // All-particle forces.  Calls wrapped->AllForce() first (uses cached rho),
    // then adds F_ex for particles inside R_f.
    virtual void AllForce(std::vector<GeometryVector> & results) override
    {
        assert(pConfig != nullptr);
        this->pWrapped->AllForce(results);   // cache hit — rho already computed by Energy()

        GeometryVector rc = BoxCenter();
        size_t N = pConfig->NumParticle();
        for (size_t i = 0; i < N; i++)
        {
            GeometryVector ric = pConfig->GetCartesianCoordinates(i) - rc;
            double dist2 = ric.Modulus2();
            double dist  = std::sqrt(dist2);
            if (dist < R_f)
                results[i].AddFrom((R_f / (dist2 * dist)) * ric);
        }
    }
};

/**
 * Push every particle that sits strictly inside R_f radially outward to
 * R_f * (1 + 1e-6), so that the first NLOPT energy/gradient evaluation is
 * finite.  Particles exactly at the centre get a random outward nudge.
 *
 * Coordinate conversion:  relative[j] = r_cart . GetReciprocalBasisVector(j)
 */
inline void PushParticlesOutside(Configuration & config, double R_f,
                                  RandomGenerator & rng)
{
    DimensionType d = config.GetDimension();
    GeometryVector rc(d);
    for (DimensionType i = 0; i < d; i++)
        rc.AddFrom(0.5 * config.GetBasisVector(i));

    const double target = R_f * (1.0 + 1e-6);
    size_t N = config.NumParticle();

    for (size_t i = 0; i < N; i++)
    {
        GeometryVector ri  = config.GetCartesianCoordinates(i);
        GeometryVector ric = ri - rc;
        double dist = std::sqrt(ric.Modulus2());

        if (dist >= R_f) continue;   // already outside, nothing to do

        GeometryVector new_cart(d);
        if (dist < 1e-12)
        {
            // Particle exactly at centre: pick a random direction in each component
            for (DimensionType k = 0; k < d; k++)
                ric.x[k] = rng.RandomDouble() - 0.5;
            double len = std::sqrt(ric.Modulus2());
            if (len < 1e-15) ric.x[0] = 1.0;   // degenerate fallback
            len = std::sqrt(ric.Modulus2());
            new_cart = rc + (target / len) * ric;
        }
        else
        {
            new_cart = rc + (target / dist) * ric;
        }

        // Convert new Cartesian position to relative coordinates:
        //   rel[j] = new_cart . b_j   where b_j = GetReciprocalBasisVector(j)
        GeometryVector rel_new(d);
        for (DimensionType j = 0; j < d; j++)
            rel_new.x[j] = new_cart.Dot(config.GetReciprocalBasisVector(j));

        config.MoveParticle(i, rel_new);
    }
}

#endif // BIASEDCCPOTENTIAL_H
