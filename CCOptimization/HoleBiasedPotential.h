/**
 * Author : Anirban
 * Date   : May 2026
 *
 * Combined potential:  Phi = Phi_s + Phi_ex(Rf)
 *
 * Phi_s  = MultiShellS0Potential  (k-space CC constraints + soft-core repulsion)
 * Phi_ex = ExclusionFieldPotential (real-space exclusion field at box center)
 *
 * EnergyS() and EnergyEx() expose the two components independently so that
 * the R_f sweep can check both Phi_s < tol AND Phi_ex < tol.
 */

#ifndef HoleBiasedPotential_H
#define HoleBiasedPotential_H

#include <Potential.h>
#include "MultiShellS0Potential.h"
#include "ExclusionFieldPotential.h"

class HoleBiasedPotential : public Potential
{
public:
    MultiShellS0Potential  * phi_s;
    ExclusionFieldPotential * phi_ex;

    HoleBiasedPotential(MultiShellS0Potential * ps, ExclusionFieldPotential * pex)
        : Potential(ps->Dimension), phi_s(ps), phi_ex(pex) {}

    virtual ~HoleBiasedPotential() {}

    double EnergyS()  { return phi_s->Energy();  }
    double EnergyEx() { return phi_ex->Energy(); }

    virtual double Energy() override
    {
        return phi_s->Energy() + phi_ex->Energy();
    }

    virtual void Force(GeometryVector & force, size_t i) override
    {
        GeometryVector f_s, f_ex;
        phi_s->Force(f_s, i);
        phi_ex->Force(f_ex, i);
        force = f_s + f_ex;
    }

    virtual void AllForce(std::vector<GeometryVector> & result) override
    {
        phi_s->AllForce(result);          // initialises result
        phi_ex->AllForce(result);         // accumulates on top
    }

    virtual void SetConfiguration(const Configuration & c) override
    {
        phi_s->SetConfiguration(c);
        phi_ex->SetConfiguration(c);
    }
};

#endif
