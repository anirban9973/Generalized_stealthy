/**
 * Author : Anirban
 * Date   : May 2026
 *
 * Real-space exclusion field that pushes particles away from the box center,
 * biasing configurations toward a large hole at the center.
 *
 * Energy:   Phi_ex = sum_i F(Rf, r_ic)
 *   where F = (Rf/r_ic - 1)^2   if r_ic < Rf
 *           = 0                  otherwise
 * and r_ic = |r_i - r_center|, r_center = 0.5*(a1 + a2 + ...).
 *
 * Force on particle i (= -grad_i Phi_ex):
 *   F_i = +2*Rf/r_ic^2 * (Rf/r_ic - 1) * r_hat_ic    if r_ic < Rf
 *       = 0                                             otherwise
 * Points away from the center (repulsive).
 */

#ifndef ExclusionFieldPotential_H
#define ExclusionFieldPotential_H

#include <Potential.h>
#include <GeometryVector.h>
#include <cmath>
#include <vector>

class ExclusionFieldPotential : public Potential
{
protected:
    const Configuration * pConfig;
    GeometryVector center_;   // box center in Cartesian coordinates
    double Rf_;               // exclusion radius

    void computeCenter()
    {
        DimensionType d = pConfig->GetDimension();
        center_ = GeometryVector(static_cast<int>(d));
        for (DimensionType i = 0; i < d; i++)
            center_ = center_ + pConfig->GetBasisVector(i) * 0.5;
    }

public:
    ExclusionFieldPotential(DimensionType dim, double Rf)
        : Potential(dim), pConfig(nullptr), center_(static_cast<int>(dim)), Rf_(Rf) {}

    void set_Rf(double Rf) { Rf_ = Rf; }
    double get_Rf() const  { return Rf_; }

    virtual void SetConfiguration(const Configuration & c) override
    {
        pConfig = &c;
        computeCenter();
    }

    virtual double Energy() override
    {
        double E = 0.0;
        size_t N = pConfig->NumParticle();
        for (size_t i = 0; i < N; i++) {
            GeometryVector diff = pConfig->GetCartesianCoordinates(i) - center_;
            double ric2 = diff.Modulus2();
            if (ric2 < 1e-24) continue;
            double ric = std::sqrt(ric2);
            if (ric < Rf_) {
                double ratio = Rf_ / ric - 1.0;
                E += ratio * ratio;
            }
        }
        return E;
    }

    virtual void Force(GeometryVector & force, size_t i) override
    {
        force = GeometryVector(static_cast<int>(this->Dimension));
        GeometryVector diff = pConfig->GetCartesianCoordinates(i) - center_;
        double ric2 = diff.Modulus2();
        if (ric2 < 1e-24) return;
        double ric = std::sqrt(ric2);
        if (ric >= Rf_) return;
        // F_i = +2*Rf/ric^2*(Rf/ric - 1)*r_hat = +2*Rf/ric^3*(Rf/ric - 1)*diff
        double coeff = 2.0 * Rf_ / (ric2 * ric) * (Rf_ / ric - 1.0);
        force = diff * coeff;
    }

    virtual void AllForce(std::vector<GeometryVector> & result) override
    {
        size_t N = pConfig->NumParticle();
        if (result.size() != N)
            result.assign(N, GeometryVector(static_cast<int>(this->Dimension)));
        for (size_t i = 0; i < N; i++) {
            GeometryVector f;
            Force(f, i);
            result[i] = result[i] + f;
        }
    }
};

#endif
