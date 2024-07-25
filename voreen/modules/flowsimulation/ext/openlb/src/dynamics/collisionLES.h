/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
 *                2021 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef DYNAMICS_COLLISION_LES_H
#define DYNAMICS_COLLISION_LES_H

namespace olb {

namespace collision {

namespace LES {

struct Smagorinsky : public descriptors::FIELD_BASE<1> { };

}

/// Implementations of meta-collisions modifying a COLLISION template argument
/**
 * This split is required due to language restrictions on template specialization resolution
 **/
namespace detail {

template <typename COLLISION, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct SmagorinskyEffectiveOmega {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  V computeEffectiveOmega(CELL& cell, PARAMETERS& parameters) any_platform {
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    const V rho = MomentaF().computeRho(cell);
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V smagorinsky = parameters.template get<collision::LES::Smagorinsky>();
    V piNeqNorm = util::sqrt(piNeqNormSqr);
    V preFactor = smagorinsky*smagorinsky
                * descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()
                * 2 * util::sqrt(2);
    /// Molecular realaxation time
    V tauMol = V{1} / omega;
    /// Turbulent realaxation time
    V tauTurb = V{0.5} * (util::sqrt(tauMol*tauMol + preFactor / rho * piNeqNorm) - tauMol);
    /// Effective realaxation time
    V tauEff = tauMol + tauTurb;
    return  V{1} / tauEff;
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    parameters.template set<descriptors::OMEGA>(
      computeEffectiveOmega(cell, parameters));
    return CollisionO().apply(cell, parameters);
  }
};

template <typename COLLISION, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct ShearSmagorinskyEffectiveOmega {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  V computeEffectiveOmega(CELL& cell, PARAMETERS& parameters) any_platform {
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    const V rho = MomentaF().computeRho(cell);
    const auto iT = parameters.template get<descriptors::LATTICE_TIME>();
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V smagorinsky = parameters.template get<collision::LES::Smagorinsky>();
    V avShear = cell.template getField<descriptors::AV_SHEAR>();
    V piNeqNorm = util::sqrt(piNeqNormSqr);
    V preFactor = smagorinsky*smagorinsky
                * descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()
                * 2 * util::sqrt(2);
    avShear = (avShear*iT + piNeqNorm) / (iT+1);
    V tauMol = V{1} / omega;
    V piNeqNormSISM = piNeqNorm - avShear;
    V tauTurb = V{0.5} * (util::sqrt(tauMol*tauMol+(preFactor*piNeqNormSISM/rho))-tauMol);
    V tauEff = tauMol + tauTurb;
    return  V{1} / tauEff;
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    parameters.template set<descriptors::OMEGA>(
      computeEffectiveOmega(cell, parameters));
    const auto iT = parameters.template get<descriptors::LATTICE_TIME>();
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    V avShear = cell.template getField<descriptors::AV_SHEAR>();
    avShear = (avShear*iT + util::sqrt(piNeqNormSqr)) / (iT+1);
    cell.template setField<descriptors::AV_SHEAR>(avShear);
    return CollisionO().apply(cell, parameters);
  }
};

template <typename COLLISION, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct ConStrainSmagorinskyEffectiveOmega {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  V computeEffectiveOmega(CELL& cell, PARAMETERS& parameters) any_platform {
    V pi[util::TensorVal<DESCRIPTOR>::n] { };
    MomentaF().computeStress(cell, pi);
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    const V rho = MomentaF().computeRho(cell);
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V smagorinsky = parameters.template get<collision::LES::Smagorinsky>();
    V piNeqNorm = util::sqrt(2*piNeqNormSqr);
    V S[util::TensorVal<DESCRIPTOR>::n];
    V cs2 = V{1} / descriptors::invCs2<V,DESCRIPTOR>();
    V tauMol = V{1} / omega;
    //Strain Tensor
    V Phi = (-0.5*(-rho*tauMol*cs2+util::sqrt(rho*rho*tauMol*tauMol*cs2*cs2+V{2}*(smagorinsky*smagorinsky)*rho*piNeqNorm))/(smagorinsky*smagorinsky*rho*piNeqNorm));
    for (int n=0; n < util::TensorVal<DESCRIPTOR>::n; ++n) {
      S[n] = Phi*pi[n];
    }
    //Strain Tensor Norm
    V SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
    }
    V SNorm = util::sqrt(2*SNormSqr);
    /// Turbulent realaxation time
    V tauTurb = smagorinsky*smagorinsky*SNorm/cs2;
    /// Effective realaxation time
    V tauEff = tauMol + tauTurb;
    return V{1} / tauEff;
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform{
    parameters.template set<descriptors::OMEGA>(
      computeEffectiveOmega(cell, parameters));
    return CollisionO().apply(cell, parameters);
  }

};

template <typename COLLISION, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct ConSmagorinskyEffectiveOmega {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  V computeEffectiveOmega(CELL& cell, PARAMETERS& parameters) any_platform {
    V pi[util::TensorVal<DESCRIPTOR>::n] { };
    MomentaF().computeStress(cell, pi);
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    const V rho = MomentaF().computeRho(cell);
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V smagorinsky = parameters.template get<collision::LES::Smagorinsky>();
    V H[util::TensorVal<DESCRIPTOR >::n];
    V conSmagoR[DESCRIPTOR::q];
    V S[util::TensorVal<DESCRIPTOR >::n];
    V tauMol = V{1} / omega;
    V cs2 = V{1} / descriptors::invCs2<V,DESCRIPTOR>();
    V piNeqNorm = util::sqrt(2*piNeqNormSqr);
    V Phi = (-0.5*(-rho*tauMol*cs2+util::sqrt(rho*rho*tauMol*tauMol*cs2*cs2+V{2}*(smagorinsky*smagorinsky)*rho*piNeqNorm))/(smagorinsky*smagorinsky*rho*piNeqNorm));
    for (int n=0; n < util::TensorVal<DESCRIPTOR>::n; ++n) {
      S[n] = Phi*pi[n];
    }
    //Strain rate Tensor Norm
    V SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
    }
    V SNorm = util::sqrt(2*SNormSqr);
    V preFactor = smagorinsky*smagorinsky
                * descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()
                * 2 * util::sqrt(2);
    //consistent Samagorinsky additional R term
    //for (int q=0; q < DESCRIPTOR::q; ++q) {
    {
      unsigned q = 0;
      V t = descriptors::t<V,DESCRIPTOR>(q);
      //Hermite-Polynom H = c*c-cs^2*kronDelta
      H[0] = descriptors::c<DESCRIPTOR>(q,0)*descriptors::c<DESCRIPTOR>(q,0)-cs2;
      H[1] = descriptors::c<DESCRIPTOR>(q,0)*descriptors::c<DESCRIPTOR>(q,1);
      H[2] = descriptors::c<DESCRIPTOR>(q,1)*descriptors::c<DESCRIPTOR>(q,1)-cs2;
      if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
        H[2] = descriptors::c<DESCRIPTOR>(q,0)*descriptors::c<DESCRIPTOR>(q,2);
        H[3] = descriptors::c<DESCRIPTOR>(q,1)*descriptors::c<DESCRIPTOR>(q,1)-cs2;
        H[4] = descriptors::c<DESCRIPTOR>(q,1)*descriptors::c<DESCRIPTOR>(q,2);
        H[5] = descriptors::c<DESCRIPTOR>(q,2)*descriptors::c<DESCRIPTOR>(q,2)-cs2;
      }
      //contraction or scalar product H*S
      V contractHS = H[0]*S[0] + 2.0*H[1]*S[1] + H[2]*S[2];
      if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
        contractHS += H[2]*S[2] + H[3]*S[3] + 2.0*H[4]*S[4] + H[5]*S[5];
      }
      //additional term
      conSmagoR[q] = t*preFactor*SNorm*contractHS;
    }
    return conSmagoR[0];
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    parameters.template set<descriptors::OMEGA>(
      computeEffectiveOmega(cell, parameters));
    return CollisionO().apply(cell, parameters);
  }
};

template <typename COLLISION, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct WaleEffectiveOmega {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  V computeEffectiveOmega(CELL& cell, PARAMETERS& parameters) any_platform {
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    const V rho = MomentaF().computeRho(cell);
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V smagorinsky = parameters.template get<collision::LES::Smagorinsky>();
    const auto velocityGradient = cell.template getField<descriptors::VELO_GRAD>();
    V piNeqNorm = util::sqrt(piNeqNormSqr);
    V preFactor = smagorinsky*smagorinsky;
    // velocity gradient tensor
    V g[3][3];
    for (unsigned i=0; i < 3; i++) {
      for (unsigned j=0; j < 3; j++) {
        g[i][j] = velocityGradient[i*3 + j];
      }
    }
    // strain rate tensor
    V s[3][3];
    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
        s[i][j] = (g[i][j] + g[j][i]) / V{2};
      }
    }
    // traceless symmetric part of the square of the velocity gradient tensor
    V G[3][3];
    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
        G[i][j] = V{0};
      }
    }
    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
        for (unsigned k = 0; k < 3; k++) {
          G[i][j] += (g[i][k]*g[k][j] + g[j][k]*g[k][i]) / V{2};
        }
      }
    }
    V trace{};
    for (unsigned i = 0; i < 3; i++) {
      trace += V{1}/V{3} * g[i][i] * g[i][i];
    }
    for (unsigned i = 0; i < 3; i++) {
      G[i][i] -= trace;
    }
    // inner product of the traceless symmetric part of the square of the velocity gradient tensor
    V G_ip{};
    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
        G_ip = G[i][j] * G[i][j];
      }
    }
    // inner product of the strain rate
    V s_ip{};
    for (unsigned i = 0; i < 3; i++) {
      for (unsigned j = 0; j < 3; j++) {
        s_ip = s[i][j] * s[i][j];
      }
    }
    // Turbulent relaxation time
    V tauTurb = V{3} * preFactor * (util::pow(G_ip,1.5) / (util::pow(s_ip,2.5) + util::pow(G_ip,1.25)));
    if ((util::pow(s_ip,2.5) + util::pow(G_ip,1.25)) == 0) {
      tauTurb = 0;
    }
    // Physical turbulent viscosity must be equal or higher that zero
    if (tauTurb < 0) {
      tauTurb = 0;
    }
    /// Molecular relaxation time
    V tauMol = V{1} / omega;
    /// Effective relaxation time
    V tauEff = tauMol + tauTurb;
    return V{1} / tauEff;
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    parameters.template set<descriptors::OMEGA>(
      computeEffectiveOmega(cell, parameters));
    return CollisionO().apply(cell, parameters);
  }
};

template <typename COLLISION, typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
struct KrauseEffectiveOmega {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    V rho, u[DESCRIPTOR::d], fNeq[DESCRIPTOR::q] { };
    MomentaF().computeRhoU(cell, rho, u);
    lbm<DESCRIPTOR>::computeFneq(cell, fNeq, rho, u); // TODO: Use EQUILIBRIUM instead of 2nd order
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V smagorinsky = parameters.template get<collision::LES::Smagorinsky>();
    const V preFactor = smagorinsky*smagorinsky
                      * 3*descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()
                      * 2*util::sqrt(2);
    FieldD<V,DESCRIPTOR,typename COLLISION::OMEGA> omegaEff;
    const V tauMol = V{1} / omega;
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V tauTurb = V{0.5}*(util::sqrt(tauMol*tauMol + preFactor/rho * util::fabs(fNeq[iPop])) - tauMol);
      omegaEff[iPop] = V{1} / (tauMol + tauTurb);
    }
    parameters.template set<typename COLLISION::OMEGA>(omegaEff);
    return CollisionO().apply(cell, parameters);
  }
};


}

/// Compute dynamics parameter OMEGA locally using Smagorinsky LES model
template <typename COLLISION>
struct SmagorinskyEffectiveOmega {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, LES::Smagorinsky
  >;

  static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                "COLLISION must be parametrized using relaxation frequency OMEGA");

  static std::string getName() {
    return "SmagorinskyEffectiveOmega<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using type = detail::SmagorinskyEffectiveOmega<COLLISION,DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
};

/// Compute dynamics parameter OMEGA locally using Shear Smagorinsky LES model
template <typename COLLISION>
struct ShearSmagorinskyEffectiveOmega {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::LATTICE_TIME, descriptors::OMEGA, LES::Smagorinsky
  >;

  static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                "COLLISION must be parametrized using relaxation frequency OMEGA");

  static std::string getName() {
    return "ShearSmagorinskyEffectiveOmega<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using type = detail::ShearSmagorinskyEffectiveOmega<COLLISION,DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
};

/// Compute dynamics parameter OMEGA locally using Consistent Strain Smagorinsky LES model
template <typename COLLISION>
struct ConStrainSmagorinskyEffectiveOmega {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, LES::Smagorinsky
  >;

  static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                "COLLISION must be parametrized using relaxation frequency OMEGA");

  static std::string getName() {
    return "ConStrainSmagorinskyEffectiveOmega<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using type = detail::ConStrainSmagorinskyEffectiveOmega<COLLISION,DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
};

/// Compute dynamics parameter OMEGA locally using Consistent Smagorinsky LES model
template <typename COLLISION>
struct ConSmagorinskyEffectiveOmega {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, LES::Smagorinsky
  >;

  static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                "COLLISION must be parametrized using relaxation frequency OMEGA");

  static std::string getName() {
    return "ConSmagorinskyEffectiveOmega<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using type = detail::ConSmagorinskyEffectiveOmega<COLLISION,DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
};

/// Compute dynamics parameter OMEGA locally using WALE
template <typename COLLISION>
struct WaleEffectiveOmega {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, LES::Smagorinsky
  >;

  static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                "COLLISION must be parametrized using relaxation frequency OMEGA");

  static std::string getName() {
    return "WaleEffectiveOmega<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using type = detail::WaleEffectiveOmega<COLLISION,DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
};

/// Compute dynamics parameter OMEGA locally using Krause LES model
template <typename COLLISION>
struct KrauseEffectiveOmega {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, LES::Smagorinsky
  >;

  static std::string getName() {
    return "KrauseEffectiveOmega<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using type = detail::KrauseEffectiveOmega<COLLISION,DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
};

}

}

#endif

#include "collisionLES.cse.h"
