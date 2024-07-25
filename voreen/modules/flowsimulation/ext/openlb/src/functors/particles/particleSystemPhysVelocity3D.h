/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
 *
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

#ifndef PARTICLE_SYSTEM_PHYS_VELOCITY_3D_H
#define PARTICLE_SYSTEM_PHYS_VELOCITY_3D_H

#include<vector>

#include "superParticleBaseF3D.h"
//#include "superCalcF3D.h"
//#include "functors/analytical/indicator/indicatorBaseF3D.h"
//
#include "blockParticleBaseF3D.h"
//#include "geometry/blockGeometry3D.h"
//#include "functors/analytical/indicator/indicatorBaseF3D.h"
//#include "indicator/blockIndicatorBaseF3D.h"
//#include "dynamics/smagorinskyBGKdynamics.h"
//#include "dynamics/porousBGKdynamics.h"
#include "particles/superParticleSystem3D.h"

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */


namespace olb {

using namespace particles; //TODO: this should be solved more elegantly!

/// functor to get pointwise phys velocity //TODO: adapt description
template <typename T, typename PARTICLETYPE>
class SuperParticleSystemPhysVelocity3D final : public SuperParticleSystemF3D<T,PARTICLETYPE> {
public:
  SuperParticleSystemPhysVelocity3D(SuperParticleSystem3D<T,PARTICLETYPE>& sParticleSystem,
                             bool print=false);
private:
  bool _print;
};





/// functor returns pointwise phys velocity //TODO: adapt description
template <typename T, typename PARTICLETYPE>
class BlockParticleSystemPhysVelocity3D final : public BlockParticleSystemF3D<T,PARTICLETYPE> {
private:
  const int  _overlap;
  const bool _print;
public:
  BlockParticleSystemPhysVelocity3D(ParticleSystem3D<T,PARTICLETYPE>& particleSystem,
                             int overlap,
                             bool print=false);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
