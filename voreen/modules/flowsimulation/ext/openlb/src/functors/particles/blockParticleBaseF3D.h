/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Lukas Baron, Mathias J. Krause
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

#ifndef BLOCK_PARTICLE_BASE_F_3D_H
#define BLOCK_PARTICLE_BASE_F_3D_H

//#include "functors/genericF.h"
//#include "core/blockData3D.h"
//#include "core/blockStructure3D.h"
//#include "core/blockLatticeStructure3D.h"
//#include "core/unitConverter.h"

#include "particles/superParticleSystem3D.h"
#include <memory>

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

using namespace particles; //TODO: this should be solved more elegantly!

/// functor to extract one component
template <typename T>
class BlockParticleSystem3D final : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
public:
  BlockParticleSystem3D(BlockF3D<T>& f);
  // access operator should not delete f, since f still has the identity as child
  bool operator() (T output[], const int input[]) override;
};

/// represents all functors that operate on a PARTICLETYPE
template <typename T, typename PARTICLETYPE>
class BlockParticleSystemF3D : public BlockF3D<T> {
protected:
  BlockParticleSystemF3D( ParticleSystem3D<T,PARTICLETYPE>& particleSystem, int targetDim);
  ParticleSystem3D<T,PARTICLETYPE>& _particleSystem;
public:
  ParticleSystem3D<T,PARTICLETYPE>& getParticleSystem();
};

/// identity functor
template <typename T, typename PARTICLETYPE>
class BlockParticleSystemIdentity3D final : public BlockParticleSystemF3D<T,PARTICLETYPE> {
protected:
  BlockParticleSystemF3D<T,PARTICLETYPE>& _f;
public:
  BlockParticleSystemIdentity3D(BlockParticleSystemF3D<T,PARTICLETYPE>& f);
  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
