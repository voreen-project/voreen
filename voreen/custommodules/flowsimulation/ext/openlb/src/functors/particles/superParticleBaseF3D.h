/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2021 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *                          Albert Mink, Benjamin Förster, Adrian Kummerlaender, Nicolas Hafen
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

#ifndef SUPER_PARTICLE_BASE_F_3D_H
#define SUPER_PARTICLE_BASE_F_3D_H

#include <memory>

//#include "functors/genericF.h"
#include "blockParticleBaseF3D.h"
//#include "indicator/superIndicatorBaseF3D.h"
//#include "communication/superStructure3D.h"
//#include "core/superData3D.h"
//
#include "particles/superParticleSystem3D.h"

/* Note: Throughout the whole source code directory functors, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


using namespace particles; //TODO: this should be solved more elegantly!

/// represents all functors that operate on a SuperLattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, typename PARTICLETYPE>
class SuperParticleSystemF3D : public SuperF3D<T,T> { //ORIG
//class SuperParticleSystemF3D : public GenericF<T,int> {
protected:
  SuperParticleSystemF3D( SuperParticleSystem3D<T,PARTICLETYPE>& sParticleSystem, int targetDim);

  SuperParticleSystem3D<T,PARTICLETYPE>& _sParticleSystem;
public:
  using identity_functor_type = SuperLatticeIdentity3D<T,PARTICLETYPE>;

  SuperParticleSystem3D<T,PARTICLETYPE>& getSuperParticleSystem();

  bool operator() (T output [], const int input []);

  using GenericF<T,int>::operator();
};

/// identity functor for memory management
template <typename T, typename PARTICLETYPE>
class SuperParticleSystemIdentity3D : public SuperParticleSystemF3D<T,PARTICLETYPE> {
protected:
  FunctorPtr<SuperParticleSystemF3D<T,PARTICLETYPE>> _f;
public:
  SuperParticleSystemIdentity3D(FunctorPtr<SuperParticleSystemF3D<T,PARTICLETYPE>>&& f);
  bool operator() (T output[], const int input[]) override;
};




} // end namespace olb

#endif
