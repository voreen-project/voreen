/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *                          Albert Mink, Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_PARTICLE_BASE_F_3D_HH
#define SUPER_PARTICLE_BASE_F_3D_HH

//#include "superBaseF3D.h"
//#include "blockBaseF3D.h"

namespace olb {

template <typename T, typename PARTICLETYPE>
SuperParticleSystemF3D<T,PARTICLETYPE>::SuperParticleSystemF3D(SuperParticleSystem3D<T,PARTICLETYPE>& superParticleSystem, int targetDim)
  : SuperF3D<T,T>(superParticleSystem, targetDim), _sParticleSystem(superParticleSystem)
{}

template <typename T, typename PARTICLETYPE>
SuperParticleSystem3D<T,PARTICLETYPE>& SuperParticleSystemF3D<T,PARTICLETYPE>::getSuperParticleSystem()
{
  return _sParticleSystem;
}

template<typename T, typename PARTICLETYPE>
bool SuperParticleSystemF3D<T, PARTICLETYPE>::operator()(
  T output[], const int input[])
{
  auto& load = this->_sParticleSystem.getLoadBalancer();

  if (load.isLocal(input[0])) {
    const int loc = load.loc(input[0]);

    return this->getBlockF(loc)(output, &input[1]);
  }
  else {
    return false;
  }
}


template <typename T, typename PARTICLETYPE>
SuperParticleSystemIdentity3D<T,PARTICLETYPE>::SuperParticleSystemIdentity3D(
  FunctorPtr<SuperParticleSystemF3D<T,PARTICLETYPE>>&& f)
  : SuperParticleSystemF3D<T,PARTICLETYPE>(f->getSuperParticleSystem(), f->getTargetDim()),
    _f(std::move(f))
{
  this->getName() = "Id(" + _f->getName() + ")";

  for (int iC = 0; iC < _f->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockParticleSystemIdentity3D<T,PARTICLETYPE>(
        static_cast<BlockParticleSystemF3D<T,PARTICLETYPE>&>(_f->getBlockF(iC))));
  }
}

template <typename T, typename PARTICLETYPE>
bool SuperParticleSystemIdentity3D<T,PARTICLETYPE>::operator()(T output[], const int input[])
{
  return _f(output, input);
}





} // end namespace olb

#endif
