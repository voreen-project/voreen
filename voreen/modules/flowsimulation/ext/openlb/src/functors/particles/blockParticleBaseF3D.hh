/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_PARTICLE_BASE_F_3D_HH
#define BLOCK_PARTICLE_BASE_F_3D_HH

//#include "blockBaseF3D.h"

namespace olb {


template <typename T, typename PARTICLETYPE>
BlockParticleSystemF3D<T,PARTICLETYPE>::BlockParticleSystemF3D
(ParticleSystem3D<T,PARTICLETYPE>& particleSystem, int targetDim)
  : BlockF3D<T>(particleSystem, targetDim), _particleSystem(particleSystem)
{ }

template <typename T, typename PARTICLETYPE>
ParticleSystem3D<T,PARTICLETYPE>& BlockParticleSystemF3D<T, PARTICLETYPE>::getParticleSystem()
{
  return _particleSystem;
}


template <typename T, typename PARTICLETYPE>
BlockParticleSystemIdentity3D<T,PARTICLETYPE>::BlockParticleSystemIdentity3D(
  BlockParticleSystemF3D<T,PARTICLETYPE>& f)
  : BlockParticleSystemF3D<T,PARTICLETYPE>(f.getBlockParticleSystem(),f.getTargetDim()),
    _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename PARTICLETYPE>
bool BlockParticleSystemIdentity3D<T,PARTICLETYPE>::operator()(T output[], const int input[])
{
  return _f(output,input);
}


} // end namespace olb

#endif
