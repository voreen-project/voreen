/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause, Davide Dapelo
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

#ifndef EUL2LAGR_OPERATION_H
#define EUL2LAGR_OPERATION_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>

#include "particles/subgrid3DLegacyFramework/particleOperations/particleOperations3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
class Eul2LagrOperation3D final : public ParticleOperation3D<T,PARTICLETYPE> {
public:
  Eul2LagrOperation3D(SuperLattice<T,DESCRIPTOR>& sLattice, int nParticles=1);
  virtual void applyParticleOperation(typename std::deque<PARTICLETYPE<T> >::iterator p, ParticleSystem3D<T, PARTICLETYPE>& pSys) override;
private:
  int _nParticles;
  SuperLattice<T,DESCRIPTOR>& _sLattice;
};

}


#endif

