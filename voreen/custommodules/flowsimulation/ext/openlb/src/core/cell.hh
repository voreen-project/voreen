/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * Definition of a LB cell -- generic implementation.
 */
#ifndef CELL_HH
#define CELL_HH

#include <algorithm>
#include "cell.h"
#include "util.h"

namespace olb {

////////////////////////// Class Cell /////////////////////////////

/** The possibility to default construct Cell objects facilitates
 * their use in various types of containers. However, they can not
 * be used directly after construction; the method defineDynamics()
 * must be used first.
 */
template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>::Cell()
  : dynamics(nullptr)
{
  iniPop();
  iniExternal();
}

/** This constructor initializes the dynamics, but not the values
 * of the distribution functions. Remember that the dynamics is not
 * owned by the Cell object, the user must ensure its proper
 * destruction and a sufficient life time.
 */
template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>::Cell(Dynamics<T,Lattice>* dynamics_)
  : dynamics(dynamics_)
{
  iniPop();
  iniExternal();
}

template<typename T, template<typename U> class Lattice>
void Cell<T,Lattice>::defineDynamics(Dynamics<T,Lattice>* dynamics_)
{
  dynamics = dynamics_;
}

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice> const* Cell<T,Lattice>::getDynamics() const
{
  OLB_PRECONDITION(dynamics);
  return dynamics;
}

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>* Cell<T,Lattice>::getDynamics()
{
  OLB_PRECONDITION(dynamics);
  return dynamics;
}

template<typename T, template<typename U> class Lattice>
void Cell<T,Lattice>::revert()
{
  for (int iPop=1; iPop<=Lattice<T>::q/2; ++iPop) {
    std::swap(this->f[iPop],this->f[iPop+Lattice<T>::q/2]);
  }
}

template<typename T, template<typename U> class Lattice>
void Cell<T,Lattice>::iniPop()
{
  for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
    this->f[iPop] = T();
  }
}

template<typename T, template<typename U> class Lattice>
void Cell<T,Lattice>::iniExternal()
{
  for (int iData=0; iData<Lattice<T>::ExternalField::numScalars; ++iData) {
    *external.get(iData) = T();
  }
}

template<typename T, template<typename U> class Lattice>
void Cell<T,Lattice>::serialize(T* data) const
{
  const int q = Lattice<T>::q;
  const int numExt = Lattice<T>::ExternalField::numScalars;
  for (int iPop=0; iPop<q; ++iPop) {
    data[iPop] = this->f[iPop];
  }
  for (int iExternal=0; iExternal < numExt; ++iExternal) {
    data[iExternal+q] = *external.get(iExternal);
  }
}

template<typename T, template<typename U> class Lattice>
void Cell<T,Lattice>::unSerialize(T const* data)
{
  const int q = Lattice<T>::q;
  const int numExt = Lattice<T>::ExternalField::numScalars;
  for (int iPop=0; iPop<q; ++iPop) {
    this->f[iPop] = data[iPop];
  }
  for (int iExternal=0; iExternal < numExt; ++iExternal) {
    *external.get(iExternal) = data[iExternal+q];
  }
}


template<typename T, template<typename U> class Lattice>
std::size_t Cell<T,Lattice>::getSerializableSize() const
{
  return sizeof(T) * Lattice<T>::q // this->f
         + sizeof(T) * Lattice<T>::ExternalField::numScalars; // externals
}


template<typename T, template<typename U> class Lattice>
bool* Cell<T,Lattice>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  this->registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->f[0], Lattice<T>::q);
  if (Lattice<T>::ExternalField::numScalars)
    this->registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *(external.get(0)),
        Lattice<T>::ExternalField::numScalars);

  return dataPtr;
}

}  // namespace olb


#endif
