/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Benjamin Förster, Adrian Kummerländer
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

#ifndef SUPER_INDICATOR_BASE_F_3D_H
#define SUPER_INDICATOR_BASE_F_3D_H

#include "functors/genericF.h"
#include "functors/lattice/superBaseF3D.h"
#include "communication/superStructure3D.h"

namespace olb {


template<typename T, typename W> class SuperF3D;
template<typename T> class SuperGeometry3D;
template<typename T> class BlockIndicatorF3D;


/// Base indicator functor (discrete)
/**
 * Provides Union, Without and Intersection arithmetic.
 * _Note: `operator()` must be overloaded by child classes._
 *
 * Implementations should maintain BlockIndicatorF3D instances as required.
 */
template <typename T>
class SuperIndicatorF3D : public SuperF3D<T,bool> {
protected:
  SuperGeometry3D<T>& _superGeometry;
public:
  using SuperF3D<T,bool>::operator();

  SuperIndicatorF3D(SuperGeometry3D<T>& geometry);
  /**
   * Get block indicator
   *
   * \returns _blockF[iCloc] cast as BlockIndicatorF3D<T>&
   **/
  BlockIndicatorF3D<T>& getBlockIndicatorF(int iCloc);
  /**
   * Get underlying super geometry
   *
   * \returns _superGeometry
   **/
  SuperGeometry3D<T>& getSuperGeometry();
  /**
   * Indicator specific function operator overload
   *
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described domain.
   **/
  virtual bool operator() (const int input[]);
  /**
   * Indicator specific function operator overload
   *
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described domain.
   **/
  virtual bool operator() (int iC, int iX, int iY, int iZ);
};


}

#endif
