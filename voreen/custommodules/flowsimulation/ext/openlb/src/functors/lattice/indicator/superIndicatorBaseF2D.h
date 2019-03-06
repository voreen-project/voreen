/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerländer
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

#ifndef SUPER_INDICATOR_BASE_F_2D_H
#define SUPER_INDICATOR_BASE_F_2D_H

#include "functors/genericF.h"
#include "functors/lattice/superBaseF2D.h"
#include "communication/superStructure2D.h"

namespace olb {

template<typename T> class BlockIndicatorF2D;

template <typename T>
class SuperIndicatorF2D : public SuperF2D<T,bool> {
public:
  using SuperF2D<T,bool>::operator();

  SuperIndicatorF2D(SuperStructure2D<T>& superStructure);
  /**
   * Get block indicator
   *
   * \returns _blockF[iCloc] cast as BlockIndicatorF2D<T>&
   **/
  BlockIndicatorF2D<T>& getBlockIndicatorF(int iCloc);
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
  virtual bool operator() (int iC, int iX, int iY);
};

} // namespace olb

#endif
