/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_VOLUME_FRACTION_APPROXIMATION_2D_H
#define LATTICE_VOLUME_FRACTION_APPROXIMATION_2D_H

#include <vector>

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "indicator/blockIndicatorF2D.h"
#include "dynamics/porousBGKdynamics.h"

namespace olb {

/// functor to get pointwise an approx. for the volume fraction
template <typename T, typename DESCRIPTOR>
class SuperLatticeVolumeFractionApproximation2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeVolumeFractionApproximation2D(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
      IndicatorF2D<T>& indicator, int refinementLevel, const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut = false);
};

/**
 *  functor returns pointwise an approximation for the volume fraction
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticeVolumeFractionApproximation2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  IndicatorF2D<T>& _indicator;
  int _refinementLevel;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  bool _insideOut;
  T _physSubGridMinPhysRshift;
  T _physSubGridDeltaX;
  T _latticeSubGridVolume;
public:
  BlockLatticeVolumeFractionApproximation2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
      BlockGeometry<T,2>& blockGeometry,
      IndicatorF2D<T>& indicator,
      int refinementLevel,
      const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut);
  bool operator() (T output[], const int input[]);
};

}
#endif
