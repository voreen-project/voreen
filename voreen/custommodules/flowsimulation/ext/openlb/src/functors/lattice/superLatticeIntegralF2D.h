/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_H
#define SUPER_LATTICE_INTEGRAL_F_2D_H

#include <vector>

#include "functors/genericF.h"
#include "blockLatticeIntegralF2D.h"
#include "superBaseF2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "indicator/superIndicatorF2D.h"
#include "functors/analytical/interpolationF2D.h"
#include "functors/lattice/reductionF2D.h"
#include "integral/superIntegralF2D.h"
#include "core/superLattice2D.h"
#include "core/vector.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry.h"
#include "superGeometryFaces2D.h"
#include "utilities/functorPtr.h"
#include "latticePhysBoundaryForce2D.h"
#include "latticePhysCorrBoundaryForce2D.h"

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// functor to get pointwise phys force acting on a indicated boundary on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDrag2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF2D<T>>              _indicatorF;
  SuperGeometryFaces2D<T>                       _facesF;
  SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR> _pBoundForceF;
  SuperSum2D<T,T>                               _sumF;

  const T _factor;
public:
  SuperLatticePhysDrag2D(SuperLattice<T,DESCRIPTOR>&      sLattice,
                         FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysDrag2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                         SuperGeometry<T,2>& superGeometry, const int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);

  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a indicated boundary on local lattice
/**
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCorrDrag2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF2D<T>>                  _indicatorF;
  SuperGeometryFaces2D<T>                           _facesF;
  SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> _pBoundForceF;
  SuperSum2D<T,T>                                   _sumF;

  const T _factor;
public:
  SuperLatticePhysCorrDrag2D(SuperLattice<T,DESCRIPTOR>&      sLattice,
                             FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysCorrDrag2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                             SuperGeometry<T,2>& superGeometry, const int material,
                             const UnitConverter<T,DESCRIPTOR>& converter);

  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
