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

#ifndef LATTICE_POROUS_MOMENTUM_LOSS_FORCE_2D_H
#define LATTICE_POROUS_MOMENTUM_LOSS_FORCE_2D_H

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

/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row (described are return values for the first particle).
 * \return output[0]-output[1] translational force - physical units
 * \return output[3] torque - physical units
 * \return output[4] number of voxels
 */
/*
template <typename T, typename DESCRIPTOR>
class SuperLatticePorousMomentumLossForce2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePorousMomentumLossForce2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                        SuperGeometry<T,2>& superGeometry,
                                        std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator,
                                        const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};
*/
/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */
/*
template <typename T, typename DESCRIPTOR>
class BlockLatticePorousMomentumLossForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  std::vector<SmoothIndicatorF2D<T,T,true>* >& _vectorOfIndicator;
public:
  BlockLatticePorousMomentumLossForce2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                        BlockGeometry<T,2>& blockGeometry,
                                        std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator,
                                        const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};
*/
}
#endif
