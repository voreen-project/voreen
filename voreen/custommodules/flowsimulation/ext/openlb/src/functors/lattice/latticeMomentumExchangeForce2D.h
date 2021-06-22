/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Robin Trunk
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

#ifndef LATTICE_MOMENTUM_EXCHANGE_FORCE_2D_H
#define LATTICE_MOMENTUM_EXCHANGE_FORCE_2D_H

#include<vector>

#include "superBaseF2D.h"
#include "superCalcF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "core/superLattice2D.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "indicator/blockIndicatorBaseF2D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "dynamics/porousBGKdynamics.h"


namespace olb {

/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticeMomentumExchangeForce2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticeMomentumExchangeForce2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry2D<T>& superGeometry,
                                      std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter,
                                      Vector<bool,2> periodicity = Vector<bool,2> (false,false) );
  bool operator() (T output[], const int input[]) override;
};

/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */

template <typename T, typename DESCRIPTOR>
class BlockLatticeMomentumExchangeForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  std::vector<SmoothIndicatorF2D<T,T,true>* >& _vectorOfIndicator;
  Vector<T,2> _cellMin;
  Vector<T,2> _cellMax;
  Vector<bool,2> _periodic;
public:
  BlockLatticeMomentumExchangeForce2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                      BlockGeometryStructure2D<T>& blockGeometry,
                                      std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter,
                                      Vector<T,2> cellMin = Vector<T,2> (0.,0.),
                                      Vector<T,2> cellMax = Vector<T,2> (0.,0.),
                                      Vector<bool,2> periodic = Vector<bool,2> (false,false) );
  void evaluate(T output[], SmoothIndicatorF2D<T,T,true>& sIndicator, int iInd);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
