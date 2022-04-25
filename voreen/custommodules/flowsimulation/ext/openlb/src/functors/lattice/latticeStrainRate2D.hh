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

#ifndef LATTICE_STRAIN_RATE_2D_HH
#define LATTICE_STRAIN_RATE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeStrainRate2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

/*template <typename T, typename DESCRIPTOR>
BlockLatticeStrainRate2D<T,DESCRIPTOR>::BlockLatticeStrainRate2D
  (BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{ this->getName() = "strainRate"; }

template <typename T, typename DESCRIPTOR>
std::vector<T> BlockLatticeStrainRate2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[4];
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();

  strainRate[0] = -pi[0]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;
  strainRate[1] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;
  strainRate[2] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;
  strainRate[3] = -pi[2]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;

  //cout << pi[0] << " " << pi[1] << " " << pi[2] << " " << descriptors::invCs2<T,DESCRIPTOR>() << std::endl;

  std::vector<T> output(strainRate, strainRate+4); // first adress, last adress
  return output;
}*/

}
#endif
