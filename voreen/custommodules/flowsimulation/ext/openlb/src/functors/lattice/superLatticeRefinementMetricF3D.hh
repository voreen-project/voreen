/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_REFINEMENT_METRIC_F_3D_HH
#define SUPER_LATTICE_REFINEMENT_METRIC_F_3D_HH

#include "superLatticeRefinementMetricF3D.h"
#include "blockLatticeRefinementMetricF3D.h"


namespace olb {


template<typename T, typename DESCRIPTOR>
SuperLatticeKnudsen3D<T, DESCRIPTOR>::SuperLatticeKnudsen3D(
  SuperLattice<T, DESCRIPTOR>& lattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(lattice, 1)
{
  this->getName() = "knudsen";

  for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeKnudsen3D<T, DESCRIPTOR>(lattice.getBlock(iC))
    );
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR>::SuperLatticeRefinementMetricKnudsen3D(
  SuperLattice<T, DESCRIPTOR>&      lattice,
  const UnitConverter<T, DESCRIPTOR>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(lattice, 1)
{
  this->getName() = "refinementMetricKnudsen";

  for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR>(
        lattice.getBlock(iC), converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR>::operator()(T output[], int glob)
{
  auto& load = this->_sLattice.getLoadBalancer();

  if (load.isLocal(glob)) {
    const int loc = load.loc(glob);

    return static_cast<BlockLatticeRefinementMetricKnudsen3D<T,DESCRIPTOR>&>(this->getBlockF(loc))(output);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR>::print()
{
  const int nC = this->_sLattice.getCuboidGeometry().getNc();

  std::vector<T> factors(nC, T{});
  T output[1] = { };

  auto& load = this->_sLattice.getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->getBlockF(iC)(output, iC);

    factors[load.glob(iC)] = output[0];
  }

  OstreamManager clout(std::cout, "refinement");

  for (int i = 0; i < nC; ++i) {
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().sendToMaster(&factors[i], 1, load.isLocal(i));
#endif

    clout << "factors[" << i << "]: " << factors[i] << std::endl;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeHighOrderKnudsen3D<T, DESCRIPTOR>::SuperLatticeHighOrderKnudsen3D(
  SuperLattice<T, DESCRIPTOR>& lattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(lattice, 1)
{
  this->getName() = "high_order_knudsen";

  for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeHighOrderKnudsen3D<T, DESCRIPTOR>(lattice.getBlock(iC))
    );
  }
}


}

#endif
