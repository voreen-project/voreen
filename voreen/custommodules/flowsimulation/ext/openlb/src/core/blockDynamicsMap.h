/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef BLOCK_DYNAMICS_MAP_H
#define BLOCK_DYNAMICS_MAP_H

#include <vector>

#include "dynamics/dynamics.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
class BlockDynamicsMap {
private:
  const std::size_t _nCells;
  std::vector<Dynamics<T,DESCRIPTOR>*> _dynamics_of_cell;

public:
  BlockDynamicsMap(std::size_t nCells, Dynamics<T,DESCRIPTOR>* default_dynamics = &instances::getNoDynamics<T,DESCRIPTOR>()):
    _nCells(nCells),
    _dynamics_of_cell(nCells, default_dynamics)
  { }

  BlockDynamicsMap(Dynamics<T,DESCRIPTOR>* default_dynamics = &instances::getNoDynamics<T,DESCRIPTOR>()):
    BlockDynamicsMap(1, default_dynamics)
  { }

  const Dynamics<T,DESCRIPTOR>& get(std::size_t iCell) const
  {
    return *_dynamics_of_cell[iCell];
  }

  Dynamics<T,DESCRIPTOR>& get(std::size_t iCell)
  {
    return *_dynamics_of_cell[iCell];
  }

  void set(std::size_t iCell, Dynamics<T,DESCRIPTOR>* dynamics)
  {
    _dynamics_of_cell[iCell] = dynamics;
  }

};


}

#endif
