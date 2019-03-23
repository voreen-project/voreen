/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#include "blockLatticeLocalF2D.h"
#include "blockLatticeLocalF2D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class BlockLatticeDissipation2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysDissipation2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeDensity2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeVelocity2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysStrainRate2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysWallShearStress2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeGeometry2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeRank2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeCuboid2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysPressure2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysVelocity2D<double,descriptors::D2Q9Descriptor>;
//template class BlockLatticePhysExternalPorosity2D<double,descriptors::D2Q9Descriptor>;
//template class BlockLatticePhysExternalVelocity2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysBoundaryForce2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysCorrBoundaryForce2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePorosity2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysPermeability2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticePhysDarcyForce2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeAverage2D<double,descriptors::D2Q9Descriptor>;
template class BlockEuklidNorm2D<double,descriptors::D2Q9Descriptor>;
template class BlockLatticeMomentumExchangeForce2D<double,descriptors::D2Q9Descriptor>;

} // end namespace olb
