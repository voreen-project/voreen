/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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

#include "offBoundaryCondition2D.h"
#include "offBoundaryCondition2D.hh"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/latticeDescriptors.hh"

namespace olb {

template class OffBoundaryConditionInstantiator2D
<double, descriptors::D2Q9Descriptor, BouzidiBoundaryManager2D < double, descriptors::D2Q9Descriptor, BGKdynamics<double,descriptors::D2Q9Descriptor> > >;

template class OffBoundaryConditionInstantiator2D
<double, descriptors::D2Q9Descriptor, BounceBackBoundaryManager2D < double, descriptors::D2Q9Descriptor> >;

template OffLatticeBoundaryCondition2D<double,descriptors::D2Q9Descriptor>*
createBouzidiBoundaryCondition2D < double,descriptors::D2Q9Descriptor,
                                 BGKdynamics<double,descriptors::D2Q9Descriptor> >
                                 (
                                   BlockLatticeStructure2D<double,descriptors::D2Q9Descriptor>& block
                                 );

template OffLatticeBoundaryCondition2D<double,descriptors::D2Q9Descriptor>*
createBounceBackBoundaryCondition2D < double,descriptors::D2Q9Descriptor>
(BlockLatticeStructure2D<double,descriptors::D2Q9Descriptor>& block);
}

