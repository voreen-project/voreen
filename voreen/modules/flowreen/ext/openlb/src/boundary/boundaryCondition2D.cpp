/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#include "boundaryCondition2D.h"
#include "boundaryCondition2D.hh"
#include "momentaOnBoundaries.h"
#include "momentaOnBoundaries.hh"
#include "momentaOnBoundaries2D.h"
#include "momentaOnBoundaries2D.hh"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/latticeDescriptors.hh"

namespace olb {

template class BoundaryConditionInstantiator2D<double, descriptors::D2Q9Descriptor,
    RegularizedBoundaryManager2D < double, descriptors::D2Q9Descriptor,
                                   RLBdynamics<double,descriptors::D2Q9Descriptor> > >;

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9Descriptor>*
createLocalBoundaryCondition2D < double,descriptors::D2Q9Descriptor,
                               RLBdynamics<double,descriptors::D2Q9Descriptor> >
                               (
                                 BlockLatticeStructure2D<double,descriptors::D2Q9Descriptor>& block
                               );

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9Descriptor>*
createLocalBoundaryCondition2D < double,descriptors::D2Q9Descriptor,
                               BGKdynamics<double,descriptors::D2Q9Descriptor> >
                               (
                                 BlockLatticeStructure2D<double,descriptors::D2Q9Descriptor>& block
                               );

template class BoundaryConditionInstantiator2D<double, descriptors::D2Q9Descriptor,
    InterpolationBoundaryManager2D < double, descriptors::D2Q9Descriptor,
                                     BGKdynamics<double,descriptors::D2Q9Descriptor> > >;

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9Descriptor>*
createInterpBoundaryCondition2D < double,descriptors::D2Q9Descriptor,
                                BGKdynamics<double,descriptors::D2Q9Descriptor> >
                                (
                                  BlockLatticeStructure2D<double,descriptors::D2Q9Descriptor>& block
                                );

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9Descriptor>*
createInterpBoundaryCondition2D < double,descriptors::D2Q9Descriptor,
                                ConstRhoBGKdynamics<double,descriptors::D2Q9Descriptor> >
                                (
                                  BlockLatticeStructure2D<double,descriptors::D2Q9Descriptor>& block
                                );
}
