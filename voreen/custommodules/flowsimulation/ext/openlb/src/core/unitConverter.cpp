/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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

#include "unitConverter.h"
#include "unitConverter.hh"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/latticeDescriptors.hh"
//#include "dynamics/advectionDiffusionLatticeDescriptors.h"
//#include "dynamics/advectionDiffusionLatticeDescriptors.hh"


namespace olb {

template class UnitConverter<double,descriptors::D2Q9Descriptor>;
template class UnitConverter<double,descriptors::D3Q19Descriptor>;
template class UnitConverter<double,descriptors::AdvectionDiffusionD3Q7Descriptor>;

template UnitConverter<double,descriptors::D2Q9Descriptor>* createUnitConverter(const olb::XMLreader&);

}  // namespace olb
