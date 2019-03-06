/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#include "momentaOnBoundaries.h"
#include "momentaOnBoundaries.hh"
#include "boundaryPostProcessors3D.h"
#include "boundaryPostProcessors3D.hh"
#include "core/postProcessing.h"
#include "core/postProcessing.hh"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/latticeDescriptors.hh"

namespace olb {

using namespace descriptors;

template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          0,1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 0,1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          0,-1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 0,-1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          1,1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 1,1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          1,-1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 1,-1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          2,1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 2,1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          2,-1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 2,-1>;

template class StraightConvectionBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          0,1>;
template class StraightConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 0,1>;
template class StraightConvectionBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          0,-1>;
template class StraightConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 0,-1>;
template class StraightConvectionBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          1,1>;
template class StraightConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 1,1>;
template class StraightConvectionBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          1,-1>;
template class StraightConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 1,-1>;
template class StraightConvectionBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          2,1>;
template class StraightConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 2,1>;
template class StraightConvectionBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          2,-1>;
template class StraightConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 2,-1>;

template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0,-1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0,-1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0,-1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0,-1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0, 1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0, 1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0, 1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0, 1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1,-1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1,-1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1,-1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1,-1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1, 1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1, 1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1, 1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1, 1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2,-1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2,-1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2,-1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2,-1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2, 1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2, 1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2, 1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2, 1, 1>;

template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1,-1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1,-1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1,-1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1,-1, 1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1, 1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1, 1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1, 1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1, 1, 1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1,-1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1,-1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1,-1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1,-1, 1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1, 1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1, 1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1, 1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1, 1, 1>;


template class SlipBoundaryProcessor3D<double, descriptors::D3Q19Descriptor>;
template class SlipBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor>;

}
