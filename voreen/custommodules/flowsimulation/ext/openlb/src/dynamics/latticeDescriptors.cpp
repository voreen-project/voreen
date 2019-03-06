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

/** \file
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- template instantiation
 */
#include "latticeDescriptors.h"
#include "latticeDescriptors.hh"

namespace olb {
namespace descriptors {
//template class DescriptorBase<1>;
//template class DescriptorBase<2>;
//template class DescriptorBase<3>;
//template class DescriptorBase<9>;
//template class DescriptorBase<19>;

template struct D2Q5DescriptorBase<double>;
template struct AdvectionDiffusionD2Q5Descriptor<double>;
template struct AdvectionDiffusionMRTD2Q5DescriptorBase<double>;
template struct AdvectionDiffusionMRTD2Q5Descriptor<double>;
template struct D2Q9DescriptorBase<double>;
template struct D2Q9Descriptor<double>;
template struct D3Q7DescriptorBase<double>;
template struct AdvectionDiffusionD3Q7Descriptor<double>;
template struct AdvectionDiffusionMRTD3Q7DescriptorBase<double>;
template struct AdvectionDiffusionMRTD3Q7Descriptor<double>;
template struct D3Q13DescriptorBase<double>;
template struct D3Q15DescriptorBase<double>;
template struct D3Q19DescriptorBase<double>;
template struct D3Q19Descriptor<double>;
template struct D3Q27DescriptorBase<double>;
}
}

