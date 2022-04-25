/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod, Albert Mink, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef SUPER_STRUCTURE_HH
#define SUPER_STRUCTURE_HH

#include "communication/superStructure.h"

namespace olb {


template<typename T, unsigned D>
SuperStructure<T,D>::SuperStructure(CuboidGeometry<T,D>& cuboidGeometry,
                                    LoadBalancer<T>& loadBalancer,
                                    int overlap)
  : _cuboidGeometry(cuboidGeometry),
    _loadBalancer(loadBalancer),
    _overlap(overlap),
    clout(std::cout, "SuperGeometry" + std::to_string(D) + "D")
{
}

template<typename T, unsigned D>
SuperStructure<T,D>::SuperStructure(int overlap)
  : SuperStructure(*(new CuboidGeometry<T,D> ()),
                   *(new LoadBalancer<T> ()),
                   overlap)
{ }

template<typename T, unsigned D>
CuboidGeometry<T,D>& SuperStructure<T,D>::getCuboidGeometry()
{
  return _cuboidGeometry;
}

template<typename T, unsigned D>
CuboidGeometry<T,D> const& SuperStructure<T,D>::getCuboidGeometry() const
{
  return _cuboidGeometry;
}

template<typename T, unsigned D>
int SuperStructure<T,D>::getOverlap()
{
  return _overlap;
}

template<typename T, unsigned D>
int SuperStructure<T,D>::getOverlap() const
{
  return _overlap;
}

template<typename T, unsigned D>
LoadBalancer<T>& SuperStructure<T,D>::getLoadBalancer()
{
  return _loadBalancer;
}

template<typename T, unsigned D>
LoadBalancer<T> const& SuperStructure<T,D>::getLoadBalancer() const
{
  return _loadBalancer;
}

}

#endif
