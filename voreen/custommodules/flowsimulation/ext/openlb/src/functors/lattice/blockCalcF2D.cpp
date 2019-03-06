/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2017 Albert Mink, Lukas Baron, Mathias J. Krause,
 *                          Adrian Kummerländer
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

#include "blockCalcF2D.h"
#include "blockCalcF2D.hh"

namespace olb {

template class BlockCalc2D<int,util::plus>;
template class BlockCalc2D<double,util::plus>;
template class BlockCalc2D<bool,util::plus>;

template class BlockCalc2D<int,util::minus>;
template class BlockCalc2D<double,util::minus>;
template class BlockCalc2D<bool,util::minus>;

template class BlockCalc2D<int,util::multiplies>;
template class BlockCalc2D<double,util::multiplies>;
template class BlockCalc2D<bool,util::multiplies>;

template class BlockCalc2D<int,util::divides>;
template class BlockCalc2D<double,util::divides>;
template class BlockCalc2D<bool,util::divides>;

} // end namespace olb
