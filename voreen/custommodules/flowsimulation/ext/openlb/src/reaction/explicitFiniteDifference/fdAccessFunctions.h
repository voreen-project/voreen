/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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
 * Access functions to access the correct element of the finite-difference external field.
 *  -- header file
 */
#ifndef FD_ACCESS_FUNCTIONS_H
#define FD_ACCESS_FUNCTIONS_H

namespace olb {

namespace fd {

template<typename T, typename DESCRIPTOR, typename FIELD>
T* accessOld(Cell<T,DESCRIPTOR> cell, std::size_t iT)
{
  return &cell.template getFieldPointer<FIELD>()[iT % 2];
}

template<typename T, typename DESCRIPTOR, typename FIELD>
T* accessNew(Cell<T,DESCRIPTOR> cell, std::size_t iT)
{
  return &cell.template getFieldPointer<FIELD>()[(iT+1) % 2];
}

} // fd

} // olb

#endif
