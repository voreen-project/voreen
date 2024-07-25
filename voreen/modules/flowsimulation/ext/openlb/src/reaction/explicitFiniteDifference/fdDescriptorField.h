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
 * Descriptor field to store the advection-diffusion field(s):
 * element[0] -> old timestep; element[1] -> current timestep
 *  -- header file
 */
#ifndef FD_DESCRIPTOR_FIELD_H
#define FD_DESCRIPTOR_FIELD_H


namespace olb {

namespace descriptors {

// Field types need to be distinct (i.e. not aliases)
//                                   Cs Ds Qs
struct AD_FIELD  : public FIELD_BASE<2, 0, 0> { };
struct AD_SOURCE : public FIELD_BASE<1, 0, 0> { };

} // descriptors

} // olb

#endif
