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
 * This object is a MRT LB dynamics as described in D.Yu et al. in
 * Progress in Aerospace Sciences 39 (2003) 329-367
 */
#ifndef MRT_DYNAMICS_H
#define MRT_DYNAMICS_H

#include "interface.h"
#include "collisionMRT.h"

namespace olb {

/**
 * Original implementation based on:
 * D'Humieres et al., "Multiple-relaxation-time lattice Boltzmann models in three dimensions",
 * Phil: Trans. R. soc. Lond. A (2002) 360, 437-451
 * and
 * Yu et al,, "LES of turbulent square jet flow using an MRT lattice Boltzmann model",
 * Computers & Fluids 35 (2006), 957-965
 **/
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using MRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::MRT
>;

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedMRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::MRT,
  forcing::LaddVerberg
>;

}

#endif
