/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_DESCRIPTOR_ALIAS_H
#define PARTICLE_DESCRIPTOR_ALIAS_H


namespace olb {

namespace descriptors {

using ResolvedParticle3D = PARTICLE_DESCRIPTOR<3, GENERAL_TMP<3>,
      MOBILITY_VERLET<3>, SURFACE_RESOLVED<3>,
      FORCING_RESOLVED<3>, PHYSPROPERTIES_RESOLVED<3> >;

using ResolvedParticle2D = PARTICLE_DESCRIPTOR<2, GENERAL_TMP<2>,
      MOBILITY_VERLET<2>, SURFACE_RESOLVED<2>,
      FORCING_RESOLVED<2>, PHYSPROPERTIES_RESOLVED<2> >;

using ResolvedSphere3D = PARTICLE_DESCRIPTOR<3, GENERAL_TMP<3>,
      MOBILITY_VERLET<3>, SURFACE_RESOLVED_CIRCULAR<3>,
      FORCING_RESOLVED<3>, PHYSPROPERTIES_RESOLVED<3> >;

using ResolvedCircle2D = PARTICLE_DESCRIPTOR<2, GENERAL_TMP<2>,
      MOBILITY_VERLET<2>, SURFACE_RESOLVED_CIRCULAR<2>,
      FORCING_RESOLVED<2>, PHYSPROPERTIES_RESOLVED<2> >;

} //namespace descriptors

} //namespace olb


#endif
