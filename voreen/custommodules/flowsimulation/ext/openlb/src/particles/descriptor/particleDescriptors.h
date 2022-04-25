/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_DESCRIPTORS_H
#define PARTICLE_DESCRIPTORS_H

#include "utilities/dimensionConverter.h"
#include "particles/descriptor/particleDescriptorUtilities.h"

namespace olb {

namespace descriptors {

//Parent fields ensuring generalized access
struct ANG_VELOCITY {};
struct ANG_ACC_STRD {};
struct ANGLE {};
struct SINDICATOR {};
struct TORQUE {};
struct MOFI {};
struct ROT_MATRIX {};

//Common
struct POSITION          : public FIELD_BASE<0,  1, 0> { };
struct YOUNG_MODULUS     : public FIELD_BASE<1,  0, 0> { };
struct SHEAR_MODULUS     : public FIELD_BASE<1,  0, 0> { };
struct POISSON_RATIO     : public FIELD_BASE<1,  0, 0> { };
struct ADHESION          : public FIELD_BASE<2,  0, 0> { };
//Resolved
struct ACCELERATION_STRD : public FIELD_BASE<0,  1, 0> { };
//Subgrid
struct FORCE_STRD        : public FIELD_BASE<0,  1, 0> { };
struct ACTIVE            : public FIELD_BASE<1,  0, 0> { };
struct CUBOID            : public FIELD_BASE<1,  0, 0> { };
struct ID                : public FIELD_BASE<1,  0, 0> { };
struct MASS_ADDED        : public FIELD_BASE<1,  0, 0> { };
struct RADIUS            : public FIELD_BASE<1,  0, 0> { };

//Dimension sensitive fields
template<unsigned D>
struct ANG_VELOCITY_XD   : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public ANG_VELOCITY { };
template<unsigned D>
struct ANG_ACC_STRD_XD   : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public ANG_ACC_STRD { };
template<unsigned D>
struct ANGLE_XD          : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public ANGLE { };
template<unsigned D>
struct TORQUE_XD         : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public TORQUE { };
template<unsigned D>
struct MOFI_XD           : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public MOFI { };
template<unsigned D>
struct ROT_MATRIX_XD     : public FIELD_BASE<utilities::dimensions::convert<D>::matrix,  0, 0>, public ROT_MATRIX { };

//Indicator
template<unsigned D>
struct SINDICATOR_XD : public SINDICATOR, public FIELD_BASE<1> {
  template <typename T>
  using value_type = typename utilities::dimensions::convert<D>::template surfaceType<T>*;

  template <typename T>
  using column_type = AbstractColumn<typename utilities::dimensions::convert<D>::template surfaceType<T>*>;

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<1>>()>{};
  }
};

//Parent descriptors ensuring generalized access
struct GENERAL {};
struct MOBILITY {};
struct SURFACE {};
struct FORCING {};
struct PHYSPROPERTIES {};
/// Mechanical properties
struct MECHPROPERTIES {};
struct NUMERICPROPERTIES {};

//Field wrapping descriptors acting as GROUP
// -like the field list itself, this list can be extended by custom GROUP descriptors
// -in order to allow generalized access via getField<GROUP,FIELD> those should always inherit from the parent descriptors
template<unsigned D>
struct GENERAL_TMP : public PARTICLE_DESCRIPTOR<D,POSITION>, public GENERAL {};

template<unsigned D>
struct MOBILITY_VERLET : public PARTICLE_DESCRIPTOR<D,VELOCITY,ACCELERATION_STRD,ANG_VELOCITY_XD<D>,ANG_ACC_STRD_XD<D>>, public MOBILITY {};

template<unsigned D>
struct MOBILITY_VERLET_NO_ANGLE : public PARTICLE_DESCRIPTOR<D,VELOCITY,ACCELERATION_STRD>, public MOBILITY {};

template<unsigned D>
struct MOBILITY_EULER_NO_ANGLE : public PARTICLE_DESCRIPTOR<D,VELOCITY>, public MOBILITY {};

template<unsigned D>
struct SURFACE_RESOLVED : public PARTICLE_DESCRIPTOR<D,ANGLE_XD<D>,ROT_MATRIX_XD<D>,SINDICATOR_XD<D>>, public SURFACE {};

template<unsigned D>
struct SURFACE_RESOLVED_CIRCULAR : public PARTICLE_DESCRIPTOR<D,ANGLE_XD<D>,SINDICATOR_XD<D>>, public SURFACE {};

struct NUMERICPROPERTIES_SUBGRID : public PARTICLE_DESCRIPTOR<1,ACTIVE>, public NUMERICPROPERTIES {};

template<unsigned D>
struct FORCING_RESOLVED : public PARTICLE_DESCRIPTOR<D,FORCE,TORQUE_XD<D>>, public FORCING {};

template<unsigned D>
struct FORCING_ADHESIVE : public PARTICLE_DESCRIPTOR<D,FORCE,TORQUE_XD<D>,ADHESION>, public FORCING {};

template<unsigned D>
struct FORCING_SUBGRID : public PARTICLE_DESCRIPTOR<D,FORCE,FORCE_STRD>, public FORCING {};

template<unsigned D>
struct PHYSPROPERTIES_RESOLVED : public PARTICLE_DESCRIPTOR<D,MASS,MOFI_XD<D>>, public PHYSPROPERTIES {};

template<unsigned D>
struct PHYSPROPERTIES_SUBGRID : public PARTICLE_DESCRIPTOR<D,MASS,MASS_ADDED,RADIUS>, public PHYSPROPERTIES {};

} //namespace descriptors

} //namespace olb


#endif
