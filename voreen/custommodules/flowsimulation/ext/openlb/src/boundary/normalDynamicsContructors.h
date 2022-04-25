/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Lennart Neukamm, Adrian Kummerlaender
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

#ifndef NORMAL_DYNAMICS_CONSTRUCTORS_H
#define NORMAL_DYNAMICS_CONSTRUCTORS_H

namespace olb {

namespace boundaryhelper {

//constructs DYNAMICS normal values and a plain momenta
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,int...> typename DYNAMICS,
  typename MOMENTA
>
struct NormalDynamicsForPlainMomenta {
  template <int... Normal>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MOMENTA,Normal...>;

  template <unsigned D>
  static auto construct(Vector<int,D> n) {
    return constructConcreteDynamicsForNormal<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with template args MixinDynamics, Momenta and normal values
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename,int...> typename DYNAMICS,
  typename MIXIN,
  typename MOMENTA
>
struct NormalMixinDynamicsForPlainMomenta {
  template <int... Normal>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA,Normal...>;

  template <unsigned D>
  static auto construct(Vector<int,D> n) {
    return constructConcreteDynamicsForNormal<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs Dynamics with template args Mixindynamics and Momenta which expects normal values
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename> typename DYNAMICS,
  typename MIXIN,
  template <int...> typename MOMENTA
>
struct PlainMixinDynamicsForNormalMomenta {
  template <int... Normal>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA<Normal...>>;

  template <unsigned D>
  static auto construct(Vector<int,D> n) {
    return constructConcreteDynamicsForNormal<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};


}//namespace boundaryhelper

}//namespace olb

#endif

