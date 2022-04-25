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

#ifndef SET_BOUNDARY_2D_H
#define SET_BOUNDARY_2D_H

#include "core/blockDynamicsMap.h"

namespace olb {

//sets boundary on indicated cells. This is a function, which can be used on many boundaries.
template<typename T, typename DESCRIPTOR>
void setBoundary(BlockLattice<T,DESCRIPTOR>& block, int iX, int iY,
                 Dynamics<T,DESCRIPTOR>* dynamics,
                 PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor)
{
  if (dynamics) {
    block.defineDynamics({iX, iY}, dynamics);
    auto cell = block.get(iX,iY);
    dynamics->initialize(cell);
  }
  if (postProcessor && !block.isPadding({iX,iY})) {
    block.addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR>
void setBoundary(BlockLattice<T,DESCRIPTOR>& block, int iX, int iY,
                 Dynamics<T,DESCRIPTOR>* dynamics)
{
  if (dynamics) {
    block.defineDynamics({iX, iY}, dynamics);
    auto cell = block.get(iX,iY);
    dynamics->initialize(cell);
  }
}

/// Adds needed Cells to the Communicator _commBC in SuperLattice
template<typename T, typename DESCRIPTOR>
void addPoints2CommBC(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,int _overlap)
{
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */

  if (_overlap == 0) {
    return;
  }

  auto& communicator = sLattice.getCommunicator(PostStream());
  communicator.template requestField<descriptors::POPULATION>();

  SuperGeometry<T,2>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        if (iX < 0 || iX > nX - 1 ||
            iY < 0 || iY > nY - 1 ) { // if within overlap
          if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY) != 0) {
            bool found = false;
            for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
              for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                const int nextX = iXo + iX;
                const int nextY = iYo + iY;
                if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY)) {
                  communicator.requestCell({iCloc, iX, iY});
                  found = true;
                }
              }
            }
          }
        }
      }
    }
  }

  communicator.exchangeRequests();
}


namespace boundaryhelper {

//instatiates DYNAMICS with normal values form discreteNormal Vector n
template <
  typename T, typename DESCRIPTOR,
  template<int...> typename DYNAMICS
>
DynamicsPromise<T,DESCRIPTOR> constructConcreteDynamicsForNormal(Vector<int,2> n)
{
  if (n == Vector<int,2> {1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,1>>{});
  }
  else if (n == Vector<int,2> {1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,-1>>{});
  }
  else if (n == Vector<int,2> {-1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<-1,1>>{});
  }
  else if (n == Vector<int,2> {-1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<-1,-1>>{});
  }
  else {
    throw std::runtime_error("Could not set Boundary.");
  }
}

//instatiates DYNAMICS with direction and orientation values from discreteNormal Vector n
template <
  typename T, typename DESCRIPTOR,
  template<int...> typename DYNAMICS
>
DynamicsPromise<T,DESCRIPTOR> constructConcreteDynamicsForDirectionOrientation(Vector<int,2> n)
{
  if (n[0] == 1) {
    return DynamicsPromise(meta::id<DYNAMICS<0,1>>{});
  }
  else if (n[0] == -1) {
    return DynamicsPromise(meta::id<DYNAMICS<0,-1>>{});
  }
  else if (n[1] == 1) {
    return DynamicsPromise(meta::id<DYNAMICS<1,1>>{});
  }
  else if (n[1] == -1) {
    return DynamicsPromise(meta::id<DYNAMICS<1,-1>>{});
  }
  else {
    throw std::runtime_error("Could not set Boundary.");
  }
}

//constructs DYNAMICS with a Momenta that expects a direction and orientation as template args
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename> typename DYNAMICS,
  template <int,int> typename MOMENTA
>
struct PlainDynamicsForDirectionOrientationMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MOMENTA<x,y>>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForDirectionOrientation<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with MixinDynamics and a Momenta that expects a direction and orientation as template args
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename> typename DYNAMICS,
  typename MIXIN,
  template <int,int> typename MOMENTA
>
struct PlainMixinDynamicsForDirectionOrientationMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA<x,y>>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForDirectionOrientation<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with MixinDynamics, direction, orientation and a Momenta that itself expects a direction and orientation
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename,int,int> typename DYNAMICS,
  typename MIXIN,
  template <int,int> typename MOMENTA
>
struct DirectionOrientationMixinDynamicsForDirectionOrientationMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA<x,y>,x,y>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForDirectionOrientation<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with MixinDynamics, Momenta, direction and orientation as template args
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename,int,int> typename DYNAMICS,
  typename MIXIN,
  typename MOMENTA
>
struct DirectionOrientationMixinDynamicsForPlainMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA,x,y>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForDirectionOrientation<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs MixinDynamics with a Momenta that expects direction and orientation as template args
template <
  typename T, typename DESCRIPTOR,
  typename MIXIN,
  template <int,int> typename MOMENTA
>
struct MixinDynamicsExchangeDirectionOrientationMomenta {
  template <int x, int y>
  using ConcreteDynamics = typename MIXIN::template exchange_momenta<MOMENTA<x,y>>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForDirectionOrientation<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with a Momenta that expects two normal values a template args
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename> typename DYNAMICS,
  template <int,int> typename MOMENTA
>
struct PlainDynamicsForNormalMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MOMENTA<x,y>>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForNormal<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with two normal values and a Momenta that itself expects two normal values
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,int,int> typename DYNAMICS,
  template <int,int> typename MOMENTA
>
struct NormalDynamicsForNormalMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MOMENTA<x,y>,x,y>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForNormal<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs DYNAMICS with template args MixinDynamics, two normal values and a momenta that expects two normal values
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename,int,int> typename DYNAMICS,
  typename MIXIN,
  template <int,int> typename MOMENTA
>
struct NormalMixinDynamicsForNormalMomenta {
  template <int x, int y>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA<x,y>,x,y>;

  static auto construct(Vector<int,2> n) {
    return constructConcreteDynamicsForNormal<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//instantiates TYPE, derived from RESULT with values form descreteNormal Vector n
//RESULT can be either Dynamics or PostProcessorGenerator2D
template <
  typename RESULT, typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE,
  typename... ARGS
>
RESULT* constructForNormal(Vector<int,2> n, ARGS&&... args)
{
  if (n == Vector<int,2> {1, 1}) {
    return new TYPE<T,DESCRIPTOR,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,2> {1, -1}) {
    return new TYPE<T,DESCRIPTOR,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,2> {-1, 1}) {
    return new TYPE<T,DESCRIPTOR,-1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,2> {-1, -1}) {
    return new TYPE<T,DESCRIPTOR,-1,-1>(std::forward<decltype(args)>(args)...);
  }
  else {
    throw std::runtime_error("Could not set Boundary.");
  }
}

template <
  typename RESULT, typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE
>
RESULT promiseForNormal(Vector<int,2> n)
{
  if (n == Vector<int,2> {1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,1>>();
  }
  else if (n == Vector<int,2> {1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,-1>>();
  }
  else if (n == Vector<int,2> {-1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,-1,1>>();
  }
  else if (n == Vector<int,2> {-1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,-1,-1>>();
  }
  else {
    throw std::runtime_error("Invalid normal");
  }
}

//constructs TYPE derived from PostProcessorGenerator2D with two normals as template args
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE,
  typename... ARGS
>
PostProcessorGenerator2D<T,DESCRIPTOR>* constructPostProcessorForNormal(Vector<int,2> n, ARGS&&... args)
{
  return constructForNormal<PostProcessorGenerator2D<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n, std::forward<decltype(args)>(args)...);
}

template <
  typename T, typename DESCRIPTOR,
  template<typename,typename,int,int> typename TYPE
>
PostProcessorPromise<T,DESCRIPTOR> promisePostProcessorForNormal(Vector<int,2> n)
{
  return promiseForNormal<PostProcessorPromise<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n);
}

//instantiates TYPE, derived from RESULT with values form descreteNormal Vector n
//RESULT can be either Dynamics or PostProcessorGenerator2D
template <
  typename RESULT, typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE,
  typename... ARGS
>
RESULT* constructForDirectionOrientation(Vector<int,2> n, ARGS&&... args)
{
  if (n[0] == 1) {
    return new TYPE<T,DESCRIPTOR,0,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n[0] == -1) {
    return new TYPE<T,DESCRIPTOR,0,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n[1] == 1) {
    return new TYPE<T,DESCRIPTOR,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n[1] == -1) {
    return new TYPE<T,DESCRIPTOR,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else {
    throw std::runtime_error("Could not set Boundary.");
  }
}

template <
  typename PROMISE,
  typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE
>
PROMISE promiseForDirectionOrientation(Vector<int,2> n)
{
  if (n[0] == 1) {
    return meta::id<TYPE<T,DESCRIPTOR,0,1>>();
  }
  else if (n[0] == -1) {
    return meta::id<TYPE<T,DESCRIPTOR,0,-1>>();
  }
  else if (n[1] == 1) {
    return meta::id<TYPE<T,DESCRIPTOR,1,1>>();
  }
  else if (n[1] == -1) {
    return meta::id<TYPE<T,DESCRIPTOR,1,-1>>();
  }
  else {
    throw std::runtime_error("Invalid normal");
  }
}

//constructs TYPE derived from PostProcessorGenerator2D with template args direction and orientation
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE,
  typename... ARGS
>
PostProcessorGenerator2D<T,DESCRIPTOR>* constructPostProcessorForDirectionOrientation(Vector<int,2> n, ARGS&&... args)
{
  return constructForDirectionOrientation<PostProcessorGenerator2D<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n, std::forward<decltype(args)>(args)...);
}

template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,int,int> typename TYPE
>
PostProcessorPromise<T,DESCRIPTOR> promisePostProcessorForDirectionOrientation(Vector<int,2> n)
{
  return promiseForDirectionOrientation<PostProcessorPromise<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n);
}

}//namespace boundaryhelper

}//namespace olb

#include "normalDynamicsContructors.h"

#endif
