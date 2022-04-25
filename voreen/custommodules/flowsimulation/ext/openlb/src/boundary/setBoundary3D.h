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

#ifndef SET_BOUNDARY_3D_H
#define SET_BOUNDARY_3D_H

namespace olb {

///sets boundary on indicated cells. This is a function, which can be used on many boundaries.
template<typename T, typename DESCRIPTOR >
void setBoundary(BlockLattice<T,DESCRIPTOR>& _block, int iX,int iY,int iZ,
                 Dynamics<T,DESCRIPTOR>* dynamics,PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor)
{
  if (dynamics) {
    _block.defineDynamics({iX, iY, iZ}, dynamics);
    auto cell = _block.get(iX,iY,iZ);
    dynamics->initialize(cell);
  }
  if (postProcessor && !_block.isPadding({iX,iY,iZ})) {
    _block.addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR >
void setBoundary(BlockLattice<T,DESCRIPTOR>& _block, int iX,int iY,int iZ,
                 Dynamics<T,DESCRIPTOR>* dynamics)
{
  if (dynamics) {
    _block.defineDynamics({iX, iY, iZ}, dynamics);
    auto cell = _block.get(iX,iY,iZ);
    dynamics->initialize(cell);
  }
}

/// Adds needed Cells to the Communicator _commBC in SuperLattice
template<typename T, typename DESCRIPTOR>
void addPoints2CommBC(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator, int _overlap)
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

  SuperGeometry<T,3>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();
    const int nZ = superGeometry.getBlockGeometry(iCloc).getNz();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        for (int iZ = -_overlap; iZ < nZ+_overlap; ++iZ) {
          if (iX < 0 || iX > nX - 1 ||
              iY < 0 || iY > nY - 1 ||
              iZ < 0 || iZ > nZ - 1 ) { // if within overlap
            if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ) != 0) {
              bool found = false;
              for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
                for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                  for (int iZo = -_overlap; iZo <= _overlap && !found; ++iZo) {
                    const int nextX = iXo + iX;
                    const int nextY = iYo + iY;
                    const int nextZ = iZo + iZ;
                    if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY, nextZ)
                        && nextX >= -_overlap && nextX < nX+_overlap
                        && nextY >= -_overlap && nextY < nY+_overlap
                        && nextZ >= -_overlap && nextZ < nZ+_overlap) {
                      communicator.requestCell({iCloc, iX, iY, iZ});
                      found = true;
                    }
                  }
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

//Instatiates DYNAMICS with three normal values from discreteNormal vector n
template <
  typename T, typename DESCRIPTOR,
  template <int...> typename DYNAMICS
>
DynamicsPromise<T,DESCRIPTOR> constructConcreteDynamicsForNormal(Vector<int,3> n)
{
  if (n == Vector<int,3> {1, 1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,1,1>>{});
  }
  else if (n == Vector<int,3> {1, -1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,-1,1>>{});
  }
  else if (n == Vector<int,3> {1, 1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,1,-1>>{});
  }
  else if (n == Vector<int,3> {1, -1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,-1,-1>>{});
  }
  else if (n == Vector<int,3> {-1, 1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<-1,1,1>>{});
  }
  else if (n == Vector<int,3> {-1, -1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<-1,-1,1>>{});
  }
  else if (n == Vector<int,3> {-1, 1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<-1,1,-1>>{});
  }
  else if (n == Vector<int,3> {-1, -1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<-1,-1,-1>>{});
  }
  else {
    throw std::runtime_error("Could not set Boundary.");
  }
}

//Instatiates DYNAMICS with a plane and two normal values from discreteNormal vector n
template <
  typename T, typename DESCRIPTOR,
  template <int...> typename DYNAMICS
>
DynamicsPromise<T,DESCRIPTOR> constructConcreteDynamicsForNormalSpecial(Vector<int,3> n)
{
  if (n == Vector<int,3> {0, 1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<0,1,1>>{});
  }
  else if (n == Vector<int,3> {0, -1, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<0,-1,1>>{});
  }
  else if (n == Vector<int,3> {0, 1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<0,1,-1>>{});
  }
  else if (n == Vector<int,3> {0, -1, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<0,-1,-1>>{});
  }
  else if (n == Vector<int,3> {1, 0, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,1,1>>{});
  }
  else if (n == Vector<int,3> {-1, 0, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,-1,-1>>{});
  }
  else if (n == Vector<int,3> {-1, 0, 1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,1,-1>>{});
  }
  else if (n == Vector<int,3> {1, 0, -1}) {
    return DynamicsPromise(meta::id<DYNAMICS<1,-1,1>>{});
  }
  else if (n == Vector<int,3> {1, 1, 0}) {
    return DynamicsPromise(meta::id<DYNAMICS<2,1,1>>{});
  }
  else if (n == Vector<int,3> {-1, 1, 0}) {
    return DynamicsPromise(meta::id<DYNAMICS<2,-1,1>>{});
  }
  else if (n == Vector<int,3> {1, -1, 0}) {
    return DynamicsPromise(meta::id<DYNAMICS<2,1,-1>>{});
  }
  else if (n == Vector<int,3> {-1, -1, 0}) {
    return DynamicsPromise(meta::id<DYNAMICS<2,-1,-1>>{});
  }
  else {
    throw std::invalid_argument("Invalid normal");
  }
}

//constructs Dynamics with a Momenta, a plane and two normals as template args
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,int,int,int> typename DYNAMICS,
  typename MOMENTA
>
struct NormalSpecialDynamicsForPlainMomenta {
  template <int x, int y, int z>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MOMENTA,x,y,z>;

  static auto construct(Vector<int,3> n) {
    return constructConcreteDynamicsForNormalSpecial<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs Dynamics with template args MixinDynamics and a Momenta that expects a plane and two normals
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename> typename DYNAMICS,
  typename MIXIN,
  template <int...> typename MOMENTA
>
struct PlainMixinDynamicsForNormalSpecialMomenta {
  template <int... PlaneAndNormal>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA<PlaneAndNormal...>>;

  static auto construct(Vector<int,3> n) {
    return constructConcreteDynamicsForNormalSpecial<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};

//constructs Dynamics with template args MixinDynamics, Momenta, plane normal1 and normal2
template <
  typename T, typename DESCRIPTOR,
  template <typename,typename,typename,typename,int,int,int> typename DYNAMICS,
  typename MIXIN,
  typename MOMENTA
>
struct NormalSpecialMixinDynamicsForPlainMomenta {
  template <int x, int y, int z>
  using ConcreteDynamics = DYNAMICS<T,DESCRIPTOR,MIXIN,MOMENTA,x,y,z>;

  static auto construct(Vector<int,3> n) {
    return constructConcreteDynamicsForNormalSpecial<T,DESCRIPTOR,ConcreteDynamics>(n);
  }
};


//Instantiates TYPE derived from RESULT with values from descreteNormal Vector n
//RESULT can be either Dynamics or PostProcessorGenerator3D
template <
  typename RESULT, typename T, typename DESCRIPTOR,
  template <typename,typename,int,int,int> typename TYPE,
  typename... ARGS
>
RESULT* constructForNormal(Vector<int,3> n, ARGS&&... args)
{
  if (n == Vector<int,3> {1, 1, 1}) {
    return new TYPE<T,DESCRIPTOR,1,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, -1, 1}) {
    return new TYPE<T,DESCRIPTOR,1,-1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, 1, -1}) {
    return new TYPE<T,DESCRIPTOR,1,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, -1, -1}) {
    return new TYPE<T,DESCRIPTOR,1,-1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, 1, 1}) {
    return new TYPE<T,DESCRIPTOR,-1,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, -1, 1}) {
    return new TYPE<T,DESCRIPTOR,-1,-1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, 1, -1}) {
    return new TYPE<T,DESCRIPTOR,-1,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, -1, -1}) {
    return new TYPE<T,DESCRIPTOR,-1,-1,-1>(std::forward<decltype(args)>(args)...);
  }
  else {
    return nullptr;
  }
}

template <
  typename RESULT, typename T, typename DESCRIPTOR,
  template <typename,typename,int,int,int> typename TYPE
>
RESULT promiseForNormal(Vector<int,3> n)
{
  if (n == Vector<int,3> {1, 1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,1,1>>();
  }
  else if (n == Vector<int,3> {1, -1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,-1,1>>();
  }
  else if (n == Vector<int,3> {1, 1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,1,-1>>();
  }
  else if (n == Vector<int,3> {1, -1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,-1,-1>>();
  }
  else if (n == Vector<int,3> {-1, 1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,-1,1,1>>();
  }
  else if (n == Vector<int,3> {-1, -1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,-1,-1,1>>();
  }
  else if (n == Vector<int,3> {-1, 1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,-1,1,-1>>();
  }
  else if (n == Vector<int,3> {-1, -1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,-1,-1,-1>>();
  }
  else {
    throw std::domain_error("Invalid normal");
  }
}

//Instantiates TYPE derived from RESULT with values from descreteNormal Vector n
//RESULT can be either Dynamics or PostProcessorGenerator3D
template <
  typename RESULT, typename T, typename DESCRIPTOR,
  template <typename,typename,int,int,int> typename TYPE,
  typename... ARGS
>
RESULT* constructForNormalSpecial(Vector<int,3> n, ARGS&&... args)
{
  if (n == Vector<int,3> {0, 1, 1}) {
    return new TYPE<T,DESCRIPTOR,0,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {0, -1, 1}) {
    return new TYPE<T,DESCRIPTOR,0,-1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {0, 1, -1}) {
    return new TYPE<T,DESCRIPTOR,0,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {0, -1, -1}) {
    return new TYPE<T,DESCRIPTOR,0,-1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, 0, 1}) {
    return new TYPE<T,DESCRIPTOR,1,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, 0, -1}) {
    return new TYPE<T,DESCRIPTOR,1,-1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, 0, 1}) {
    return new TYPE<T,DESCRIPTOR,1,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, 0, -1}) {
    return new TYPE<T,DESCRIPTOR,1,-1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, 1, 0}) {
    return new TYPE<T,DESCRIPTOR,2,1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, 1, 0}) {
    return new TYPE<T,DESCRIPTOR,2,-1,1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {1, -1, 0}) {
    return new TYPE<T,DESCRIPTOR,2,1,-1>(std::forward<decltype(args)>(args)...);
  }
  else if (n == Vector<int,3> {-1, -1, 0}) {
    return new TYPE<T,DESCRIPTOR,2,-1,-1>(std::forward<decltype(args)>(args)...);
  }
  else {
    return nullptr;
  }
}

template <
  typename PROMISE,
  typename T, typename DESCRIPTOR,
  template <typename,typename,int,int,int> typename TYPE
>
PROMISE promiseForNormalSpecial(Vector<int,3> n)
{
  if (n == Vector<int,3> {0, 1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,0,1,1>>();
  }
  else if (n == Vector<int,3> {0, -1, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,0,-1,1>>();
  }
  else if (n == Vector<int,3> {0, 1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,0,1,-1>>();
  }
  else if (n == Vector<int,3> {0, -1, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,0,-1,-1>>();
  }
  else if (n == Vector<int,3> {1, 0, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,1,1>>();
  }
  else if (n == Vector<int,3> {-1, 0, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,-1,-1>>();
  }
  else if (n == Vector<int,3> {-1, 0, 1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,1,-1>>();
  }
  else if (n == Vector<int,3> {1, 0, -1}) {
    return meta::id<TYPE<T,DESCRIPTOR,1,-1,1>>();
  }
  else if (n == Vector<int,3> {1, 1, 0}) {
    return meta::id<TYPE<T,DESCRIPTOR,2,1,1>>();
  }
  else if (n == Vector<int,3> {-1, 1, 0}) {
    return meta::id<TYPE<T,DESCRIPTOR,2,-1,1>>();
  }
  else if (n == Vector<int,3> {1, -1, 0}) {
    return meta::id<TYPE<T,DESCRIPTOR,2,1,-1>>();
  }
  else if (n == Vector<int,3> {-1, -1, 0}) {
    return meta::id<TYPE<T,DESCRIPTOR,2,-1,-1>>();
  }
  else {
    throw std::domain_error("Invalid normal");
  }
}

//constructs TYPE derived from PostProcessorGenerator3D with three normals as template args
template <
  typename T, typename DESCRIPTOR,
  template<typename,typename,int,int,int> typename TYPE,
  typename... ARGS
>
PostProcessorGenerator3D<T,DESCRIPTOR>* constructPostProcessorForNormal(Vector<int,3> n, ARGS&&... args)
{
  return constructForNormal<PostProcessorGenerator3D<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n, std::forward<decltype(args)>(args)...);
}

template <
  typename T, typename DESCRIPTOR,
  template<typename,typename,int,int,int> typename TYPE
>
PostProcessorPromise<T,DESCRIPTOR> promisePostProcessorForNormal(Vector<int,3> n)
{
  return promiseForNormal<PostProcessorPromise<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n);
}

//constructs TYPE derived from PostProcessorGenerator3D with a plane and two normals as template args
template <
  typename T, typename DESCRIPTOR,
  template<typename,typename,int,int,int> typename TYPE,
  typename... ARGS
>
PostProcessorGenerator3D<T,DESCRIPTOR>* constructPostProcessorForNormalSpecial(Vector<int,3> n, ARGS&&... args)
{
  return constructForNormalSpecial<PostProcessorGenerator3D<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n, std::forward<decltype(args)>(args)...);
}

template <
  typename T, typename DESCRIPTOR,
  template<typename,typename,int,int,int> typename TYPE
>
PostProcessorPromise<T,DESCRIPTOR> promisePostProcessorForNormalSpecial(Vector<int,3> n)
{
  return promiseForNormalSpecial<PostProcessorPromise<T,DESCRIPTOR>,T,DESCRIPTOR,TYPE>(n);
}

}//namespace boundaryhelper

}//namespace olb

#include "normalDynamicsContructors.h"

#endif
