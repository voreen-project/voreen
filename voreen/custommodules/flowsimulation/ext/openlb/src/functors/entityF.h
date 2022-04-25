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

#ifndef ENTITY_F_H
#define ENTITY_F_H

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template <typename T, unsigned D, typename ENTITY>
class EntityF : public GenericF<T,int> {
protected:
  EntityF(int targetDim, int sourceDim, ENTITY& entity);
  ENTITY& _entity;
public:
  static constexpr int d = D;
  static constexpr bool isSuper = false;
  /// virtual destructor for defined behaviour
  ~EntityF() override {};
  /// get extent
  auto getExtent(){ return _entity.getExtent();}
  /// expose entity reference
  ENTITY& getEntityF(){ return _entity; }
  /// operator
  using GenericF<T,int>::operator();
};


template <typename T, unsigned D, typename ENTITY, typename W=T>
class SuperEntityF : public GenericF<W,int> {
protected:
  SuperEntityF(LoadBalancer<T>& loadBalancer, int targetDim, int sourceDim);
  LoadBalancer<T>& _loadBalancer;
  /**
   * By convention: If entity level functors are used at all, they should
   * number exactly LoadBalancer<T>::size per process.
   **/
  std::vector<std::unique_ptr<EntityF<W,D,ENTITY>>> _entityF;
public:
  static constexpr int d = D;
  static constexpr bool isSuper = true;
  /// expose load balancer
  LoadBalancer<T>& getLoadBalancer();
  /// return number of entity functors
  int getNumOfEntityF() const;
  /// get specific entity functor
  EntityF<W,D,ENTITY>& getEntityF(int iCloc);
  /// operator
  bool operator() (W output[], const int input []);
  using GenericF<W,int>::operator();
};

} //namespace olb

#endif
