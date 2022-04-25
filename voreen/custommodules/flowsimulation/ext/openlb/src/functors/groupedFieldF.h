/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause,
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

#ifndef GROUPED_FIELD_F_H
#define GROUPED_FIELD_F_H


#include "functors/genericF.h"
#include "core/container.h"

namespace olb {

/**
 *  ContainerF is a NON-PARALLELIZED (no block/super differentiation) functor intended
 *  to extract data from Container objects as used e.g. in the particle framework.
 */
template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W=T>
class ContainerF : public GenericF<W,int> {
public:
  ContainerF( Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& container, int targetDim );
  
  Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& _container;

public:
  static constexpr bool isSuper = false;
  static constexpr unsigned d = DESCRIPTOR::d;

  /// \return container
  Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& getContainer();

  /// \return size of container 
  int getContainerSize() const;

  bool operator() (W output[], const int input []);

  using GenericF<W,int>::operator();
};



/**
 *  GroupedFieldF is a NON-PARALLELIZED (no block/super differentiation) functor.
 *  Motivated by enabling easy access to quantities in the particle system, this funcor
 *  can be used for any grouped field frameworks.
 *  When a proper parrallelization methodology for the particle framework exists,
 *  this might have to be extended.
 */

template <typename T, typename DESCRIPTOR, typename GROUP, typename FIELD>
class GroupedFieldF : public ContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>,T> {
public:
  GroupedFieldF( Container<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>>& container );
  
public:
  bool operator() (T output[], const int input []);

};



//TODO: Prototype edition for now!

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W=T>
class SuperContainerF : public GenericF<W,int> {
protected:
  SuperContainerF( LoadBalancer<T>& loadBalancer, int targetDim );

  LoadBalancer<T>& _loadBalancer;

  std::vector<std::unique_ptr<ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE>>> _containerF;
public:
  static constexpr bool isSuper = true;
  static constexpr unsigned d = DESCRIPTOR::d;

  LoadBalancer<T>& getLoadBalancer();

  int getContainerSize() const;
  
  ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& getContainerF(int iCloc);

  bool operator() (W output[], const int input []);

  using GenericF<W,int>::operator();
};


//Forward declaration
namespace particles{
template<typename T, typename DESCRIPTOR>
class SuperParticleSystem;
}

template <typename T, typename DESCRIPTOR, typename GROUP, typename FIELD, typename W=T>
class SuperParticleGroupedFieldF : public SuperContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>,W> { 
private:
  particles::SuperParticleSystem<T,DESCRIPTOR>& _sParticleSystem;
public:
  SuperParticleGroupedFieldF( particles::SuperParticleSystem<T,DESCRIPTOR>& sParticleSystem );
};




} //namespace olb

#endif
