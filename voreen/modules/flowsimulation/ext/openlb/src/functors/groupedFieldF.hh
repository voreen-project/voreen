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

#ifndef GROUPED_FIELD_F_HH
#define GROUPED_FIELD_F_HH

#include "functors/groupedFieldF.h"

namespace olb {


template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::ContainerF( 
  Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& container, int targetDim )
  : GenericF<W,int>( targetDim, 1), _container(container)
{
  this->getName() = "ContainerF";
}

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& 
  ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::getContainer()
{
  return _container;
}

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
int ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::getContainerSize() const
{
  OLB_ASSERT(_container.size() < INT32_MAX,
             "cast from std::size_t to int unsafe");
  return _container.size();
}

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
bool ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::operator()(W output[], const int input[])
{
  return true; //TODO: add Error handling using true/false
}


template <typename T, typename DESCRIPTOR, typename GROUP, typename FIELD>
GroupedFieldF<T,DESCRIPTOR,GROUP,FIELD>::GroupedFieldF( 
  Container<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>>& container )
  : ContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>,T>( container,
      DESCRIPTOR::template size<typename meta::derived_type_in_nested<DESCRIPTOR,GROUP,FIELD>::type>() )
{
  this->getName() = std::string(typeid(FIELD).name());
}

template <typename T, typename DESCRIPTOR, typename GROUP, typename FIELD>
bool GroupedFieldF<T,DESCRIPTOR,GROUP,FIELD>::operator()(T output[], const int input[])
{
  using namespace olb::descriptors;
  using GROUP_EVAL = typename DESCRIPTOR::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL::template derivedField<FIELD>;
  auto fieldValue = this->getContainer().data().template get<GROUP_EVAL>()
                  .template getFieldPointer<FIELD_EVAL>(input[0]);

  for (int iDim=0; iDim<this->getTargetDim(); ++iDim){
    output[iDim] = fieldValue[iDim]; 
  }
  return true; //TODO: add Error handling using true/false
}




//TODO: Prototype edition for now!
//

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
SuperContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::SuperContainerF(LoadBalancer<T>& loadBalancer, int targetDim)
  : GenericF<W,int>(targetDim,1), _loadBalancer(loadBalancer)
{ }

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
LoadBalancer<T>& SuperContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::getLoadBalancer()
{
  return _loadBalancer;
}

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
int SuperContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::getContainerSize() const
{
  OLB_ASSERT(_containerF.size() < INT32_MAX,
             "cast from std::size_t to int unsafe");
  return _containerF.size();
}

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
ContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE>& SuperContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::getContainerF(int iCloc)
{
  OLB_ASSERT(iCloc < int(_containerF.size()) && iCloc >= 0,
             "block functor index outside bounds");
  return *(_containerF[iCloc]);
}

template <typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE, typename W>
bool SuperContainerF<T,DESCRIPTOR,FIELD_ARRAY_TYPE,W>::operator()(W output[], const int input[])
{
  if (_loadBalancer.isLocal(input[0])) {
    const int loc = _loadBalancer.loc(input[0]);
    return this->getContainerF(loc)(output, &input[1]);
  }
  else {
    return false;
  }
}


template <typename T, typename DESCRIPTOR, typename GROUP, typename FIELD, typename W>
SuperParticleGroupedFieldF<T,DESCRIPTOR,GROUP,FIELD,W>::SuperParticleGroupedFieldF( 
  particles::SuperParticleSystem<T,DESCRIPTOR>& sParticleSystem )
  : SuperContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>,W>(
      sParticleSystem.getLoadBalancer(),
      DESCRIPTOR::template size<typename meta::derived_type_in_nested<DESCRIPTOR,GROUP,FIELD>::type>() ),
    _sParticleSystem(sParticleSystem)
{
  this->getName() = std::string(typeid(FIELD).name());

  int maxC = this->_sParticleSystem.getLoadBalancer().size();
  this->_containerF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {

    auto bParticleSystems = _sParticleSystem.getBlockParticleSystems();
    auto& particleSystem = *bParticleSystems[iC];
    auto& container = particleSystem.get();

    this->_containerF.emplace_back(new GroupedFieldF<T,DESCRIPTOR,GROUP,FIELD>( container ) );
  }
}





} //namespace olb

#endif
