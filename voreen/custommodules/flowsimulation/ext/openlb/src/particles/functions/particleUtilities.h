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



#ifndef PARTICLE_UTILITIES_H
#define PARTICLE_UTILITIES_H

#include <cassert>

namespace olb {

namespace particles {

namespace sorting {


//Sort particleSystem by provided nested fields
template<typename T, typename PARTICLETYPE, typename ...NESTED_FIELDS>
std::size_t partitionParticleSystem( ParticleSystem<T,PARTICLETYPE>& particleSystem ){
  using namespace descriptors;
  //Find first occurence of element belonging to partition B
  std::size_t iPfirstB=particleSystem.size(); //Fallback (simple to check)
  bool valB = true;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
    auto particle = particleSystem.get(iP);
    bool active = particle.template getField<NESTED_FIELDS...>();
    if (active==valB){
      iPfirstB=iP; 
      break;
    }
  }
  //Find succeeding elements belonging to partition A and move them to othe beginning
  for (std::size_t iP=iPfirstB; iP<particleSystem.size(); ++iP){
    auto particle = particleSystem.get(iP);
    bool active = particle.template getField<NESTED_FIELDS...>();
    if (active!=valB){
      particleSystem.swapParticles(iP,iPfirstB);
      ++iPfirstB;
    } 
  }
  //Return splitpoint between partitions
  return iPfirstB;
}

} //namespace sorting

} //namespace particles

} //namespace olb


#endif
