/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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


/* This file contains particle creator functions.
 * Those are generally meant for users to be able to create particles
 * in a comfortable way. As this section might be subject to frequent changes,
 * those should, however, be avoided inside header files. Here, direct field access
 * as shown below should be preferred!
 *
*/


#ifndef PARTICLE_CREATOR_FUNCTIONS_2D_H
#define PARTICLE_CREATOR_FUNCTIONS_2D_H


#include "particles/particleSystem.h"
#include "particles/functions/particleCreatorHelperFunctions.h"


namespace olb {

namespace particles {

namespace creators {

template<typename T, typename PARTICLETYPE, bool ROTATION_IS_OPTIONAL=false>
void setResolvedObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle, std::size_t idxSurface,
  const Vector<T,2>& position, T density, T angle,
  const Vector<T,2>& velocity)
{
  setResolvedObject<T,PARTICLETYPE,ROTATION_IS_OPTIONAL>( particleSystem, idxParticle, idxSurface, position, density,
      Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>( angle ), velocity );
}

/// Add resolved object as new particle with existing surface
template<typename T, typename PARTICLETYPE>
void addResolvedObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem, std::size_t idxSurface,
  const Vector<T,2>& position, T density=0., T angle=0.,
  const Vector<T,2>& velocity = Vector<T,2> (0.))
{
  /// Set resolved object 2D at given index
  addResolvedObject( particleSystem, idxSurface, position, density,
                     Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>( angle ), velocity );
}


//// CIRCLE 2D

/// Set resolved circle for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCircle2D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  Vector<T,2> position, T radius, T epsilon, T density=0.,
  Vector<T,2> velocity = Vector<T,2>(0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF2D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCircle2D<T, T, true>(Vector<T,2>(0.), radius, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  // Safety mechanism for wrong particle type - It is better not to use rotation matrix!
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>() ) {
    OstreamManager clout(std::cout, "creatorCircle2D");
    clout << "WARNING: A rotation matrix is provided but is not necessary for a circle." << std::endl;
  }

  /// Set resolved circle 2D with given suface
  setResolvedObject<T,PARTICLETYPE,true>(particleSystem, idxParticle, idxSurface,
                                         position, density, 0., velocity);
}

/// Add resolved circle as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCircle2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          Vector<T,2> position, T radius, T epsilon, T density=0.,
                          Vector<T,2> velocity = Vector<T,2>(0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved circle 2D at given index
  setResolvedCircle2D( particleSystem, idxParticle, position, radius, epsilon, density, velocity );
}


//// CUBOID 2D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCuboid2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          std::size_t idxParticle,
                          Vector<T,2> position, Vector<T,2> extend, T epsilon, T density=0.,
                          T angle = 0.,
                          Vector<T,2> velocity = Vector<T,2> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF2D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCuboid2D<T, T, true>(Vector<T,2>(0.), extend[0], extend[1], epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  /// Set resolved cuboid 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angle, velocity);
}


/// Add resolved cuboid as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCuboid2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          Vector<T,2> position, Vector<T,2> extend, T epsilon, T density=0.,
                          T angle = 0.,
                          Vector<T,2> velocity = Vector<T,2> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cuboid 2D at given index
  setResolvedCuboid2D( particleSystem, idxParticle, position,
                       extend, epsilon, density, angle, velocity );
}


//// TRIANGLE 2D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedTriangle2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                            std::size_t idxParticle,
                            Vector<T,2> position, T radius, T epsilon, T density=0.,
                            T angle = 0.,
                            Vector<T,2> velocity = Vector<T,2> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF2D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorTriangle2D<T, T, true>(Vector<T,2>(0.), radius, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  /// Set resolved triangle 2D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angle, velocity);
}


/// Add resolved triangle as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedTriangle2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                            Vector<T,2> position, T radius, T epsilon, T density=0.,
                            T angle = 0.,
                            Vector<T,2> velocity = Vector<T,2> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved triangle 2D at given index
  setResolvedTriangle2D( particleSystem, idxParticle, position,
                         radius, epsilon, density, angle, velocity );
}


//// ARBITRARY SHAPE 2D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedArbitraryShape2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                  std::size_t idxParticle,
                                  Vector<T,2> position,
                                  T latticeSpacing,
                                  std::shared_ptr<IndicatorF2D<T>> indPtr,
                                  T epsilon, T density=0., T angle = 0.,
                                  Vector<T,2> velocity = Vector<T,2> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF2D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCustom2D<T, T, true>(latticeSpacing, indPtr, Vector<T,2> (0.), epsilon, 0.));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  /// Set resolved arbitrary shape with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angle, velocity);
}

/// Add resolved arbitrary shape as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedArbitraryShape2D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                  Vector<T,2> position,
                                  T latticeSpacing,
                                  std::shared_ptr<IndicatorF2D<T>> indPtr,
                                  T epsilon, T density=0., T angle = 0.,
                                  Vector<T,2> velocity = Vector<T,2> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved arbitrary shape 2D at given index
  setResolvedArbitraryShape2D( particleSystem, idxParticle, position, latticeSpacing, indPtr,
      epsilon, density, angle, velocity );
}

} //namespace creators

} //namespace particles

} //namespace olb

#endif
