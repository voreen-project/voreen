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


#ifndef PARTICLE_CREATOR_FUNCTIONS_3D_H
#define PARTICLE_CREATOR_FUNCTIONS_3D_H


#include "particles/particleSystem.h"
#include "particles/functions/particleCreatorHelperFunctions.h"


namespace olb {

namespace particles {

namespace creators {

//// SPHERE 3D

/// Set resolved sphere for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedSphere3D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,3>& position, T radius, T epsilon, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorSphere3D<T, T, true>(Vector<T,3> (0.), radius, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  // Safety mechanism for wrong particle type - It is better not to use rotation matrix!
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>() ) {
    OstreamManager clout(std::cout, "creatorSphere3D");
    clout << "WARNING: A rotation matrix is provided but is not necessary for a sphere." << std::endl;
  }

  /// Set resolved sphere 3D with given suface
  setResolvedObject<T,PARTICLETYPE,true>(particleSystem, idxParticle, idxSurface,
                                         position, density, Vector<T,3>(0.), velocity);
}

/// Add resolved sphere as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedSphere3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          const Vector<T,3>& position, T radius, T epsilon, T density=0.,
                          const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved sphere 3D at given index
  setResolvedSphere3D( particleSystem, idxParticle, position, radius, epsilon, density, velocity );
}


//// CYLINDER 3D

/// Set resolved cylinder for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCylinder3D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,3>& position, const Vector<T,3>& normal,
  T height, T radius, T epsilon, T density=0.,
  const Vector<T,3>& angle = Vector<T,3> (0.),
  const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCylinder3D<T, T, true>(Vector<T,3>(0.), normal, radius, height, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  /// Set resolved cylinder 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angle, velocity);
}

/// Add resolved cylinder as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCylinder3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                            const Vector<T,3>& position, const Vector<T,3>& normal,
                            T height, T radius, T epsilon, T density=0.,
                            const Vector<T,3>& angle = Vector<T,3> (0.),
                            const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cylinder 3D at given index
  setResolvedCylinder3D( particleSystem, idxParticle, position,
                         normal, height, radius, epsilon, density, angle, velocity );
}


//// CUBOID 3D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCuboid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          std::size_t idxParticle,
                          const Vector<T,3>& position, const Vector<T,3>& extend, T epsilon, T density=0.,
                          const Vector<T,3>& angle = Vector<T,3> (0.),
                          const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCuboid3D<T, T, true>(extend[0], extend[1], extend[2], Vector<T,3>(0.), epsilon ));

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
void addResolvedCuboid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          const Vector<T,3>& position, const Vector<T,3>& extend, T epsilon, T density=0.,
                          const Vector<T,3>& angle = Vector<T,3> (0.),
                          const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cuboid 3D at given index
  setResolvedCuboid3D( particleSystem, idxParticle, position,
                       extend, epsilon, density, angle, velocity );
}


//// CONE 3D

/// Set resolved cone for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCone3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                        std::size_t idxParticle,
                        const Vector<T,3>& position1, const Vector<T,3>& position2,
                        T radius1, T radius2, T epsilon, T density=0,
                        const Vector<T,3>& angle = Vector<T,3> (0.),
                        const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCone3D<T, T, true>(position1, position2, radius1, radius2, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );
  const Vector<T,3> position = sIndicatorSPtr->calcCenterOfMass();

  /// Set resolved cone 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angle, velocity);
}


/// Add resolved cone as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCone3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                        const Vector<T,3>& position1, const Vector<T,3>& position2,
                        T radius1, T radius2, T epsilon, T density=0,
                        const Vector<T,3>& angle = Vector<T,3> (0.),
                        const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cone 3D at given index
  setResolvedCone3D( particleSystem, idxParticle,
                     position1, position2, radius1, radius2, epsilon, density, angle, velocity );
}


//// Ellipsoid 3D

/// Set resolved ellipsoid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedEllipsoid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                             std::size_t idxParticle,
                             const Vector<T,3>& position,
                             const Vector<T,3>& radius,
                             T epsilon, T density=0.,
                             const Vector<T,3>& angle = Vector<T,3> (0.),
                             const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorEllipsoid3D<T, T, true>(Vector<T,3>(0.), radius, epsilon, Vector<T,3>(0.) ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( sIndicatorSPtr );

  /// Set resolved ellipsoid 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angle, velocity);
}


/// Add resolved ellipsoid as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedEllipsoid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                             const Vector<T,3>& position,
                             const Vector<T,3>& radius,
                             T epsilon, T density=0.,
                             const Vector<T,3>& angle = Vector<T,3> (0.),
                             const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved ellipsoid 3D at given index
  setResolvedEllipsoid3D( particleSystem, idxParticle, position,
                          radius, epsilon, density, angle, velocity );
}


//// ARBITRARY SHAPE 3D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedArbitraryShape3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                  std::size_t idxParticle,
                                  Vector<T,3> position,
                                  T latticeSpacing,
                                  std::shared_ptr<IndicatorF3D<T>> indPtr,
                                  T epsilon, T density=0.,
                                  Vector<T,3> angle = Vector<T,3> (0.),
                                  Vector<T,3> velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::shared_ptr<SIndicatorBaseType> sIndicatorSPtr(
    new SmoothIndicatorCustom3D<T, T, true>(latticeSpacing, indPtr, Vector<T,3> (0.), epsilon, Vector<T,3> (0.)));

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
void addResolvedArbitraryShape3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                  Vector<T,3> position,
                                  T latticeSpacing,
                                  std::shared_ptr<IndicatorF3D<T>> indPtr,
                                  T epsilon, T density=0.,
                                  Vector<T,3> angle = Vector<T,3> (0.),
                                  Vector<T,3> velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved arbitrary shape 3D at given index
  setResolvedArbitraryShape3D( particleSystem, idxParticle, position,
      latticeSpacing, indPtr, epsilon, density, angle, velocity );
}


} //namespace creators

} //namespace particles

} //namespace olb

#endif
