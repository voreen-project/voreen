/*  DESCRIPTOR Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2016 Thomas Henn, Fabian Klemens, Robin Trunk, Davide Dapelo
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

#ifndef PARTICLEDYNAMICS_2D_H
#define PARTICLEDYNAMICS_2D_H

#include "core/superLattice2D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF2D.h"


namespace olb {

template<typename T, typename DESCRIPTOR>
class ParticleDynamics2D {
private:
  SuperLattice2D<T, DESCRIPTOR>& _sLattice;
  UnitConverter<T,DESCRIPTOR> const& _converter;
  SuperGeometry2D<T>& _superGeometry;
  std::shared_ptr<IndicatorF2D<T> > _indicatorF;
  std::vector<SmoothIndicatorF2D<T,T,true>* > _vectorOfIndicator;
  T _lengthX;
  T _lengthY;
  Vector<T,2> _accExt;
  bool _escapeFromDomain=false;
  bool _oldWallCollision=true;
  Vector<bool,2> _periodicity;
  bool _useZetaField=false;

public:
  ParticleDynamics2D(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry2D<T>& superGeometry,
                     T lengthX, T lengthY, Vector<T,2> accExt = Vector<T,2> (0.,0.),
                     Vector<bool,2> periodicity = Vector<bool,2> (false,false) )
    : _sLattice(sLattice),
      _converter(converter),
      _superGeometry(superGeometry),
      _lengthX(lengthX),
      _lengthY(lengthY),
      _accExt(accExt),
      _periodicity(periodicity)
  {}

  ParticleDynamics2D(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry2D<T>& superGeometry, std::shared_ptr<IndicatorF2D<T> > indicatorF,
                     bool escapeFromDomain, bool oldWallCollision=false,
                     Vector<T,2> accExt = Vector<T,2> (0.,0.),
                     Vector<bool,2> periodicity = Vector<bool,2> (false,false) )
    : _sLattice(sLattice),
      _converter(converter),
      _superGeometry(superGeometry),
      _indicatorF(indicatorF),
      _accExt(accExt),
      _escapeFromDomain(escapeFromDomain),
      _oldWallCollision(oldWallCollision),
      _periodicity(periodicity)
  {}

  void addCircle(Vector<T,2> center, T radius, T density, T epsilon,
                 Vector<BaseType<T>,2> vel = Vector<BaseType<T>,2> (0.,0.));

  void addCuboid(Vector<T,2> center, T xLength, T yLength, T density,
                 T epsilon, T theta=0, Vector<BaseType<T>,2> vel = Vector<BaseType<T>,2> (0.,0.));

  void addTriangle(Vector< T, 2> center, T radius, T density, T epsilon,
                   T theta, Vector<BaseType<T>,2> vel = Vector<BaseType<T>,2> (0.,0.));

  void addParticle(SmoothIndicatorF2D<T, T, true>& indicator);

  void computeBoundaryForce(std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator);

  void addWallColl(SmoothIndicatorF2D<T,T,true>& indicator, T delta);

  void verletIntegration(SmoothIndicatorF2D<T,T,true>& indicator);

  void updateParticleDynamics(std::string name, SmoothIndicatorF2D<T,T,true>& indicator);

  void checkAndRemoveEscaped();

  void addParticleField(SmoothIndicatorF2D<T,T,true>& indicator);

  void simulateTimestep(std::string name);

  // currently only for single particle / uniform movement
  void simulatePrescribedVelocity(std::string name, const Vector<T,2>& vel, const T angularVel);

  void calculateDragLiftCoefficients(SmoothIndicatorF2D<T, T, true>& indicator, Vector<T,2>& coeff, const Vector<T,2>& fluidVel = Vector<T,2> (0.,0.) );

  void print();

  void load(std::string filename, T epsilon);

  void save(std::string filename);

//  void lennardJonesColl(SmoothIndicatorF2D<T, T>& indicator, const std::vector<T>& pos2, T delta = 0.)

//  void addCollisionModel(SmoothIndicatorF2D<T,T,true>& indicator, SmoothIndicatorF2D<T,T,true>& indicator2);

//  void eulerIntegration(SmoothIndicatorF2D<T,T,true>& indicator);

//  void addParticleColl(SmoothIndicatorF2D<T,T,true>& indicator, SmoothIndicatorF2D<T,T,true>& indicator2, T delta);

  void setZetaFieldUsage(bool useZetaField);

};

}

#endif
