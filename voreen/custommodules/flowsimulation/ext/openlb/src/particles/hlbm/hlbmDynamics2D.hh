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



#ifndef PARTICLEDYNAMICS_2D_HH
#define PARTICLEDYNAMICS_2D_HH

#include "functors/analytical/indicator/indicCalc2D.h"
#include "functors/lattice/latticeIndicatorSmoothIndicatorIntersection2D.h"
#include "hlbmDynamics2D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::addCircle(Vector< T, 2> center, T radius, T density, T epsilon, Vector<BaseType<T>,2> vel)
{
  _vectorOfIndicator.push_back (
    new SmoothIndicatorCircle2D<T,T,true>(center, radius, epsilon, density, vel) );
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::addCuboid(Vector< T, 2> center, T xLength, T yLength, T density, T epsilon, T theta, Vector<BaseType<T>,2> vel)
{
  _vectorOfIndicator.push_back (
    new SmoothIndicatorCuboid2D<T,T,true>(center, xLength, yLength, epsilon, theta, density, vel) );
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::addTriangle(Vector< T, 2> center, T radius, T density, T epsilon, T theta, Vector<BaseType<T>,2> vel)
{
  _vectorOfIndicator.push_back (
    new SmoothIndicatorTriangle2D<T,T,true>(center, radius, epsilon, theta, density, vel) );
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::addParticle(SmoothIndicatorF2D<T,T,true>& indicator)
{
  _vectorOfIndicator.push_back(&indicator);
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::computeBoundaryForce(std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator)
{
#ifdef FEATURE_HLBM_MLA
  SuperLatticePorousMomentumLossForce2D<T, DESCRIPTOR> force(_sLattice, _superGeometry, indicator,_converter);
#else
  SuperLatticeMomentumExchangeForce2D<T, DESCRIPTOR> force(_sLattice, _superGeometry, indicator, _converter, _periodicity);
#endif

  T sumF[force.getTargetDim()];
  for (int i=0; i<force.getTargetDim(); i++) {
    sumF[i]=0.;
  }
  int input[1];
  force(sumF, input);
  for (typename std::vector<SmoothIndicatorF2D<T,T,true>* >::size_type iInd=0; iInd!=indicator.size(); iInd++) {
    /// get particle acceleration through boundary force and gravity (and buoyancy)
    Vector<T,2> acceleration2;
    Vector<T,2> force;
    force[0] = sumF[0+4*iInd];
    force[1] = sumF[1+4*iInd];
    acceleration2[0] = sumF[0+4*iInd] / indicator[iInd]->getMass() + _accExt[0];
    acceleration2[1] = sumF[1+4*iInd] / indicator[iInd]->getMass() + _accExt[1];
    indicator[iInd]->setHydrodynamicForce( force );
    indicator[iInd]->setAcc2( acceleration2 );
    indicator[iInd]->setAlpha2( sumF[2+4*iInd] / indicator[iInd]->getMofi() );
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::addWallColl(SmoothIndicatorF2D<T,T,true>& indicator, T delta)
{
  T w1 = 1e-7 / 2.;
  T w = 1e-7 / 2.;

  T rad = indicator.getRadius();
  T massInv = 1. / indicator.getMass();

  std::vector<T> dx(2, T());
  dx[0] = _lengthX - indicator.getPos()[0];
  dx[1] = _lengthY - indicator.getPos()[1];

  for (int i = 0; i < 2; i++) {
    if (dx[i] <= rad) {
      indicator.getAcc2()[i] += massInv * -dx[i] * (rad - dx[i]) / w1;
    }
    if (indicator.getPos()[i] <= rad) {
      indicator.getAcc2()[i] += massInv * indicator.getPos()[i] * (rad - indicator.getPos()[i]) / w1;
    }
    if (dx[i] > rad && dx[i] <= rad + delta) {
      indicator.getAcc2()[i] += massInv * -dx[i] * std::pow((rad + delta - dx[i]), 2) / w;
    }
    if (indicator.getPos()[i] > rad && indicator.getPos()[i] <= rad + delta) {
      indicator.getAcc2()[i] += massInv * indicator.getPos()[i] * std::pow((rad + delta - indicator.getPos()[i]), 2) / w;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::verletIntegration(SmoothIndicatorF2D<T,T,true>& indicator)
{
  T time = _converter.getConversionFactorTime();
  T time2 = time*time;

  Vector<T,2> position, velocity;
  Vector<T,4> rotationMatrix;
  for (int i=0; i<2; i++) {
    position[i] = indicator.getPos()[i] + indicator.getVel()[i] * time + (0.5 * indicator.getAcc()[i] * time2);
    T avgAcc = (indicator.getAcc()[i] + indicator.getAcc2()[i]) * 0.5;
    velocity[i] = indicator.getVel()[i] + avgAcc * time;
  }
  indicator.setPos(position);
  indicator.setVel(velocity);
  indicator.setAcc(indicator.getAcc2());

  indicator.setTheta( std::fmod( indicator.getTheta() + indicator.getOmega() * time + (0.5 * indicator.getAlpha() * time2), 2.*M_PI ) );
  T avgAlpha = (indicator.getAlpha() + indicator.getAlpha2()) * 0.5;
  indicator.setOmega( indicator.getOmega() + avgAlpha * time);
  indicator.setAlpha( indicator.getAlpha2() );

  T cos = std::cos(indicator.getTheta());
  T sin = std::sin(indicator.getTheta());

  rotationMatrix[0] = cos;
  rotationMatrix[1] = sin;
  rotationMatrix[2] = -sin;
  rotationMatrix[3] = cos;
  indicator.setRotationMatrix(rotationMatrix);
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::updateParticleDynamics(std::string name, SmoothIndicatorF2D<T,T,true>& indicator)
{
//    if (name == "euler")
//      this->eulerIntegration(indicator);
//    else if (name == "verlet")
  this->verletIntegration(indicator);
//    else
//      std::cout << "ERROR: no valid integration...use 'euler' or 'verlet'"
//                << std::endl;
}

// TODO: Needs work as one evaluation of "layer" requires 3-5 indicator
// evaluations (each a virtual function calls), this may be costly
template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::checkAndRemoveEscaped()
{
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    auto internal = std::make_shared<IndicatorLayer2D<T> > (*_indicatorF.get(), -_converter.getConversionFactorLength()+(**i).getEpsilon() );
    IndicMinus2D<T> layer(_indicatorF, internal);
    SuperLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,true> intersection(
      _sLattice, _superGeometry, layer, **i );
    int input[1];
    T output[1] = {0.};
    intersection(output, input);
    if (output[0] == 1) {
      _vectorOfIndicator.erase(i);
      if (i==_vectorOfIndicator.end() ) {
        break;
      }
    }
  }
}


template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::addParticleField(SmoothIndicatorF2D<T,T,true>& indicator)
{
  /// Analytical2D functor for particle motion (trans+rot)
  ParticleU2D<T,T,DESCRIPTOR> velocity(indicator, _converter);
  if (_useZetaField) {
    setSuperZetaParticleField(_superGeometry, velocity, indicator, _sLattice);
  } else if (_periodicity[0] || _periodicity[1]) {
    setSuperExternalParticleField(_superGeometry, velocity, indicator, _sLattice, _periodicity);
  } else {
    setSuperExternalParticleField(_superGeometry, velocity, indicator, _sLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::simulateTimestep(std::string name)
{
  // Remove particles from domain if _escapeFromDomain is toggled
  if (_escapeFromDomain) {
    checkAndRemoveEscaped();
  }

  // Compute force acting on particles boundary
  computeBoundaryForce(_vectorOfIndicator);

  // Update particle dynamics and porous particle field
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    updateParticleDynamics(name, **i);
    addParticleField(**i);
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::simulatePrescribedVelocity(std::string name, const Vector<T,2>& vel, const T angularVel)
{
  computeBoundaryForce(_vectorOfIndicator);

  Vector<T,2> zero = {0.,0.};
  // Update particle dynamics and porous particle field
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    updateParticleDynamics(name, **i);
    (**i).setOmega(angularVel);
    (**i).setVel(vel);
    (**i).setAcc(zero);
    (**i).setAcc2(zero);
    (**i).setAlpha(0.);
    (**i).setAlpha2(0.);
    addParticleField(**i);
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::calculateDragLiftCoefficients(SmoothIndicatorF2D<T, T, true>& indicator, Vector<T,2>& coeff, const Vector<T,2>& fluidVel)
{
  // calculate area of indicator object in direction of main axes
  SuperGeometryFacesIndicator2D<T,true> f(_superGeometry, indicator, 1, _converter.getPhysDeltaX());
  T area[5];
  int input[3];
  f(area, input);

  // calculate relative velocity
  for (int i=0; i<2; i++) {
    T relativeVelSqr = std::pow( indicator.getVel()[i]-fluidVel[i], 2.);
    coeff[i] = 2.*indicator.getHydrodynamicForce()[i] / (_converter.getPhysDensity()*area[i]*relativeVelSqr);
  }

  std::cout << "////// " << indicator.getHydrodynamicForce()[0] << " -- " << area[0] << " -- " << std::pow( indicator.getVel()[0]-fluidVel[0], 2.) << std::endl;

}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::print()
{
  OstreamManager clout(std::cout, "ParticleDynamics2D");
  clout << "Number of particles=" << _vectorOfIndicator.size() << std::endl;
  int count = 0;
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    clout << "Particle " << ++count << " (" << (**i).name() << "):" << std::endl;
    clout << " |Circum radius(m)=     " << std::setw(13) << (**i).getCircumRadius() <<std::endl;
    clout << " |Mass(kg)=             " << std::setw(13) << (**i).getMass() << std::endl;
    clout << " |Position(m)=         (" << std::setw(13) << (**i).getPos()[0] << ", " << std::setw(13) << (**i).getPos()[1] << ")" << std::endl;
    clout << " |Angle(°)=             " << std::setw(13) << (**i).getTheta()*(180/M_PI) << std::endl;
    clout << " |Velocity(m/s)=       (" << std::setw(13) <<(**i).getVel()[0] << ", " << std::setw(13) << (**i).getVel()[1] << ")"  << std::endl;
    clout << " |Ang. velocity(°/s)=   " << std::setw(13) << (**i).getOmega()*(180/M_PI) << std::endl;
    clout << " |Hydro. Force(N)=     (" << std::setw(13) <<(**i).getHydrodynamicForce()[0] << ", " << std::setw(13) << (**i).getHydrodynamicForce()[1] << ")"  << std::endl;
    clout << " |Acceleration(m/s^2)= (" << std::setw(13) << (**i).getAcc()[0] << ", " << std::setw(13) << (**i).getAcc()[1] << ")"  << std::endl;
    clout << " |Ang. acc.(°/s^2)=     " << std::setw(13) << (**i).getAlpha()*(180/M_PI) << std::endl;
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::load(std::string filename, T epsilon)
{
  std::ifstream iout( (singleton::directories().getLogOutDir() + filename).c_str() );

  while (iout) {
    std::string name = "";
    iout >> name;

    if (name == "circle") {
      T radius = T(), mass = T();
      Vector<T, 2> center = {T(), T()};
      Vector<T, 2> vel = {T(), T()};
      iout >> center[0] >> center[1] >> radius >> mass >> vel[0] >> vel[1];
      addCircle(center, radius, mass, epsilon, vel);
    }
  }

  iout.close();
}

// WORKS ONLY FOR CIRCLES FOR NOW!!
template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::save(std::string filename)
{
  std::ofstream fout( (singleton::directories().getLogOutDir() + filename).c_str() );

  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    if ( (**i).name() == "circle") {
      fout << (**i).name()            << " "
           << (**i).getPos()[0]       << " " << (**i).getPos()[1] << " "
           << (**i).getCircumRadius() << " " << (**i).getMass()   << " "
           << (**i).getVel()[0]       << " " << (**i).getVel()[1]
           << std::endl;
    }
  }

  fout.close();
}

//template<typename T, typename DESCRIPTOR>
//void ParticleDynamics2D<T, DESCRIPTOR>::addCollisionModel(SmoothIndicatorF2D<T,T,true>& indicator,
//    SmoothIndicatorF2D<T,T,true>& indicator2) {
//this->addParticleColl(indicator, indicator2, _converter.getLatticeL());
//this->addWallColl(indicator, 10. * _converter.getLatticeL());
//}

//template<typename T, typename DESCRIPTOR>
//void ParticleDynamics2D<T, DESCRIPTOR>::eulerIntegration(SmoothIndicatorF2D<T,T,true>& indicator) {
//    T time = _converter.physTime();
//
//    for (int i = 0; i < 2; i++) {
//      _vel[i] += _A[i] * time;
//      _pos[i] += _vel[i] * time;
//
//    }
//    _omega += _A[2] * time;   //angular velocity
//    _theta += _omega * time;  //angle
//  }

//template<typename T, typename DESCRIPTOR>
//void ParticleDynamics2D<T, DESCRIPTOR>::addParticleColl(SmoothIndicatorF2D<T,T,true>& indicator,
//    SmoothIndicatorF2D<T,T,true>& indicator2, T delta) {
//  Vector<T, 2> pos2;
//
//  T e1 = 1e-7;
//  T e = 1e-7;
//  T massInv = 1. / indicator.getMass();
//
//  T rad = indicator.getRadius() + delta;
//  pos2 = indicator2.getPos();
//  T dist2 = std::pow(indicator.getPos()[0] - pos2[0], 2) + std::pow(indicator.getPos()[1] - pos2[1], 2);
//
//  if (dist2 <= std::pow(2. * rad, 2)) {
//    indicator.getAcc2()[0] += massInv * (indicator.getPos()[0] - pos2[0]) * (2. * rad - std::sqrt(dist2))
//        / e1;
//    indicator.getAcc2()[1] += massInv * (indicator.getPos()[1] - pos2[1]) * (2. * rad - std::sqrt(dist2))
//        / e1;
//  }
//  else if (std::pow(2 * rad, 2) < dist2
//      && dist2 <= std::pow(2 * rad + delta, 2)) {
//    indicator.getAcc2()[0] += massInv * (indicator.getPos()[0] - pos2[0])
//        * std::pow((2. * rad + delta - std::sqrt(dist2)), 2) / e;
//    indicator.getAcc2()[1] += massInv * (indicator.getPos()[1] - pos2[1])
//        * std::pow((2. * rad + delta - std::sqrt(dist2)), 2) / e;
//  }
//}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics2D<T, DESCRIPTOR>::setZetaFieldUsage(bool useZetaField)
{
  _useZetaField = useZetaField;
}




}

#endif
