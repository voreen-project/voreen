/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Simon Berg
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

#include "ParticleDistribution.h"
#include <utility>
#ifndef PARTICLE_DISTRIBUTION_H
#define PARTICLE_DISTRIBUTION_H

namespace olb {
  template <typename T, typename S>
  LogNormalDistribution<T, S>::LogNormalDistribution(T standardDeviation, T geometricMean, T scale) :AnalyticalF1D<T, S>(1), _standardDeviation(standardDeviation), _geometricMean(geometricMean), _scale(scale) {
    this->getName() = "LogNDistribution";
  };

  template <typename T, typename S>
  bool LogNormalDistribution<T, S>:: operator()(T output[1], const S x[1]) {
    output[0] = _scale / (x[0] * _standardDeviation * util::pow((2 * M_PI), 0.5)) * util::exp(-util::pow(util::log(x[0] - _geometricMean), 2.) / (2 * util::pow((_standardDeviation), 2.)));

    return true;
  };

  template <typename T, typename S>
  GaussDistribution<T, S>::GaussDistribution(T standardDeviation, T mean, T scale) :AnalyticalF1D<T, S>(1), _standardDeviation(standardDeviation), _mean(mean), _scale(scale) {
    this->getName() = "GaussDistribution";
  };

  template <typename T, typename S>
  bool GaussDistribution<T, S>::operator()(T output[1], const S x[1]) {
    output[0] = _scale / (_standardDeviation * util::pow((2 * M_PI), 0.5)) * util::exp(-0.5 * util::pow((x[0] - _mean) / (_standardDeviation), 2.));
    return true;
  };

  /* When creating a ParticleDistribution a matrix ('_sizeCalculation') is created with evenly spaced  
  * particle categories ranging from 'minRadius' to 'maxRadius' and each concentration is calculated  
  * with the given distribution
  */
  template <typename T, typename S>
  template <typename C>
  ParticleDistribution<T, S>::ParticleDistribution(FunctorPtr<AnalyticalF1D<T, S>>&& distribution, C numberOfParticles, C numberOfParticleSizeGroups, T minRadius, T maxRadius, T lengthConversion) :
    _distribution(std::move(distribution)), _numberOfParticles(numberOfParticles), _numberOfParticleSizeGroups(numberOfParticleSizeGroups), _minRadius(minRadius), _maxRadius(maxRadius) {
    _sizeCalculation = new T * [_numberOfParticleSizeGroups];
    const T delta_d = (_maxRadius - _minRadius) / _numberOfParticleSizeGroups;
    for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
      _sizeCalculation[i] = new T[4];
    }
    T output[1] = {};
    T input[1] = {};
    for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
      // calculate size groups
      _sizeCalculation[i][0] = (_minRadius + delta_d * i) * lengthConversion;
      input[0] = _sizeCalculation[i][0]  / lengthConversion;
      // calculate concentration
      _distribution->operator()(output, input);
      _sizeCalculation[i][1] = output[0];
    }
    calculateParticleArray();
  };
  /* For the calculation of the particleArray, the total concentration is determined
  * and for every size group the fraction (relative to all particles) and the number 
  * of particles is calculated.
  * If the exact particle number has to be met, the function will be called again
  * with the difference of particles in the particle array to the specified total  
  * number as an offset.
  */
  template <typename T, typename S>
  bool ParticleDistribution<T, S>::calculateParticleArray(bool getExactParticleNumber) {
    T sum_np = 0;
    for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
      // calculate total particle conc.
      sum_np += _sizeCalculation[i][1];
    }

    for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
      // calculate fraction
      _sizeCalculation[i][2] = _sizeCalculation[i][1] / sum_np;
      // calculate number of whole particles
      // ofset is added to correct difference in _numberOfParticles to _particleArray..size() from rounding errors
      _sizeCalculation[i][3] = round(_sizeCalculation[i][2] * (_numberOfParticles + _particleOffset));
    }

    long long p_count = 0;
    if (getExactParticleNumber) {
      for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
        p_count += _sizeCalculation[i][3];
      }
    }
    if (p_count == _numberOfParticles || getExactParticleNumber == false) {
      p_count = 0;
      for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
        if (_sizeCalculation[i][3] > 0) {
          for (long long j = 1; j <= _sizeCalculation[i][3]; j++) {
            // write particle radii to array, where every entry represents 1 particle
            _particleArray.push_back(_sizeCalculation[i][0]);
            p_count++;
          }
        }
      }
      return true;
    }
    else {
      _particleArray.clear();
      // calculate offset and calls function again with offset for convergence, choose 
      // different (larger) divisor if function does not converge
      _particleOffset += (_numberOfParticles - p_count)/7.;

      return this->calculateParticleArray(true);
    }
  };
  template <typename T, typename S>
  void ParticleDistribution<T, S>::getParticleArray(std::vector<T>& particleArray) {
    particleArray = _particleArray;
  };
  /* returns the next particle radius in particleArray, loops back to the first entry and/or 
  * shuffles the array before returning first particle if specified.
  */
  template <typename T, typename S>
  T ParticleDistribution<T, S>::nextParticleRadius(bool loopAround, bool shuffleParticles) {
    if (_currentParticle == 0 && shuffleParticles) {
      shuffleParticleArray();
      T Radius = _particleArray[_currentParticle];
      _currentParticle++;
      return Radius;
    }
    else if (_currentParticle == (_numberOfParticles)) {
      if (loopAround) {
        T Radius = _particleArray[_currentParticle];
        _currentParticle = 0;
        if (shuffleParticles) {
          shuffleParticleArray();
        }
        return Radius;
      }
      else {
        return 0;
      }
    }
    else {
      T Radius = _particleArray[_currentParticle];
      _currentParticle++;
      return Radius;
    }
  };
  /// shuffles particle array using random_shuffle
  template <typename T, typename S>
  void ParticleDistribution<T, S>::shuffleParticleArray() {
    random_shuffle(&_particleArray[0], &_particleArray[_particleArray.size() - 1]);
  };

  template <typename T, typename S>
  /// prints the calculation matrix to the console
  void ParticleDistribution<T, S>::printSizeMatrix() {
    OstreamManager clout(std::cout, "printSizeMatrix");
    clout << std::setw(7) << "GroupNo" << std::setw(12) << "Radius" << std::setw(12) << "concentr." << std::setw(12) << "partiton" << std::setw(12) << "p_number" << std::endl;
    for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
      clout << std::setw(7) << i + 1;
      for (long long j = 0; j <= 3; j++) {
        clout << std::setw(12) << _sizeCalculation[i][j];
      }
      clout << std::endl;
    }
  };

  template <typename T, typename S>
  /// prints the particle Radii to the console
  void ParticleDistribution<T, S>::printParticleArray() {
    OstreamManager clout(std::cout, "printParticleArray");
    clout << std::setw(7) << "p_number" << std::setw(12) << "radius" << std::endl;
    for (long long i = 0; i < _particleArray.size(); i++) {
      clout << std::setw(7) << i << std::setw(12) << _particleArray[i] << std::endl;
    }
  };
  /* To spread out the particle creation according to a timeDistribution a time array
  * is calculated by computing the time discrete integral in the interval ['begin', 'end']
  * and the timeDistribution is evaluated for time-step and added up until a complete particle
  * is reached. Every time-step a complete particle should be spawned is stored in '_timeArray'.
  */
  template <typename T, typename S>
  template <typename DESCRIPTOR>
  bool ParticleDistribution<T, S>::calculateTimeArray(AnalyticalF1D<T, S>& timeDistribution, UnitConverter<T, DESCRIPTOR>& converter, T begin, T end) {
    _timeArray.clear();
    T input[1] = {  };
    T nSum = 0;
    T output[1] = { };
    T integral = 0;
    // calculate discrete integral 
    for (std::size_t iT = converter.getLatticeTime(begin); iT <= converter.getLatticeTime(end); iT++) {
      input[0] = converter.getPhysTime(iT);
      timeDistribution(output, input);

      integral += output[0];
    }
    // evaluate timeDistribution for every time-step insde the time intervall
    for (std::size_t iT = converter.getLatticeTime(begin); iT <= converter.getLatticeTime(end); iT++) {
      input[0] = converter.getPhysTime(iT);
      timeDistribution(output, input);
      output[0] *= ((_numberOfParticles + _particleOffsetTime) / integral);
      if (output[0] > 0) {
        // adding up particle parts until a complete particle can be added
        nSum += output[0];
        if (nSum >= 1.) {
          for (int i = 1; i <= nSum; ++i) {
            // add the value of the time-step for one particle to _timeArray
            _timeArray.push_back(iT);
          }
          nSum = 0;
        }
      }
    }
    // check if all calculated particles match total number of particles
    if ((long long)(_timeArray.size()) == _numberOfParticles) {
      return true;
    }
    else {
      // calculate offset and calls function again with offset for convergence, choose 
      // different (larger) divisor if function does not converge
      _particleOffsetTime += (_numberOfParticles - _timeArray.size())/7.;
      calculateTimeArray(timeDistribution, converter, begin, end);
      return false;
    }
  };

  /// prints out all time-steps a new particle should appear
  template <typename T, typename S>
  void ParticleDistribution<T, S>::printTimeArray() {
    OstreamManager clout(std::cout, "printTimeArray");
    clout << std::setw(7) << "p_number" << std::setw(12) << "time" << std::endl;
    for (long long i = 0; i < _timeArray.size(); i++) {
      clout << std::setw(7) << i << std::setw(12) << _timeArray[i] << std::endl;
    }
  };
  
  /* this function can be called every time-step to add a new spherical particle when
  *  one should be added to the material number specified by the indicator. 
  * 'loopAround' only affects particle array after reaching last particle
  * 'shuffleParticles=true' shuffles particle array before adding the first particle and after last 
  */
  template <typename T, typename S>
  template <typename C>
  void ParticleDistribution<T, S>::spawnSphericalParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, IndicatorF3D<S>& indicator, T density, C iT, bool loopAround, bool shuffleParticles) {
    while (_timeArray[_currentParticle] == iT) {
      T p_radius = nextParticleRadius(loopAround, shuffleParticles);
      if (p_radius > 0) {
        supParticleSystem.addParticle(indicator, 4. / 3. * M_PI * util::pow(p_radius, 3) * density, p_radius, 1);
      }
    }
  }; 

  /* This function needs to be called at every time-step iT. It checks if according to the pre-calculated time-array a 
  * particle has to be added. It then activates the next particle in line.
  */
  template <typename T, typename S>
  template <typename C>
  void ParticleDistribution<T, S>::timeActivateParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, C iT, C idOfset) {
    // check if a new particle has to be activated 
    while (_timeArray[_currentParticle] == iT) {
      // Searches all sub-particle systems
      for (unsigned int pSystems = 0; pSystems < supParticleSystem.getParticleSystems().size(); pSystems++) {
        // checks id of particle in sub-particle system
        for (int particles = 0; particles < supParticleSystem.getParticleSystems()[pSystems]->size(); particles++) {
          unsigned int currentPiD = supParticleSystem.getParticleSystems()[pSystems]->getParticlesPointer()[particles]->getID();
          // checks id of particle in sub-particle system
          if (currentPiD == _currentParticle + 1 + idOfset) {
            // activates particles
            supParticleSystem.getParticleSystems()[pSystems]->getParticlesPointer()[particles]->setActive(true);
            }
        }
      }
      _currentParticle++;
      // starts at beginning when last particle is reached 
      if(_currentParticle == _numberOfParticles){
        _currentParticle = 0;
      }
    }
  };

  /// Sets the particle state 'active' to 'false' for particles with ids from 'beginID' to 'endID'  
  template <typename T, typename S>
  void ParticleDistribution<T, S>::deactivateParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, int beginID, int endID) {
    // Searches all sub-particle systems
      for (int pSystems = 0; pSystems < supParticleSystem.numOfPSystems(); pSystems++) {
        // searches all particles in sub-particle systems
        for (int particles = 0; particles < supParticleSystem.getParticleSystems()[pSystems]->size(); particles++) {
          // checks id of particle in sub-particle system
          int currentPiD = supParticleSystem.getParticleSystems()[pSystems]->getParticlesPointer()[particles]->getID();
          if (currentPiD <= endID && currentPiD >= beginID) {
            // deactivates particles
            supParticleSystem.getParticleSystems()[pSystems]->getParticlesPointer()[particles]->setActive(false);
          }
        }
      }
  };
  /*This function adds all particles that are specified by the particle array to the ParticleSystem inside the indicator as spherical particles with
  * their density. It sets their ID from 'beginID' to 'beginID+_numberOfParticles'. Use different IDs for every particle as they get activated by their
  * ID. If shuffle=true, the particle array gets shuffled before adding.
  */
  template <typename T, typename S>
  void ParticleDistribution<T, S>::preSpawnSphericalParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, IndicatorF3D<S>& indicator, T density, int beginID, bool shuffle) {
    // call nextParticleRadius to get particle radius
    T p_radius = nextParticleRadius(false, shuffle);
    do {  
      // adds a spherical particle with specified density and '_currentParticle' as its ID
      supParticleSystem.addTracerParticle(indicator, (T)(_currentParticle), 4. / 3. * M_PI * util::pow(p_radius, 3) * density, p_radius, 1);
      p_radius = nextParticleRadius(false, shuffle);
    } while (p_radius > 0);
    _currentParticle = 0;
    deactivateParticles(supParticleSystem, beginID, beginID + _numberOfParticles);
  };
  /// destructor
  template <typename T, typename S>
  ParticleDistribution<T, S>::~ParticleDistribution() {
    for (long long i = 0; i < _numberOfParticleSizeGroups; i++) {
      delete[] _sizeCalculation[i];
    }
    delete[] _sizeCalculation;
  };
};
#endif
