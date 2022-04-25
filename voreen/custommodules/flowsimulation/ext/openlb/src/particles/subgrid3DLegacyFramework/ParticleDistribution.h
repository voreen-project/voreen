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

#ifndef PARTICLE_DISTRIBUTION_HH
#define PARTICLE_DISTRIBUTION_HH

namespace olb {
  /// Functor for a logarithmic normal distribution
  template <typename T, typename S>
  class LogNormalDistribution :public AnalyticalF1D<T, S> {
  private:
    T _standardDeviation;
    T _geometricMean;
    T _scale;

  public:

    LogNormalDistribution(T standardDeviation, T geometricMean, T scale = 1);

    bool operator()(T output[1], const S x[1]);
  };

  /// Functor for a Gaussian (normal) distribution
  template <typename T, typename S>
  class GaussDistribution :public AnalyticalF1D<T, S> {
  private:
    T _standardDeviation;
    T _mean;
    T _scale;

  public:
    GaussDistribution(T standardDeviation, T mean, T scale = 1);
    bool operator()(T output[1], const S x[1]);
  };

  /// Particle distribution for time and size discretization
  /**
  * size discretization
  * This class uses a given particle distribution ('distribution') and discretizes it in
  * 'numberOfParticleSizeGroups' Radius-size categories, ranging from 'minRadius' to 'maxRadius'.
  * The particles are stored in '_particleArray', where every entry equals the particle Radius (in meter) of one
  * particle according to the distribution.
  * To get the exact No of particles use 'calculateParticleArray(true)'.
  * WARNING: leaving getExactParticleNumber = false may result in fewer or more particles than specified.
  * BUT when putting getExactParticleNumber = true it may not converge.
  * If the particles should appear random according to the distribution use 'shuffleParticleArray'.  *
  *
  * time discretization
  * To calculate when a new particle shall be spawned, use 'calculateTimeArray' with the corresponding time
  * distribution and the time interval ranging from 'begin' to 'end'.
  * Every time-step ('iT') when a new particle should be spawned is stored in '_timeArray'.
  * 'spawnSphericalParticles' can be called every time-step to automatically spawn a spherical particle according
  * to the time and size discretization.
  * Use 'loopAround = true', to start at the beginning of '_particleArray' after reaching its end and
  * 'reshuffleAfterLoop = true' to reshuffle it at the end
  * 
  * time activation
  * When dealing with large number of particles adding them will get more time consuming because of the needed
  * communication between sub-grids of all threads. To counteract this, the particles can be added as passive 
  * particles in bulk with 'preSpawnSphericalParticles' and calling 'timeActivateParticles' will activate the
  * particles according to the time distribution.
  **/
  template <typename T, typename S>
  class ParticleDistribution {
  private:
    FunctorPtr <AnalyticalF1D<T, S>> _distribution;
    std::vector<T> _particleArray;
    std::vector<T> _timeArray;
    long long _numberOfParticles;
    long long _numberOfParticleSizeGroups;
    T _particleOffset = 0;
    T _particleOffsetTime = 0;
    long long _currentParticle = 0;
    T _minRadius;
    T _maxRadius;

    T** _sizeCalculation;   // particle size calculation matrix: Radius group | concentration | fraction | number of particles in group
  public:
    template <typename C>
    ParticleDistribution(FunctorPtr<AnalyticalF1D<T, S>>&& distribution, C numberOfParticles, C numberOfParticleSizeGroups, T minRadius, T maxRadius, T lengthConversion = 1);
    bool calculateParticleArray(bool getExactParticleNumber = false);
    template <typename DESCRIPTOR>
    bool calculateTimeArray(AnalyticalF1D<T, S>& timeDistribution, UnitConverter<T, DESCRIPTOR>& converter, T begin, T end);
    void printTimeArray();
    template <typename C>
    void spawnSphericalParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, IndicatorF3D<S>& indicator, T density, C iT, bool loopAround = false, bool reshuffleAfterLoop = false);
    template<typename C>
    void timeActivateParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, C iT, C idOfset = 0);
    void deactivateParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, int beginID, int endID);
    void preSpawnSphericalParticles(SuperParticleSystem3D<T, Particle3D>& supParticleSystem, IndicatorF3D<S>& indicator, T density, int beginID = 1, bool shuffle = false);
    void getParticleArray(std::vector<T>& particleArray);
    T nextParticleRadius(bool loopAround = false, bool reshuffleAfterLoop = false);
    void shuffleParticleArray();
    void printSizeMatrix();
    void printParticleArray();
    ~ParticleDistribution();
  };
};
#endif 