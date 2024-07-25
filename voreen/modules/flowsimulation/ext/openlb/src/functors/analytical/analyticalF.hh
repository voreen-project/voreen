/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_F_HH
#define ANALYTICAL_F_HH

#include<vector>
#include<cmath>     // for lpnorm
#include<stdlib.h>  // for random
#include <iostream>
#include<time.h>

#include "analyticalF.h"
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "utilities/calc.h"
#include "utilities/vectorHelpers.h"
#include "core/radiativeUnitConverter.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {


template <unsigned D, typename T, typename S>
AnalyticalComposed<D,T,S>::AnalyticalComposed( std::vector<AnalyticalF<D,T,S>>& f)
  : AnalyticalF<D,T,S>(f.size())
{
  if (D != f.size()) {
    throw std::length_error("D does not match with f.size().");
  }
  this->getName() = "composed";
  for (auto&& i : f) {
    _f.push_back(i);
  }
}

template <unsigned D, typename T, typename S>
bool AnalyticalComposed<D,T,S>::operator()(T output[], const S x[])
{
  for (unsigned int i=0; i<D; ++i) {
    T outputTmp[_f[i].get().getTargetDim()];
    _f[i].get()(outputTmp,x);
    output[i]=outputTmp[0];
  }
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalConst<D,T,S>::AnalyticalConst(const std::vector<T>& value)
  : AnalyticalF<D,T,S>(value.size()), _c(value)
{
  this->getName() = "const";
}

template <unsigned D, typename T, typename S>
AnalyticalConst<D,T,S>::AnalyticalConst(T value)
  : AnalyticalF<D,T,S>(1)
{
  _c.push_back(value);
  this->getName() = "const";
}

template <unsigned D, typename T, typename S>
AnalyticalConst<D,T,S>::AnalyticalConst(T value0, T value1)
  : AnalyticalF<D,T,S>(2)
{
  _c.push_back(value0);
  _c.push_back(value1);
  this->getName() = "const";
}

template <unsigned D, typename T, typename S>
AnalyticalConst<D,T,S>::AnalyticalConst(T value0, T value1, T value2)
  : AnalyticalF<D,T,S>(3)
{
  _c.push_back(value0);
  _c.push_back(value1);
  _c.push_back(value2);
  this->getName() = "const";
}

template <unsigned D, typename T, typename S>
AnalyticalConst<D,T,S>::AnalyticalConst(const Vector<T,2>& value)
  : AnalyticalF<D,T,S>(2)
{
  _c.reserve(2);
  _c.push_back(value[0]);
  _c.push_back(value[1]);
  this->getName() = "const";
}

template <unsigned D, typename T, typename S>
AnalyticalConst<D,T,S>::AnalyticalConst(const Vector<T,3>& value)
  : AnalyticalF<D,T,S>(3)
{
  _c.reserve(3);
  _c.push_back(value[0]);
  _c.push_back(value[1]);
  _c.push_back(value[2]);
  this->getName() = "const";
}

template <unsigned D, typename T, typename S>
bool AnalyticalConst<D,T,S>::operator()(T output[], const S x[])
{
  for ( unsigned i = 0; i < _c.size(); ++i) {
    output[i] = _c[i];
  }
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalNormal<D,T,S>::AnalyticalNormal(std::vector<T> mean, T stdDev)
  : AnalyticalF<D,T,S>(1), _mean(mean), _stdDev(stdDev)
{
  this->getName() = "normal";
}

template <unsigned D, typename T, typename S>
bool AnalyticalNormal<D,T,S>::operator()(T output[], const S x[])
{
  T r2 = T();
  for ( unsigned i = 0; i < D; ++i) {
    r2 += (x[i] - _mean[i]) * (x[i] - _mean[i]);
  }
  output[0] = util::exp(-r2/(2*_stdDev*_stdDev)) / (_stdDev*util::sqrt(2.*M_PI));
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalRandomBase<D,T,S>::AnalyticalRandomBase()
  : AnalyticalF<D,T,S>(1),
    gen{rd()}
{
  this->getName() = "random";
}


template <unsigned D, typename T, typename S>
AnalyticalRandomUniform<D,T,S>::AnalyticalRandomUniform(T minVal, T maxVal)
  : AnalyticalRandomBase<D,T,S>(),
    distro{minVal, maxVal}
{ }

template <unsigned D, typename T, typename S>
bool AnalyticalRandomUniform<D,T,S>::operator()(T output[], const S x[])
{
  output[0] = distro(this->gen);
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalRandomNormal<D,T,S>::AnalyticalRandomNormal(T mean, T stdDev)
  : AnalyticalRandomBase<D,T,S>(),
    distro{mean, stdDev}
{ }

template <unsigned D, typename T, typename S>
bool AnalyticalRandomNormal<D,T,S>::operator()(T output[], const S x[])
{
  output[0] = distro(this->gen);
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalRandomTruncatedNormal<D,T,S>::AnalyticalRandomTruncatedNormal(T mean, T stdDev, T n)
  : AnalyticalRandomNormal<D,T,S>(mean, stdDev),
    _min(mean - n*stdDev),
    _max(mean + n*stdDev)
{ }

template <unsigned D, typename T, typename S>
bool AnalyticalRandomTruncatedNormal<D,T,S>::operator()(T output[], const S x[])
{
  do {
    output[0] = this->distro(this->gen);
  }
  while (output[0] < _min || output[0] > _max);
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalRandomOld<D,T,S>::AnalyticalRandomOld()
  : AnalyticalF<D,T,S>(1)
{
#ifdef PARALLEL_MODE_MPI
  int  nameLen, numProcs, myID;
  char processorName[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_rank(MPI_COMM_WORLD,&myID);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Get_processor_name(processorName,&nameLen);
  srand(time(0)+myID*numProcs + nameLen);
#else
  srand(time(nullptr));
#endif
  this->getName() = "random";
}

template <unsigned D, typename T, typename S>
bool AnalyticalRandomOld<D,T,S>::operator()(T output[], const S x[])
{
  output[0]=(rand()%RAND_MAX)/(T)RAND_MAX;
  return true;
}


template <typename T, typename S, typename DESCRIPTOR>
EccentricVelocityField<T,S,DESCRIPTOR>::EccentricVelocityField(
  Vector<T,DESCRIPTOR::d> position,
  Vector<T,DESCRIPTOR::d> velocity,
  Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation> angularVelocity )
  : AnalyticalF<DESCRIPTOR::d,T,S>(DESCRIPTOR::d),
    _position(position),
    _velocity(velocity),
    _angularVelocity(angularVelocity)
{
  this->getName() = "EccentricVelocityField";
}

template <typename T, typename S, typename DESCRIPTOR>
bool EccentricVelocityField<T,S,DESCRIPTOR>::operator()(T output[], const S input[])
{
  constexpr unsigned D = DESCRIPTOR::d;

  Vector<T,D> localVel = util::calculateLocalVelocity(_position, _velocity, _angularVelocity, Vector<T,D>(input));
  for (unsigned iDim=0; iDim<D; ++iDim) {
    output[iDim] = localVel[iDim];
  }

  return true;
}

template <typename T, typename S, typename DESCRIPTOR>
EccentricLatticeVelocityField<T,S,DESCRIPTOR>::EccentricLatticeVelocityField(
  Vector<T,DESCRIPTOR::d> position,
  Vector<T,DESCRIPTOR::d> velocity,
  Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation> angularVelocity,
  UnitConverter<T,DESCRIPTOR> const& converter )
  : AnalyticalF<DESCRIPTOR::d,T,S>(DESCRIPTOR::d),
    _position(position),
    _velocity(velocity),
    _angularVelocity(angularVelocity),
    _converter(converter)
{
  this->getName() = "EccentricLatticeVelocityField";
}

template <typename T, typename S, typename DESCRIPTOR>
bool EccentricLatticeVelocityField<T,S,DESCRIPTOR>::operator()(T output[], const S input[])
{
  constexpr unsigned D = DESCRIPTOR::d;

  Vector<T,D> localVel = util::calculateLocalVelocity(_position, _velocity, _angularVelocity, Vector<T,D>(input));
  for (unsigned iDim=0; iDim<D; ++iDim) {
    output[iDim] = _converter.getLatticeVelocity( localVel[iDim] );
  }

  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalSquareWave<D,T,S>::AnalyticalSquareWave(
  T period, T amplitude, T difference)
  : AnalyticalF<D,T,S>(1), _period(period),
    _amplitude(amplitude), _difference(difference)
{ }

template <unsigned D, typename T, typename S>
bool AnalyticalSquareWave<D,T,S>::operator()(T output[], const S x[])
{
  if ( util::fmod_pos(x[0],_period) <= (_difference * _period)) {
    output[0] = _amplitude;
  }
  else {
    output[0] = - _amplitude;
  }
  return true;
}


template <unsigned D, typename T, typename S>
AnalyticalSmoothedSquareWave<D,T,S>::AnalyticalSmoothedSquareWave(
  T period, T amplitude, T difference, T epsilon)
  : AnalyticalSquareWave<D,T,S>(period, amplitude, difference),
    _epsilon(epsilon)
{
  OLB_ASSERT(epsilon < difference, "Smoothing region is to large for given phase difference");
}

template <unsigned D, typename T, typename S>
bool AnalyticalSmoothedSquareWave<D,T,S>::operator()(T output[], const S x[])
{
  if ( util::fmod_pos(x[0],this->_period) <= 0.5 * _epsilon) {
    output[0] = this->_amplitude * util::sin(M_PI * util::fmod_pos(x[0],this->_period) / _epsilon);
  }
  else if ( util::fmod_pos(x[0],this->_period) <= this->_difference * this->_period - 0.5 * _epsilon ) {
    output[0] = this->_amplitude;
  }
  else if ( util::fmod_pos(x[0],this->_period) <= this->_difference * this->_period + 0.5 * _epsilon ) {
    output[0] = - this->_amplitude
       * util::sin(M_PI * (util::fmod_pos(x[0],this->_period) - this->_difference * this->_period) / _epsilon);
  }
  else if ( util::fmod_pos(x[0],this->_period) <= this->_period - 0.5 * _epsilon ) {
    output[0] = - this->_amplitude;
  }
  else {
    output[0] = this->_amplitude
       * util::sin(M_PI * (util::fmod_pos(x[0],this->_period) - this->_period) / _epsilon);
  }
  return true;
}




////////////// OLD IMPLEMENTATION //////////////////////////

//////////////////////////////////1D////////////////////////////////////////////
template <typename T, typename S>
AnalyticalLinear1D<T,S>::AnalyticalLinear1D(T a, T b)
  : AnalyticalF1D<T,S>(1), _a(a), _b(b)
{
  this->getName() = "linear";
}

template <typename T, typename S>
AnalyticalLinear1D<T,S>::AnalyticalLinear1D(S x0, T v0, S x1, T v1)
  : AnalyticalF1D<T,S>(1)
{
  if ( util::nearZero(x1-x0) ) {
    std::cout << "Error: x1-x2=0" << std::endl;
  }
  else {
    _a = ( v1-v0 ) / ( x1-x0 );
    _b = v0 - _a*x0;
  }
  this->getName() = "linear";
}

template <typename T, typename S>
bool AnalyticalLinear1D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_a*x[0] + _b;
  return true;
}


template <typename T, typename S>
AnalyticalSquare1D<T,S>::AnalyticalSquare1D(S cp, S r, T maxi)
  : AnalyticalF1D<T,S>(1), _cp(cp), _r(r), _maxi(maxi)
{
  this->getName() = "square";
}

template <typename T, typename S>
bool AnalyticalSquare1D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_maxi*(1-(x[0]-_cp) * (x[0]-_cp) / (T)_r / (T)_r);
  return true;
}



/////////////////////////////////someOtherFunctors//////////////////////////////
template <typename T, typename S>
PolynomialStartScale<T,S>::PolynomialStartScale(S numTimeSteps, T maxValue)
  : AnalyticalF1D<T,S>(1), _numTimeSteps(numTimeSteps), _maxValue(maxValue)
{
  this->getName() = "polyStartScale";
}

template <typename T, typename S>
bool PolynomialStartScale<T,S>::operator()(T output[], const S x[])
{
  output[0]=(T) x[0] / (T) _numTimeSteps;
  output[0]=_maxValue * output[0]*output[0]*output[0] * ( 10.0 + output[0] * (6.0*output[0] - 15.0) );
  return true;
}


template <typename T, typename S>
SinusStartScale<T,S>::SinusStartScale(int numTimeSteps, T maxValue)
  : AnalyticalF1D<T,S>(1), _numTimeSteps(numTimeSteps), _maxValue(maxValue)
{
  this->getName() = "sinusStartScale";
}

template <typename T, typename S>
bool SinusStartScale<T,S>::operator()(T output[], const S x[])
{
  output[0]=(_maxValue * (util::sin(-M_PI / 2.0 + (T)x[0] / (T)_numTimeSteps * M_PI) + 1.0)) / 2.0;
  return true;
}


template <typename T, typename S>
Sinus<T,S>::Sinus(T period, T amplitude)
  : AnalyticalF1D<T,S>(1), _period(period), _amplitude(amplitude)
{
  this->getName() = "Sinus";
}

template <typename T, typename S>
bool Sinus<T,S>::operator()(T output[], const S x[])
{
  output[0] = _amplitude * util::sin((x[0] / _period) * 2 * M_PI);
  return true;
}


template <typename T, typename S>
Cosinus<T,S>::Cosinus(T period, T amplitude)
  : AnalyticalF1D<T,S>(1),
    _period(period),
    _amplitude(amplitude)
{
  this->getName() = "Cosinus";
}

template <typename T, typename S>
bool Cosinus<T,S>::operator()(T output[], const S x[])
{
  output[0] = _amplitude * util::cos((x[0] / _period) * T(2) * M_PI);
  return true;
}

template <typename T, typename S>
CosinusComposite<T,S>::CosinusComposite(T period, T amplitude, T difference)
  : AnalyticalF1D<T,S>(1),
    _period(period),
    _difference(difference),
    _amplitude(amplitude)
{
  this->getName() = "CosinusComposite";
}

template <typename T, typename S>
bool CosinusComposite<T,S>::operator()(T output[], const S x[])
{
  if ( util::fmod_pos(x[0],_period) <= (_difference * _period)) {
    output[0] = _amplitude
     * util::cos(((util::fmod_pos(x[0], _period))
                   / (_difference * _period)) * M_PI);
  }
  else {
    output[0] = _amplitude
     * util::cos(((util::fmod_pos(x[0], _period) - (_difference * _period))
                   / (_period - (_difference * _period))) * M_PI + M_PI);
  }
  return true;
}




///////////////////////////////////////2D///////////////////////////////////////
template <typename T, typename S>
AnalyticalLinear2D<T,S>::AnalyticalLinear2D(T a, T b, T c)
  : AnalyticalF2D<T,S>(1), _a(a), _b(b), _c(c)
{
  this->getName() = "linear";
}

template <typename T, typename S>
AnalyticalLinear2D<T,S>::AnalyticalLinear2D(S x0, S y0, T v0, S x1, S y1,
    T v1, S x2, S y2, T v2)
  : AnalyticalF2D<T,S>(1)
{
  this->getName() = "linear";
  T n2= (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
  if ( util::nearZero(n2) ) {
    std::cout << "Error function" << std::endl;
  }
  else {
    T n0 = (y1-y0)*(v2-v0) - (v1-v0)*(y2-y0);
    T n1 = (v1-v0)*(x2-x0) - (x1-x0)*(v2-v0);
    _a = -n0 / n2;
    _b = -n1 / n2;
    _c = (x0*n0 + y0*n1 + v0*n2) / n2;
  }
}

template <typename T, typename S>
bool AnalyticalLinear2D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_a*x[0] + _b*x[1] + _c;
  return true;
}

template <typename T, typename S>
AnalyticalParticleAdsorptionLinear2D<T,S>::AnalyticalParticleAdsorptionLinear2D(T center[], T radius, T maxValue) : AnalyticalF2D<T,S>(2), _radius(radius), _maxValue(maxValue)
{
  _center[0] = center[0];
  _center[1] = center[1];
  this->getName() = "particleAdsorptionLinear2D";
}

template <typename T, typename S>
bool AnalyticalParticleAdsorptionLinear2D<T,S>::operator()(T output[], const S input[])
{
  T dist = util::sqrt((input[0]-_center[0])*(input[0]-_center[0]) + (input[1]-_center[1])*(input[1]-_center[1]));

  if (dist > _radius) {
    output[0] = T();
    output[1] = T();
    return true;
  }
  else {
    output[0] = _maxValue*(T(1) - dist/_radius)*(_center[0]-input[0])/_radius;
    output[1] = _maxValue*(T(1) - dist/_radius)*(_center[1]-input[1])/_radius;
    return true;
  }
}


//TODO: to be removed due to eccentricVelocityField
template <typename T, typename S, typename DESCRIPTOR>
ParticleU2D<T,S,DESCRIPTOR>::ParticleU2D(SmoothIndicatorF2D<T,T,true>& indicator, UnitConverter<T,DESCRIPTOR> const& converter)
  :AnalyticalF2D<T,S>(2), _indicator(indicator), _converter(converter)
{
  this->getName() = "ParticleU";
}

template <typename T, typename S, typename DESCRIPTOR>
bool ParticleU2D<T,S,DESCRIPTOR>::operator()(T output[], const S input[])
{
  Vector<T,2> tmp = _indicator.getLocalVel( Vector<S,2>(input) );
  output[0] = _converter.getLatticeVelocity(tmp[0]);
  output[1] = _converter.getLatticeVelocity(tmp[1]);

  return true;
}


///////////////////////////////////////3D///////////////////////////////////////
template <typename T, typename S>
AnalyticalLinear3D<T,S>::AnalyticalLinear3D(T a, T b, T c, T d)
  : AnalyticalF3D<T,S>(1), _a(a), _b(b), _c(c), _d(d)
{
  this->getName() = "linear";
}

template <typename T, typename S>
AnalyticalLinear3D<T,S>::AnalyticalLinear3D(S x0, S y0, S z0, T v0, S x1,
    S y1, S z1, T v1, S x2, S y2, S z2, T v2, S x3, S y3, S z3, T v3)
  : AnalyticalF3D<T,S>(1)
{
  this->getName() = "linear";
  T n = ( (y3-y0)*(x1-x0)-(x3-x0)*(y1-y0) ) * ( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) )
        +( (z3-z0)*(x1-x0)-(x3-x0)*(z1-z0) ) * ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
  if ( util::nearZero(n) ) {
    std::cout << "Error function" << std::endl;
  }
  else {
    T w = ( (y1-y0)*(x3-x0)-(x1-x0)*(y3-y0) ) * ( (v2-v0)-(x2-x0)*(v1-v0) / (x1-x0) )
          /( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) ) + (v3-v0) - (x3-x0)*(v1-v0) / (x1-x0);
    T zx = (y1-y0)*( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) )
           -(z1-z0)*( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T rx = (v1-v0)/(x1-x0) - (y1-y0)*(v2-v0) / ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) )
           +(y1-y0)*(x2-x0)*(v1-v0) / ( (y2-y0)*(x1-x0)*(x1-x0)-(x2-x0)*(y1-y0)*(x1-x0) );
    T zy = (x1-x0)*( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) );
    T ry = ( (x1-x0)*(v2-v0)-(x2-x0)*(v1-v0) ) / ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T zz = (x1-x0)*( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T h = w/n;
    _a = rx + zx*h;
    _b = ry + zy*h;
    _c = zz*h;
    _d = v0 - x0*_a - y0*_b - z0*_c;
  }
}

template <typename T, typename S>
bool AnalyticalLinear3D<T,S>::operator()(T output[], const S x[])
{
  output[0]=_a*x[0] + _b*x[1] + _c*x[2] + _d;
  return true;
}


template <typename T, typename S>
AnalyticalScaled3D<T,S>::AnalyticalScaled3D(AnalyticalF3D<T,S>& f, T scale)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f), _scale(scale)
{
  this->getName() = "scaled";
}

template <typename T, typename S>
bool AnalyticalScaled3D<T,S>::operator()(T output[], const S x[])
{
  T outputTmp[_f.getTargetDim()];
  _f(outputTmp,x);
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = _scale*outputTmp[iDim];
  }
  return true;
}


// see Mink et al. 2016 in Sec.3.1.
template <typename T, typename S, typename DESCRIPTOR>
PLSsolution3D<T,S,DESCRIPTOR>::PLSsolution3D(RadiativeUnitConverter<T,DESCRIPTOR> const& converter)
  : AnalyticalF3D<T,S>(1),
    _physSigmaEff(util::sqrt( converter.getPhysAbsorption() / converter.getPhysDiffusion() )),
    _physDiffusionCoefficient(converter.getPhysDiffusion())
{
  this->getName() = "PLSsolution3D";
}

template <typename T, typename S, typename DESCRIPTOR>
bool PLSsolution3D<T,S,DESCRIPTOR>::operator()(T output[1], const S x[3])
{
  double r = util::sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
  output[0] = 1. / (4.0*M_PI*_physDiffusionCoefficient*r) *util::exp(-_physSigmaEff * r);
  return true;
}


template <typename T, typename S, typename DESCRIPTOR>
LightSourceCylindrical3D<T,S,DESCRIPTOR>::LightSourceCylindrical3D(RadiativeUnitConverter<T,DESCRIPTOR> const& converter, Vector<T,3> center)
  : AnalyticalF3D<T,S>(1),
    _physSigmaEff(util::sqrt( converter.getPhysAbsorption() / converter.getPhysDiffusion() )),
    _physDiffusionCoefficient(converter.getPhysDiffusion()),
    _center(center)
{
  this->getName() = "LightSourceCylindrical3D";
}

template <typename T, typename S, typename DESCRIPTOR>
bool LightSourceCylindrical3D<T,S,DESCRIPTOR>::operator()(T output[1], const S x[3])
{
  double r = util::sqrt( (x[0]-_center[0])*(x[0]-_center[0]) + (x[1]-_center[1])*(x[1]-_center[1]) );
  if ( util::nearZero(r) ) {
    std::cout << "Warning: evaluation of \"LightSourceCylindrical3D\" functor close to singularity." << std::endl;
  }
  output[0] = 1. / (4.0*M_PI*_physDiffusionCoefficient*r) *util::exp(-_physSigmaEff * r);
  return true;
}


template <typename T, typename S>
Spotlight<T,S>::Spotlight(Vector<T,3> position, Vector<T,3> direction, T falloff)
  : AnalyticalF3D<T,S>(1), _position(position), _orientation(direction[0]/norm(direction), direction[1]/norm(direction), direction[2]/norm(direction)), _falloff(falloff)
{
  this->getName() = "Spotlight";
}

template <typename T, typename S>
bool Spotlight<T,S>::operator()(T output[1], const S x[3])
{
  Vector<T,3> w(x[0] -_position[0], x[1] -_position[1], x[2] -_position[2]);
  w = normalize(w);
  T cosPhi = w*_orientation;
  output[0] = 1. / (norm(w)*norm(w)) * util::pow(util::max(0., cosPhi), _falloff);
  return true;
}


template <typename T, typename S>
GaussianHill2D<T,S>::GaussianHill2D(T sigma, Vector<T,2> x0, T c0)
  : AnalyticalF2D<T,S>(1), _sigma(sigma), _x0(x0[0],x0[1]), _c0(c0)
{
  this->getName() = "GaussianHill";
}

template <typename T, typename S>
bool GaussianHill2D<T,S>::operator()(T output[1], const S x[2])
{
  output[0] =  _c0 * util::exp(- ((x[0]-_x0[0])*(x[0]-_x0[0]) + (x[1]-_x0[1])*(x[1]-_x0[1])) / (2*_sigma*_sigma) );
  return true;
}


template <typename T, typename S>
GaussianHillTimeEvolution2D<T,S>::GaussianHillTimeEvolution2D(T sigma0, T D, T t, Vector<T,2> x0, Vector<T,2> u, T c0)
  : AnalyticalF2D<T,S>(1), _sigma02(sigma0*sigma0), _D(D), _t(t), _x0(x0[0],x0[1]), _u(u[0],u[1]), _c0(c0)
{
  this->getName() = "GaussianHillTimeEvolution";
}

template <typename T, typename S>
bool GaussianHillTimeEvolution2D<T,S>::operator()(T output[1], const S x[2])
{
  output[0] = _sigma02 / (_sigma02+2*_D*_t) * _c0 * util::exp(- ((x[0]-_x0[0]-_u[0]*_t)*(x[0]-_x0[0]-_u[0]*_t) + (x[1]-_x0[1]-_u[1]*_t)*(x[1]-_x0[1]-_u[1]*_t)) / (2*(_sigma02+2*_D*_t) ));
  return true;
}


//TODO: to be removed due to eccentricVelocityField
template <typename T, typename S, typename DESCRIPTOR>
ParticleU3D<T,S,DESCRIPTOR>::ParticleU3D(SmoothIndicatorF3D<T,T,true>& indicator, UnitConverter<T,DESCRIPTOR> const& converter)
  :AnalyticalF3D<T,S>(3), _indicator(indicator), _converter(converter)
{
  this->getName() = "ParticleU";
}

template <typename T, typename S, typename DESCRIPTOR>
bool ParticleU3D<T,S,DESCRIPTOR>::operator()(T output[], const S input[])
{
  Vector<T,3> tmp = _indicator.getLocalVel( Vector<S,3>(input) );
  output[0] = _converter.getLatticeVelocity(tmp[0]);
  output[1] = _converter.getLatticeVelocity(tmp[1]);
  output[2] = _converter.getLatticeVelocity(tmp[2]);

  return true;
}




} // end namespace olb

#endif
