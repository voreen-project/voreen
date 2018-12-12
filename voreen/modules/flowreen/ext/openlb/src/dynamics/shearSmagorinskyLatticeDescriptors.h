/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathen
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

#ifndef SHEAR_SMAGORINSKY_LATTICE_DESCRIPTOR_H
#define SHEAR_SMAGORINSKY_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {

// 2D Descriptors for flow with Shear-Improved Smagorinsky

struct ShearSmagorinsky2dDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
};

struct ShearSmagorinsky2dDescriptorBase {
  typedef ShearSmagorinsky2dDescriptor ExternalField;
};

template <typename T> struct ShearSmagorinskyD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ShearSmagorinsky2dDescriptorBase {
};


/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Smagorinsky

struct ShearSmagorinsky3dDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
};

struct ShearSmagorinsky3dDescriptorBase {
  typedef ShearSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct ShearSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ShearSmagorinsky3dDescriptorBase {
};


/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Forced Shear-Improved Smagorinsky

struct ShearSmagorinskyForced3dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
  static const int forceBeginsAt    = 1;
  static const int sizeOfForce      = 3;
};

struct ShearSmagorinskyForced3dDescriptorBase {
  typedef ShearSmagorinskyForced3dDescriptor ExternalField;
};

template <typename T> struct ShearSmagorinskyForcedD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ShearSmagorinskyForced3dDescriptorBase {
};

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Forced Shear-Improved Smagorinsky

struct ForcedShearWallSmagorinsky3dDescriptor {
  static const int numScalars = 5;
  static const int numSpecies = 3;
  static const int avShearIsAt = 0;
  static const int sizeOfAvShear = 1;
  static const int forceBeginsAt    = 1;
  static const int sizeOfForce      = 3;
  static const int tauWIsAt   = 4;
  static const int sizeOfTauW      = 1;
};

struct ForcedShearWallSmagorinsky3dDescriptorBase {
  typedef ForcedShearWallSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct ForcedShearWallSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ForcedShearWallSmagorinsky3dDescriptorBase {
};

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Kalman Finitie Difference Smagorinsky

struct FDKalmanShearSmagorinsky3dDescriptor {
  static const int numScalars = 23;
  static const int numSpecies = 4;
  static const int ErrorCovarianceIsAt = 0;
  static const int sizeOfErrorCovariance = 1;
  static const int VarianceIsAt = 1;
  static const int sizeOfVariance = 1;
  static const int velocityIsAt = 2;
  static const int sizeOfVelocity = 3;
  static const int FilteredvelGradIsAt = 5;
  static const int sizeOfFilteredVelGrad = 9;
  static const int velGradIsAt = 14;
  static const int sizeOfVelGrad = 9;
};

struct FDKalmanShearSmagorinsky3dDescriptorBase {
  typedef FDKalmanShearSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct FDKalmanShearSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public FDKalmanShearSmagorinsky3dDescriptorBase {
};

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Kalman Finitie Difference Smagorinsky

struct FDKalmanShearSmagorinskyForced3dDescriptor {
  static const int numScalars = 26;
  static const int numSpecies = 5;
  static const int ErrorCovarianceIsAt = 0;
  static const int sizeOfErrorCovariance = 1;
  static const int VarianceIsAt = 1;
  static const int sizeOfVariance = 1;
  static const int velocityIsAt = 2;
  static const int sizeOfVelocity = 3;
  static const int FilteredvelGradIsAt = 5;
  static const int sizeOfFilteredVelGrad = 9;
  static const int velGradIsAt = 14;
  static const int sizeOfVelGrad = 9;
  static const int forceBeginsAt    = 23;
  static const int sizeOfForce      = 3;
};

struct FDKalmanShearSmagorinskyForced3dDescriptorBase {
  typedef FDKalmanShearSmagorinskyForced3dDescriptor ExternalField;
};

template <typename T> struct FDKalmanShearSmagorinskyForcedD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public FDKalmanShearSmagorinskyForced3dDescriptorBase {
};

////////////////////////////////////////////////////
// Kalman filter : Adaptive exponential smoothing //
////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Smagorinsky - Kalman Filter
// Boudet et al. (2016) A Kalman filter adapted of the estimation of mean gradients
//   in the a large-eddy simulation of unsteady turbulent flows.

struct KalmanShearSmagorinsky3dDescriptor {
  static const int numScalars = 22;
  static const int numSpecies = 4;
  static const int ErrorCovarianceIsAt = 1;
  static const int sizeOfErrorCovariance = 1;
  static const int VarianceIsAt = 2;
  static const int sizeOfVariance = 1;
  static const int TauSgsIsAt = 3;
  static const int sizeOfTauSgs = 1;
  static const int FilteredPopulationIsAt = 4;
  static const int sizeOfFilteredPopulation = 19;
};

struct KalmanShearSmagorinsky3dDescriptorBase {
  typedef KalmanShearSmagorinsky3dDescriptor ExternalField;
};

template <typename T> struct KalmanShearSmagorinskyD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public KalmanShearSmagorinsky3dDescriptorBase {
};


} // namespace descriptors

} // namespace olb

#endif
