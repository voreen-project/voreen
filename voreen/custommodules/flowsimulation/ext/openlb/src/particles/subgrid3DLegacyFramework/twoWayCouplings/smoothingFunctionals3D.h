/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Smoothings functionals for Lagrangian two-way coupling methods -- header file.
 */

#ifndef LB_SMOOTHING_FUNCTIONALS_3D_H
#define LB_SMOOTHING_FUNCTIONALS_3D_H

#include "functors/lattice/reductionF3D.h"

namespace olb {

/** Data structure for smoothing functionals.
  * Stores the lattice position of a cell within smoothing kernel length
  * and the related multiplicative weight and continuous phase volume fraction.
  */
template<typename T>
struct LatticePosAndWeight {
  // Global IC
  int globic = 0;
  // Lattice position
  int latticePos[3] = {0, 0, 0};
  // Multiplicative weight for kernel smoothing
  T weight = T();
  // Continuous phase volume fraction
  T continuousPhaseFraction = 1.;
};

/** Abstact class for all the smoothing functionals.
  */
template<typename T, typename Lattice>
class SmoothingFunctional {
public:
  // Rebuilds _latticePosAndWeight with the new cells within _kernelLength from the bubble's position
  virtual bool update(T physPosP[], int globic);
  // Returns read-only access to internal data
  const std::deque<LatticePosAndWeight<T>> getData() const
  {
    return _latticePosAndWeight;
  }
  int getNvoxelInterpPoints()
  {
    return _nVoxelInterpPoints;
  }
protected:
  /// Constructor
  SmoothingFunctional (
    T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, int nVoxelInterpPoints=2 );
  /// The actual smoothing function
  virtual T smoothingFunction(T delta)=0;
  /// Returns the weight for smoothing.
  virtual T compute(T physPosP[], T physPosL[])=0;
  T _kernelLength; // Kernel's smoothing length.
  UnitConverter<T, Lattice>& _converter;  // reference to a UnitConverter
  SuperLattice<T, Lattice>& _sLattice; // reference to a lattice
  // positions and weights of the cells within _kernelLength from bubble's position
  std::deque<LatticePosAndWeight<T> > _latticePosAndWeight;
  int _nVoxelInterpPoints;
};

/** Abstact class for all the linear-averaging smoothing functionals.
  */
template<typename T, typename Lattice>
class LinearAveragingSmoothingFunctional : public SmoothingFunctional<T, Lattice> {
protected:
  /// Constructor
  LinearAveragingSmoothingFunctional (
    T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, int nVoxelInterpPoints=2 );
  /// Returns the weight for smoothing.
  virtual T compute(T physPosP[], T physPosL[]) override;
};

/** Abstact class for all the volume-averaging smoothing functionals.
  */
template<typename T, typename Lattice>
class VolumeAveragingSmoothingFunctional : public SmoothingFunctional<T, Lattice> {
protected:
  /// Constructor
  VolumeAveragingSmoothingFunctional (
    T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, int nVoxelInterpPoints=2 );
  /// Returns the weight for smoothing.
  virtual T compute(T physPosP[], T physPosL[]) override;
};

/** Smoothing functional as in Deen et al (2004), Chem. Eng. Sci 59.
  */
template<typename T, typename Lattice>
class DeenSmoothingFunctional : public LinearAveragingSmoothingFunctional<T, Lattice> {
public:
  /// Constructor
  DeenSmoothingFunctional (
    T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
protected:
  /// The actual smoothing function
  virtual T smoothingFunction(T delta) override;
};

/** Smoothing functional as in Evrard, Denner and van Wachem (2019).
  */
template<typename T, typename Lattice>
class vanWachemSmoothingFunctional : public LinearAveragingSmoothingFunctional<T, Lattice> {
public:
  /// Constructor
  vanWachemSmoothingFunctional (
    T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, T radius, int nVoxelInterpPoints );
  // Rebuilds _latticePosAndWeight with the new cells within _kernelLength from the bubble's position
  virtual bool update(T physPosP[], int globic) override;
protected:
  /// The actual smoothing function
  virtual T smoothingFunction(T delta) override;
  /// Updates _latticePosAndWeight with contribution from continuous phase fraction. To be called AFTER compute(...) where required.
  void updateContinuousPhaseFraction(T physPosP[], int globic);
  T _radius; // particle radius (I assume all the particles in the system have the same radius!)
};

/** Stepwise smoothing functional.
  */
template<typename T, typename Lattice>
class StepSmoothingFunctional : public VolumeAveragingSmoothingFunctional<T, Lattice> {
public:
  /// Constructor
  StepSmoothingFunctional (
    T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
protected:
  /// The actual smoothing function
  virtual T smoothingFunction(T delta) override;
};

}

#endif
