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

/* Helper functionals for Lagrangian two-way coupling methods -- header file.
 */

#ifndef LB_TWO_WAY_HELPER_FUNCTIONALS_H
#define LB_TWO_WAY_HELPER_FUNCTIONALS_H

#include "functors/lattice/latticeInterpPhysVelocity3D.h"
#include "functors/lattice/reductionF3D.h"

namespace olb {

/** Abstact class for all the local forward-coupling models,
  * viz., momentum coupling from fluid to particle.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice>
class TwoWayHelperFunctional {
public:
  /// Computes the momentum transfer from fluid to particle.
  virtual bool operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic )=0;
  virtual ~TwoWayHelperFunctional();
protected:
  /// Constructor
  TwoWayHelperFunctional ( UnitConverter<T, Lattice>& converter,
                           SuperLattice<T, Lattice>& sLattice );
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
  SuperLattice<T, Lattice>& _sLattice; // reference to a lattice
  std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > _interpLatticeDensity;
  std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > _interpLatticeVelocity ;
};

/** Naive way
  */
template<typename T, typename Lattice>
class NaiveMomentumExchange : public TwoWayHelperFunctional<T, Lattice> {
public:
  /// Constructor
  NaiveMomentumExchange ( UnitConverter<T, Lattice>& converter,
                          SuperLattice<T, Lattice>& sLattice,
                          std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity );
  /// Computes the momentum transfer from fluid to particle.
  virtual bool operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic ) override;
};

/** Using Ladd mechanism
  */
template<typename T, typename Lattice>
class LaddMomentumExchange : public TwoWayHelperFunctional<T, Lattice> {
public:
  /// Constructor
  LaddMomentumExchange ( UnitConverter<T, Lattice>& converter,
                         SuperLattice<T, Lattice>& sLattice,
                         std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity,
                         std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > interpLatticeVelocity );
  /// Computes the momentum transfer from fluid to particle.
  virtual bool operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic ) override;
};

}

#endif
