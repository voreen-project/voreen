/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_PSM_PHYS_FORCE_2D_HH
#define LATTICE_PSM_PHYS_FORCE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePSMPhysForce2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticePSMPhysForce2D<T, DESCRIPTOR>::SuperLatticePSMPhysForce2D(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "PSMPhysForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePSMPhysForce2D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePSMPhysForce2D<T,DESCRIPTOR>::BlockLatticePSMPhysForce2D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2)
{
  this->getName() = "physPSMForce";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePSMPhysForce2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  T epsilon = 1. - this->_blockLattice.get(input[0], input[1]).template getField<descriptors::POROSITY>();

  //if ((epsilon > 1e-5 && epsilon < 1 - 1e-5)) 
  if ((epsilon > 1e-5)) {
    T rho, u[DESCRIPTOR::d], u_s[DESCRIPTOR::d];

    for (int i = 0; i < DESCRIPTOR::d; i++) {
      u_s[i] = this->_blockLattice.get(input[0], input[1]).template getFieldComponent<descriptors::VELOCITY_SOLID>(i);
    }
    T paramA = this->_converter.getLatticeRelaxationTime() - 0.5;
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);

    T omega_s;
    T omega = 1. / this->_converter.getLatticeRelaxationTime();

    this->_blockLattice.get(input[0], input[1]).computeRhoU(rho, u);

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s);
    T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      //switch (mode) {
      //case M2:
        omega_s = (lbm< DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s)
                   - this->_blockLattice.get(input[0], input[1])[iPop])
                  + (1 - omega)
                  * (this->_blockLattice.get(input[0], input[1])[iPop]
                     - lbm< DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
      // break;
      //case M3:
      //  omega_s =
      //    (this->_blockLattice.get(input[0], input[1])[descriptors::opposite<DESCRIPTOR>(iPop)]
      //     - lbm< DESCRIPTOR>::equilibrium(
      //       descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s, uSqr_s))
      //    - (this->_blockLattice.get(input[0], input[1])[iPop]
      //       - lbm< DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s));
      //}

      for (int i = 0; i < this->getTargetDim(); ++i) {
        output[i] -= descriptors::c<DESCRIPTOR>(iPop,i) * omega_s;
      }
    }

    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = this->_converter.getPhysForce(output[i] * paramB);
    }
  }
  return true;
}



//TODO: 200116 preliminary version: proper indicatorcheck to be included
template<typename T, typename DESCRIPTOR>
SuperLatticePSMPhysForce2DMod<T, DESCRIPTOR>::SuperLatticePSMPhysForce2DMod(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  SuperGeometry<T,2>& superGeometry,
  const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF2D<T>& shapeIndicator )
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2) 
{
  this->getName() = "PSMPhysForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePSMPhysForce2DMod<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        superGeometry.getBlockGeometry(iC),
        this->_converter,
        shapeIndicator));
  }
}


//TODO: 200116 preliminary version: proper indicatorcheck to be included
template <typename T, typename DESCRIPTOR>
BlockLatticePSMPhysForce2DMod<T,DESCRIPTOR>::BlockLatticePSMPhysForce2DMod(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
                             BlockGeometry<T,2>& blockGeometry,
  const UnitConverter<T,DESCRIPTOR>&     converter,
  IndicatorF2D<T>& shapeIndicator )
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2), _shapeIndicator(shapeIndicator), _blockGeometry(blockGeometry)
{
  this->getName() = "physPSMForce";
}

//TODO: 200116 preliminary version: proper indicatorcheck to be included
template<typename T, typename DESCRIPTOR>
bool BlockLatticePSMPhysForce2DMod<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  T epsilon = 1. - this->_blockLattice.get(input[0], input[1]).template getField<descriptors::POROSITY>();
  if ((epsilon > 1e-5)) {

    T posBlockGeo[2] = {0.};
    bool insideShape[1];
    _blockGeometry.getPhysR( posBlockGeo, {input[0], input[1]});
    _shapeIndicator( insideShape, posBlockGeo );

    if ( insideShape[0] == true ) {

      T rho, u[DESCRIPTOR::d], u_s[DESCRIPTOR::d];

      for (int i = 0; i < DESCRIPTOR::d; i++) {
        u_s[i] = this->_blockLattice.get(input[0], input[1]).template getFieldComponent<descriptors::VELOCITY_SOLID>(i);
      }
      T paramA = this->_converter.getLatticeRelaxationTime() - 0.5;
      // speed up paramB
      T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);

      T omega_s;
      T omega = 1. / this->_converter.getLatticeRelaxationTime();

      this->_blockLattice.get(input[0], input[1]).computeRhoU(rho, u);

      const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s);
      T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        //switch (mode) {
        //  case M2:
            omega_s = (lbm< DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s)
                - this->_blockLattice.get(input[0], input[1])[iPop])
                + (1 - omega)
                    * (this->_blockLattice.get(input[0], input[1])[iPop]
                        - lbm< DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
        //  break;
        //  case M3:
        //    omega_s =
        //        (this->_blockLattice.get(input[0], input[1])[descriptors::opposite<DESCRIPTOR>(iPop)]
        //            - lbm< DESCRIPTOR>::equilibrium(
        //                descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s, uSqr_s))
        //            - (this->_blockLattice.get(input[0], input[1])[iPop]
        //                - lbm< DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s));
        //}

        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] -= descriptors::c<DESCRIPTOR>(iPop,i) * omega_s;
        }
      }

      for (int i = 0; i < this->getTargetDim(); ++i) {
        output[i] = this->_converter.getPhysForce(output[i] * paramB);
      }
    }
  return true;
  }
}



}
#endif
