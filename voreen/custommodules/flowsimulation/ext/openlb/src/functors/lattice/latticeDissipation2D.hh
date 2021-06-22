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

#ifndef LATTICE_DISSIPATION_2D_HH
#define LATTICE_DISSIPATION_2D_HH

#include <vector>
#include <cmath>
#include <limits>

#include "latticeDissipation2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry2D.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "core/blockLattice2D.h"
#include "communication/mpiManager.h"
#include "core/blockLatticeStructure2D.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticeDissipation2D<T,DESCRIPTOR>::SuperLatticeDissipation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDissipation2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeDissipation2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{

  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  ////  int lociz= input[3];
  //  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
  //    // local coords are given, fetch local cell and compute value(s)
  //    T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  //    int overlap = this->sLattice.getOverlap();
  //    this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeAllMomenta(rho, uTemp, pi);

  //    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  //    if (util::TensorVal<DESCRIPTOR >::n == 6)
  //      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  //    T nuLattice = converter.getLatticeNu();
  //    T omega = converter.getOmega();
  //    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*descriptors::invCs2<T,DESCRIPTOR>(),2)/rho/2.;

  //    return std::vector<T>(1,finalResult);
  //  } else {
  //    return std::vector<T>(); // empty vector
  //  }
  return false;
}

template <typename T, typename DESCRIPTOR>
BlockLatticeDissipation2D<T,DESCRIPTOR>::BlockLatticeDissipation2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _converter(converter)
{
  this->getName() = "dissipation";
}

// todo: get functor working
template <typename T, typename DESCRIPTOR>
bool BlockLatticeDissipation2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    // local coords are given, fetch local cell and compute value(s)
  //    T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  //    int overlap = this->_blockLattice.getOverlap();
  //    this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeAllMomenta(rho, uTemp, pi);

  //    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  //    if (util::TensorVal<DESCRIPTOR >::n == 6)
  //      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  //    T nuLattice = converter.getLatticeNu();
  //    T omega = converter.getOmega();
  //    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*descriptors::invCs2<T,DESCRIPTOR>(),2)/rho/2.;

  //    return std::vector<T>(1,finalResult);
  //  } else {
  //    return std::vector<T>(); // empty vector
  //  }

  return false;
}

}
#endif
