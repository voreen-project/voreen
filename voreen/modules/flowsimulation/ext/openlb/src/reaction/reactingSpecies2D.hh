/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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

/** \file
 * Classes to provide access to a reacting species - may be a LB density, or an external field.
 *  -- generic implementation
 */
#ifndef REACTING_SPECIES_2D_HH
#define REACTING_SPECIES_2D_HH

namespace olb {


///////////////////////////////////// class ReactingSpeciesBase2D /////////////////////////////////////

template <typename T>
ReactingSpeciesBase2D<T>::ReactingSpeciesBase2D(T stoichioCoeff)
  : _stoichioCoeff(stoichioCoeff)
{}

template <typename T>
T ReactingSpeciesBase2D<T>::getStoichioCoeff()
{
  return _stoichioCoeff;
}


///////////////////////////////////// class ReactingSpecies2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE>
ReactingSpecies2D<T,DESCRIPTOR,SOURCE>::
ReactingSpecies2D(SuperGeometry<T,2>& superGeometry, SuperLattice<T, DESCRIPTOR>& sLattice, T stoichioCoeff)
  : ReactingSpeciesBase2D<T>(stoichioCoeff)
{
  std::vector<BlockLattice<T,DESCRIPTOR>*> blockLattices;
  for (int iC = 0; iC < superGeometry.getLoadBalancer().size(); ++iC) {
    blockLattices.push_back ( &sLattice.getBlock(iC) );
  }
  _blockLattice = blockLattices[0];
}

template <typename T, typename DESCRIPTOR, typename SOURCE>
T ReactingSpecies2D<T,DESCRIPTOR,SOURCE>::
getSource(int iX, int iY)
{
  return this->_blockLattice->get(iX, iY).template getFieldPointer<SOURCE>()[0];
}

template <typename T, typename DESCRIPTOR, typename SOURCE>
void ReactingSpecies2D<T,DESCRIPTOR,SOURCE>::
setSource(T val, int iX, int iY)
{
  this->_blockLattice->get(iX, iY).template setField<SOURCE>(val);
}


///////////////////////////////////// class FiniteDifferenceReactingSpecies2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,FIELD,SOURCE>::
FiniteDifferenceReactingSpecies2D(SuperGeometry<T,2>& superGeometry, SuperLattice<T, DESCRIPTOR>& sLattice, T stoichioCoeff, std::size_t& iT)
  : ReactingSpecies2D<T,DESCRIPTOR,SOURCE>(superGeometry, sLattice, stoichioCoeff),
    _iT(iT)
{
  static_assert(DESCRIPTOR::template size<FIELD>()  == 2, "FIELD must have size 2." );
  static_assert(DESCRIPTOR::template size<SOURCE>() == 1, "SOURCE must have size 1.");
}

template <typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
T FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,FIELD,SOURCE>::
getField(int iX, int iY)
{
  return *fd::accessNew<T,DESCRIPTOR,FIELD>(this->_blockLattice->get(iX, iY), this->_iT);
}


///////////////////////////////////// class LatticeBoltzmannReactingSpecies2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE>
LatticeBoltzmannReactingSpecies2D<T,DESCRIPTOR,SOURCE>::
LatticeBoltzmannReactingSpecies2D(SuperGeometry<T,2>& superGeometry, SuperLattice<T, DESCRIPTOR>& sLattice, T stoichioCoeff)
  : ReactingSpecies2D<T,DESCRIPTOR,SOURCE>(superGeometry, sLattice, stoichioCoeff)
{}

template <typename T, typename DESCRIPTOR, typename SOURCE>
T LatticeBoltzmannReactingSpecies2D<T,DESCRIPTOR,SOURCE>::
getField(int iX, int iY)
{
  return this->_blockLattice->get(iX, iY).computeRho();
}


}  // namespace olb

#endif
