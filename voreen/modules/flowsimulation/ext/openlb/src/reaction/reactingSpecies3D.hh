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
#ifndef REACTING_SPECIES_3D_HH
#define REACTING_SPECIES_3D_HH

namespace olb {


///////////////////////////////////// class ReactingSpeciesBase3D /////////////////////////////////////

template <typename T>
ReactingSpeciesBase3D<T>::ReactingSpeciesBase3D(T stoichioCoeff)
  : _stoichioCoeff(stoichioCoeff)
{}

template <typename T>
T ReactingSpeciesBase3D<T>::getStoichioCoeff()
{
  return _stoichioCoeff;
}


///////////////////////////////////// class ReactingSpecies3D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE>
ReactingSpecies3D<T,DESCRIPTOR,SOURCE>::
ReactingSpecies3D(SuperGeometry<T,3>& superGeometry, SuperLattice<T, DESCRIPTOR>& sLattice, T stoichioCoeff)
  : ReactingSpeciesBase3D<T>(stoichioCoeff)
{
  std::vector<BlockLattice<T,DESCRIPTOR>*> blockLattices;
  for (int iC = 0; iC < superGeometry.getLoadBalancer().size(); ++iC) {
    blockLattices.push_back ( &sLattice.getBlock(iC) );
  }
  _blockLattice = blockLattices[0];
}

template <typename T, typename DESCRIPTOR, typename SOURCE>
T ReactingSpecies3D<T,DESCRIPTOR,SOURCE>::
getSource(int iX, int iY, int iZ)
{
  return this->_blockLattice->get(iX, iY, iZ).template getFieldPointer<SOURCE>()[0];
}

template <typename T, typename DESCRIPTOR, typename SOURCE>
void ReactingSpecies3D<T,DESCRIPTOR,SOURCE>::
setSource(T val, int iX, int iY, int iZ)
{
  this->_blockLattice->get(iX, iY, iZ).template setField<SOURCE>(val);
}


///////////////////////////////////// class FiniteDifferenceReactingSpecies3D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FiniteDifferenceReactingSpecies3D<T,DESCRIPTOR,FIELD,SOURCE>::
FiniteDifferenceReactingSpecies3D(SuperGeometry<T,3>& superGeometry, SuperLattice<T, DESCRIPTOR>& sLattice, T stoichioCoeff, std::size_t& iT)
  : ReactingSpecies3D<T,DESCRIPTOR,SOURCE>(superGeometry, sLattice, stoichioCoeff),
    _iT(iT)
{
  static_assert(DESCRIPTOR::template size<FIELD>()  == 2, "FIELD must have size 2." );
  static_assert(DESCRIPTOR::template size<SOURCE>() == 1, "SOURCE must have size 1.");
}

template <typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
T FiniteDifferenceReactingSpecies3D<T,DESCRIPTOR,FIELD,SOURCE>::
getField(int iX, int iY, int iZ)
{
  return *fd::accessNew<T,DESCRIPTOR,FIELD>(this->_blockLattice->get(iX, iY, iZ), this->_iT);
}


///////////////////////////////////// class LatticeBoltzmannReactingSpecies3D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE>
LatticeBoltzmannReactingSpecies3D<T,DESCRIPTOR,SOURCE>::
LatticeBoltzmannReactingSpecies3D(SuperGeometry<T,3>& superGeometry, SuperLattice<T, DESCRIPTOR>& sLattice, T stoichioCoeff)
  : ReactingSpecies3D<T,DESCRIPTOR,SOURCE>(superGeometry, sLattice, stoichioCoeff)
{}

template <typename T, typename DESCRIPTOR, typename SOURCE>
T LatticeBoltzmannReactingSpecies3D<T,DESCRIPTOR,SOURCE>::
getField(int iX, int iY, int iZ)
{
  return this->_blockLattice->get(iX, iY, iZ).computeRho();
}


}  // namespace olb

#endif
