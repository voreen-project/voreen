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
 *  -- header file
 */
#ifndef REACTING_SPECIES_3D_H
#define REACTING_SPECIES_3D_H

namespace olb {

/*
 * Base virtual class. Its raison d'etre is that it is not templeted.
 */
template <typename T>
class ReactingSpeciesBase3D {
protected:
  // Hidden constructor
  ReactingSpeciesBase3D(T stoichioCoeff);
public:
  // Get the stoichiometric coefficient
  T getStoichioCoeff();
  // Read the field value at a input lattice site within the input BlockLattice.
  virtual T getField(int iX, int iY, int iZ) =0;
  // Read the source value at a input lattice site within the input BlockLattice.
  virtual T getSource(int iX, int iY, int iZ) =0;
  // Write the source value at a input lattice site within the input BlockLattice.
  virtual void setSource(T val, int iX, int iY, int iZ) =0;
private:
  // Stoichiometric coefficient. If negative, the species is a reagent; if positive, a product.
  T _stoichioCoeff;
};

/*
 * This class provides a unified set of reading-writing methods to reacting specie fields.
 * The methods work the same to the user regardless of whether the specie field is solved through explicit Finite-Difference
 * (thereby being defined as an external field), or through Lattice-Boltzmann (thereby being defined as the 0-th momenum of
 * Lattice-Boltzmann one-particle density functions). This is achieved through ad-hoc specializations.
 * The Finite-Difference implementation is flexible enough to allow the external field being defined on both the main lattice,
 * and on a coupled lattice. Conversely, the Lattice-Boltzmann implementation allows only a coupled lattice.
 */
template <typename T, typename DESCRIPTOR, typename SOURCE>
class ReactingSpecies3D : public ReactingSpeciesBase3D<T> {
public:
  // Constructor
  ReactingSpecies3D(SuperGeometry<T,3>& superGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, T stoichioCoeff);
  // Read the source value at a input lattice site within the input BlockLattice.
  virtual T getSource(int iX, int iY, int iZ) override;
  // Write the source value at a input lattice site within the input BlockLattice.
  virtual void setSource(T val, int iX, int iY, int iZ) override;
protected:
  BlockLattice<T,DESCRIPTOR>* _blockLattice;
};

template <typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
class FiniteDifferenceReactingSpecies3D final : public ReactingSpecies3D<T,DESCRIPTOR,SOURCE> {
public:
  // Constructor
  FiniteDifferenceReactingSpecies3D(SuperGeometry<T,3>& superGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, T stoichioCoeff, std::size_t& iT);
  // Read the field value at a input lattice site within the input BlockLattice.
  // Finite-Difference: access external field.
  virtual T getField(int iX, int iY, int iZ) override;
private:
  // Reference to the simulation's timestep in order to take track of odd / even timesteps
  std::size_t& _iT;
};

template <typename T, typename DESCRIPTOR, typename SOURCE>
class LatticeBoltzmannReactingSpecies3D final : public ReactingSpecies3D<T,DESCRIPTOR,SOURCE> {
public:
  // Constructor
  LatticeBoltzmannReactingSpecies3D(SuperGeometry<T,3>& superGeometry, SuperLattice<T,DESCRIPTOR>& sLattice, T stoichioCoeff);
  // Read the field value at a input lattice site within the input BlockLattice.
  // Lattice-Boltzmann: compute the density according to the dynamics specified to the lattice.
  virtual T getField(int iX, int iY, int iZ) override;
};

};  // namespace olb

#endif
