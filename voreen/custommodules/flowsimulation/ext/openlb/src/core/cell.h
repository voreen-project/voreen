/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, 2015 Mathias J. Krause
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
 * Definition of a LB cell -- header file.
 */
#ifndef CELL_H
#define CELL_H

#include "olbDebug.h"
#include "serializer.h"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/dynamics.h"

namespace olb {

/// A LB lattice cell.
/** A cell contains the q values of the distribution functions f on
 * one lattice point, as well as a pointer to the dynamics of the
 * cell. Thanks to this pointer, one can have a space dependend de-
 * finition of the dynamics. This mechanism is useful e.g. for the
 * implementation of boundary conditions, or an inhomogeneous body
 * force.
 *
 * The dynamics object is not owned by the class, it is not
 * destructed in the Cell destructor.
 *
 * This class is not intended to be derived from.
 */
template<typename T, class Descriptor>
class CellBase {

protected:
  /// The lattice populations are defined as a q-element C-array.
  T f[Descriptor::q];   ///< distribution functions

public:
  /// Read-write access to distribution functions.
  /** \param iPop index of the accessed distribution function */
  T& operator[](int const& iPop)
  {
    OLB_PRECONDITION( iPop < Descriptor::q );
    return f[iPop];
  }
  /// Read-only access to distribution functions.
  /** \param iPop index of the accessed distribution function */
  T const& operator[](int const& iPop) const
  {
    OLB_PRECONDITION( iPop < Descriptor::q );
    return f[iPop];
  }
};


template<typename T, template<typename U> class Lattice>
class Cell : public CellBase<T, typename Lattice<T>::BaseDescriptor>, public Serializable {
public:
  /// Additional per-cell scalars for external fields, e.g. forces
  typedef descriptors::ExternalFieldArray <
  T, typename Lattice<T>::ExternalField > External;
private:
  External             external;  ///< external scalars
  Dynamics<T,Lattice>* dynamics;  ///< local LB dynamics

public:
  /// Default constructor.
  Cell();
  /// Constructor, to be used whenever possible.
  Cell(Dynamics<T,Lattice>* dynamics_);
public:
  /// Get a pointer to an external field
  T* getExternal(int offset)
  {
    OLB_PRECONDITION( offset < Lattice<T>::ExternalField::numScalars );
    return external.get(offset);
  }
  /// Get a const pointer to an external field
  T const* getExternal(int offset) const
  {
    OLB_PRECONDITION( offset < Lattice<T>::ExternalField::numScalars );
    return external.get(offset);
  }
  /// Define or re-define dynamics of the cell.
  /** \param dynamics_ a pointer to the dynamics object, whos
    *    memory management falls under the responsibility of the
    *    user */
  void defineDynamics(Dynamics<T,Lattice>* dynamics_);
  /// Get a non-modifiable pointer to the dynamics
  Dynamics<T,Lattice> const* getDynamics() const;
  /// Get a non-modifiable pointer to the dynamics
  Dynamics<T,Lattice>* getDynamics();

  // The following helper functions forward the function call
  // to the Dynamics object
public:
  /// Apply LB collision to the cell according to local dynamics.
  void collide(LatticeStatistics<T>& statistics)
  {
    OLB_PRECONDITION( dynamics );
    dynamics->collide(*this, statistics);
  }
  /// Apply LB collision with fixed velocity to the cell.
  void staticCollide(const T u[Lattice<T>::d], LatticeStatistics<T>& statistics)
  {
    OLB_PRECONDITION( dynamics );
    dynamics->staticCollide(*this, u, statistics);
  }

  /// Compute particle density on the cell.
  /** \return particle density
   */
  T computeRho() const
  {
    OLB_PRECONDITION( dynamics );
    return dynamics->computeRho(*this);
  }
  /// Compute fluid velocity on the cell.
  /** \param u fluid velocity
   */
  void computeU(T u[Lattice<T>::d]) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeU(*this, u);
  }
  /// Compute components of the stress tensor on the cell.
  /** \param pi stress tensor */
  void computeStress (
    T pi[util::TensorVal<Lattice<T> >::n]) const
  {
    OLB_PRECONDITION( dynamics );
    T rho, u[Lattice<T>::d];
    dynamics->computeRhoU(*this, rho, u);
    dynamics->computeStress(*this, rho, u, pi);
  }
  /// Compute fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  void computeRhoU(T& rho, T u[Lattice<T>::d]) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeRhoU(*this, rho, u);
  }
  /// Compute all momenta on the celll, up to second order.
  /** \param rho particle density
   *  \param u fluid velocity
   *  \param pi stress tensor
   */
  void computeAllMomenta (
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeAllMomenta(*this, rho, u, pi);
  }
  /// Access external fields through the dynamics object.
  /** This method is similar to getExternal(): it delivers the
   * value of the external fields. This time, those values
   * are however computed through a virtual call to the dynamics
   * object.
   */
  void computeExternalField(int pos, int size, T* ext) const
  {
    OLB_PRECONDITION(pos+size <= Lattice<T>::ExternalField::numScalars);
    T const* externalData = this->getExternal(pos);
    for (int iExt=0; iExt<size; ++iExt) {
      ext[iExt] = externalData[iExt];
    }
  }
  /// Set particle density on the cell.
  /** \param rho particle density
   */
  void defineRho(T rho)
  {
    OLB_PRECONDITION( dynamics );
    dynamics->defineRho(*this, rho);
  }
  /// Set fluid velocity on the cell.
  /** \param u fluid velocity
   */
  void defineU(const T u[Lattice<T>::d])
  {
    OLB_PRECONDITION( dynamics );
    dynamics->defineU(*this, u);
  }
  /// Define fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  void defineRhoU(T rho, const T u[Lattice<T>::d])
  {
    OLB_PRECONDITION( dynamics );
    dynamics->defineRhoU(*this, rho, u);
  }
  /// Define particle populations through the dynamics object.
  /** This method is similar to operator[]: it modifies the
   * value of all the particle populations.
   */
  void definePopulations(const T* f_)
  {
    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
    {
      this->f[iPop] = f_[iPop];
    }
  }
  /// Define external fields through the dynamics object.
  /** This method is similar to getExternal(): it accesses the
   * value of the external fields.
   */
  void defineExternalField(int pos, int size, const T* ext)
  {
    OLB_PRECONDITION(pos+size <= Lattice<T>::ExternalField::numScalars);
    T* externalData = this->getExternal(pos);
    for (int iExt=0; iExt<size; ++iExt) {
      externalData[iExt] = ext[iExt];
    }
  }
  /// Add external fields through the dynamics object.
  /** Similar to defineExternalField(),but instead of replacing existing values
   *  ext is added to existing values.
   */
  inline void addExternalField(int pos, int size, const T* ext)
  {
    OLB_PRECONDITION(pos+size <= Lattice<T>::ExternalField::numScalars);
    T* externalData = this->getExternal(pos);
    for (int iExt=0; iExt<size; ++iExt) {
      externalData[iExt] += ext[iExt];
    }
  }
  /// Add external fields through the dynamics object.
  /** Similar to defineExternalField(),but instead of replacing existing values
   *  ext is multiplied to existing values.
   */
  inline void multiplyExternalField(int pos, int size, const T* ext)
  {
    OLB_PRECONDITION(pos+size <= Lattice<T>::ExternalField::numScalars);
    T* externalData = this->getExternal(pos);
    for (int iExt=0; iExt<size; ++iExt) {
      externalData[iExt] *= ext[iExt];
    }
  }
  /// Initialize all f values to their local equilibrium
  void iniEquilibrium(T rho, const T u[Lattice<T>::d])
  {
    OLB_PRECONDITION( dynamics );
    dynamics->iniEquilibrium(*this, rho, u);
  }
  /// Revert ("bounce-back") the distribution functions.
  void revert();
  void serialize(T* data) const;
  void unSerialize(T const* data);

  /// \return the number of data blocks for the serializable interface
  std::size_t getNblock() const override
  {
    if (Lattice<T>::ExternalField::numScalars)
      return 2;
    else
      return 1;
  };
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// \return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

private:
  void iniPop();
  void iniExternal();
};

template<typename T, template<typename U> class Lattice>
struct WriteCellFunctional {
  virtual ~WriteCellFunctional() { };
  virtual void apply(Cell<T,Lattice>& cell, int pos[Lattice<T>::d]) const =0;
};

}  // namespace olb

#endif
