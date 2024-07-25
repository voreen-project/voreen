/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
 *                2021 Clara Schragmann, Adrian Kummerlaender
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
 * Representation of the 2D block geometry view -- header file.
 */

#ifndef BLOCK_GEOMETRY_H
#define BLOCK_GEOMETRY_H

#include <vector>
#include <list>

#include "core/blockStructure.h"
#include "core/fieldArrayD.hh"
#include "geometry/blockGeometryStatistics2D.h"
#include "geometry/cuboid2D.h"
#include "communication/communicatable.h"

// All OpenLB code is contained in this namespace.
namespace olb {


/// Representation of a block geometry
/**
 * This class is derived from block geometry structure. It
 * holds the actual data with the materials. It stores pointers
 * to all dependent block geometry views.
 * It presents a volume of voxels where different types are
 * given my material numbers which is important e.g. to work
 * with different boundaries (like for inflow/output regions).
 */
template <typename T, unsigned D>
class BlockGeometry final : public BlockStructureD<D> {
private:
  /// Material number storage
  FieldArrayD<T,descriptors::SPATIAL_DESCRIPTOR<2>,Platform::CPU_SISD,descriptors::MATERIAL> _data;
  /// Material communicatable
  ConcreteCommunicatable<ColumnVector<cpu::sisd::Column<int>,1>> _communicatable;
  /// Cuboid which charaterizes the block geometry
  Cuboid<T,D> _cuboid;
  /// Number of the cuboid, default=-1
  int _iCglob;
  /// Statistic class
  BlockGeometryStatistics<T, D> _statistics;
  /// class specific output stream
  mutable OstreamManager clout;
  /// List to all depending statistic status objects
  std::list<bool*> _statisticsUpdateNeeded;

public:
  static constexpr Platform platform = Platform::CPU_SISD;

  BlockGeometry(Cuboid<T,D>& cuboid, int padding, int iCglob=-1);

  Platform getPlatform() const {
    return platform;
  }

  /// Write access to the associated block statistic
  BlockGeometryStatistics<T,D>& getStatistics(bool verbose=true);
  /// Read only access to the associated block statistic
  BlockGeometryStatistics<T,D> const& getStatistics(bool verbose=true) const;

  bool hasCommunicatable(std::type_index field) const {
    return field == typeid(descriptors::MATERIAL);
  }
  auto& getCommunicatable(std::type_index field) {
    OLB_ASSERT(field == typeid(descriptors::MATERIAL),
               "BlockGeometry only offers MATERIAL for communication");
    return _communicatable;
  }

  /// Read only access to the global iC number which is given !=-1 if the block geometries are part of a super geometry
  int const& getIcGlob() const;
  /// Returns the extend of the block in lattice units
  Vector<int,D> getExtent() const;

  /// Read only access to the origin position given in SI units (meter)
  Vector<T,D> getOrigin() const;
  /// Read only access to the voxel size given in SI units (meter)
  T getDeltaR() const;

  /// Write access to a material number
  int& get(LatticeR<D> latticeR);
  int& get(const int latticeR[D]);
  int& get(std::size_t iCell);

  template <typename... L>
  std::enable_if_t<sizeof...(L) == D, int&>
  get(L... latticeR) {
    return this->get(LatticeR<D>{latticeR...});
  }

  /// Read only access to a material number
  int get(LatticeR<D> latticeR) const;
  int get(const int latticeR[D]) const;
  int get(std::size_t iCell) const;

  template <typename... L>
  std::enable_if_t<sizeof...(L) == D, int>
  get(L... latticeR) const {
    return this->get(LatticeR<D>{latticeR...});
  }

  /// returns the (iX,iY) entry in the 2D scalar field
  int getMaterial(LatticeR<D> latticeR) const; // TODO old
  int getMaterial(const int latticeR[D]) const;

  template <typename... L>
  std::enable_if_t<sizeof...(L) == D, int>
  getMaterial(L... latticeR) const {
    return this->getMaterial(LatticeR<D>{latticeR...});
  }

  Vector<T,D> getPhysR(LatticeR<D> latticeR) {
    T physR[D];
    getPhysR(physR, latticeR);
    return Vector<T,D>(physR);
  }
  /// Transforms lattice to physical coordinates (wrapped from cuboid geometry)
  void getPhysR(T physR[D], const int latticeR[D]) const;
  void getPhysR(T physR[D], LatticeR<D> latticeR) const;

  template <typename... L>
  std::enable_if_t<sizeof...(L) == D, void>
  get(T physR[D],L... latticeR) const {
    return this->getPhysR(physR[D],LatticeR<D>{latticeR...});
  }

  /// Changes all cell materials which are not in bulkMaterials to 0 if
  /// there is no neighbour from bulkMaterials
  int clean(bool verbose=true, std::vector<int> bulkMaterials={1});
  /// Changes all cell materials from bulkMaterials to 0 if there is a neighbour with material 0
  int outerClean(bool verbose=true, std::vector<int> bulkMaterials={1});
  /// Changes all cell materials which are not 0 or 1 to 1 if there is a non robust constiallation
  int innerClean(bool verbose=true);
  /// Changes all cells with material fromM to 1 if there is a non robust constiallation
  int innerClean(int fromM, bool verbose=true);

  /// Resets all cell materials inside of a domain to 0
  void reset(IndicatorF<T,D>& domain);

  /// Returns the coordinates (iX,iY) of a voxel with a given material number (material) if there exists an neighbourhood of size (offsetX,offsetY) only with voxels of the  given material number
  bool find(int material, std::vector<unsigned> offset, std::vector<int> var);

  /// Returns true if at position (iX,iY) and in a neighbourhood of size (offsetX,offsetY) only voxels with a given material number (material) are there
  bool check(int material, std::vector<int> var, std::vector<unsigned> offset);
  /// Checks for errors (searches for all outer voxels (=0) with an inner voxel (=1) as a direct neighbour)
  bool checkForErrors(bool verbose=true) const;

  /// Replaces all material numbers (fromM) to another (toM)
  void rename(int fromM, int toM);
  /// Replaces all material numbers (fromM) to another (toM) if an indicator functor condition is fulfilled
  void rename(int fromM, int toM, IndicatorF<T,D>& condition);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX-offsetX,..,iX,..,ix+offsetX), .. are of the original material number (fromM)
  void rename(int fromM, int toM, LatticeR<D> offset);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX+1,iX+2,..,ix+testDirection[0]), .. are of another material number (testM)
  void rename(int fromM, int toM, int testM, std::vector<int> testDirection);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX+discreteNormal[0],iX+2*discreteNormal[0]), .. are of another material number (testM) and if an indicator functor condition is fulfilled
  void rename(int fromM, int toM, int fluidM, IndicatorF<T,D>& condition, Vector<int,D> discreteNormal);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX+discreteNormal[0],iX+2*discreteNormal[0]), .. are of another material number (fluidM) and if an indicator functor condition is fulfilled, the discreteNormal is computed from all fromM which fulfill the indicator functor condition
  void rename(int fromM, int toM, int fluidM, IndicatorF<T,D>& condition);

  /// Copy a layer of material numbers inside an indicator in a discrete normal direction
  void copyMaterialLayer(IndicatorF3D<T>& condition, int discreteNormal[D], int numberOfLayers);

  /// Replaces all material numbers (fromM) to another (toM) using a seed point and max. directions indicated by offsetX,Y != 0
  void regionGrowing(int fromM, int toM, LatticeR<D> seed, std::vector<int> offset, std::map<std::vector<int>, int >* tmp=nullptr);

  /// Prints a chosen part of the block geometry
  void printLayer(std::vector<int> min, std::vector<int> max, bool linenumber = false);
  /// Prints a chosen part of the block geometry
  void printLayer(int direction, int layer, bool linenumber = false);
  /// Prints a chosen node and its neighbourhood
  void printNode(std::vector<int> loc);

  /// Adds a pointer to the list of dependent statistic classes
  void addToStatisticsList(bool* statisticStatus);
  /// Removes a pointer from the list of dependent statistic classes if existing
  void removeFromStatisticsList(bool* statisticStatus);

private:
  /// Resets all depending statistic flags
  void resetStatistics();

};

template <typename T>
using BlockGeometry2D = BlockGeometry<T,2>;

template <typename T>
using BlockGeometry3D = BlockGeometry<T,3>;

/// Curried BlockGeometry template for use in callUsingConcretePlatform
template<typename T, unsigned D>
struct ConcretizableBlockGeometry {

using base_t = BlockGeometry<T,D>;

template <Platform PLATFORM>
using type = BlockGeometry<T,D>;

};

} // namespace olb

#endif
