/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause, Peter Weisbrod
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
 * Representation of a parallel 2D geometry -- header file.
 */

/// A super geometry represents a parrallel voxel mesh
/** A super geometry consits of a number of block geometries,
 * where the material numbers are stored. It is constructed
 * from a cuboid geometry. All coboids of the cuboid geometry
 * are asigned to block geometries which are extended by an
 * overlap in order to enable efficient parallelisation.
 *
 * By the class access is provied to the material numbers of
 * the mesh. Methods for renaming materials are provided as
 * well as a statistic class.
 *
 * This class is not intended to be derived from.
 */


#ifndef SUPER_GEOMETRY_H
#define SUPER_GEOMETRY_H

#include <vector>
#include <iostream>
#include <string>

#include "geometry/cuboidGeometry2D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/superGeometryStatistics2D.h"
#include "geometry/superGeometryStatistics3D.h"
#include "geometry/blockGeometry.h"
#include "communication/superStructure.h"
#include "communication/loadBalancer.h"
#include "communication/superCommunicator.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"


// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T, unsigned D>
class SuperGeometry : public SuperStructure<T,D> {
private:
  /// Vector of block geometries with overlap
  std::vector<std::unique_ptr<BlockGeometry<T,D>>> _block{};
  /// Communicator for propagation of populations between blocks
  std::unique_ptr<SuperCommunicator<T,SuperGeometry<T,D>>> _communicator{};
  /// True if changes need to be communicated
  bool _communicationNeeded{};
  /// Statistic class
  SuperGeometryStatistics<T,D> _statistics{};
  /// class specific output stream
  mutable OstreamManager clout;

public:
  constexpr static unsigned d = D;

  using block_t = ConcretizableBlockGeometry<T,D>;

  SuperGeometry(CuboidGeometry<T,D>& cuboidGeometry,
                LoadBalancer<T>& loadBalancer,
                int overlap = 3);

  /// Read only access to the material numbers, error handling: returns 0 if data is not available
  int const& get(int iCglob, LatticeR<D> latticeR) const;
  int const& get(const int latticeR[D+1]) const;
  int const& get(LatticeR<D+1> latticeR) const;

  template <typename... L>
  std::enable_if_t<sizeof...(L) == (D+1), int>
  get(L... latticeR) {
    return this->get(LatticeR<D+1>{latticeR...});
  }

  /// Read only access to the material numbers with global communication to all ranks
  int getAndCommunicate(int iCglob, LatticeR<D> latticeR) const;
  int getAndCommunicate(LatticeR<D+1> latticeR) const;
  /// Write access to the material numbers, error handling: stops the program if data is not available
  int& set(std::vector<int> latticeR); //TODO to be removed set->get, problem: with get calling wrong function


  /// Transforms a lattice to physical position (SI unites)
  std::vector<T> getPhysR(int iCglob, LatticeR<D> latticeR) const;
  /// Transforms a lattice to physical position (SI unites)
  std::vector<T> getPhysR(LatticeR<D+1> latticeR) const;
  /// Transforms a lattice to physical position (SI unites)
  void getPhysR(T output[D], const int latticeR[D+1]) const;
  void getPhysR(T output[D], const int iCglob, LatticeR<D> latticeR) const;

  /// Read and write access to a single block geometry
  BlockGeometry<T,D>& getBlockGeometry(int locIC);
  /// Read only access to a single block geometry
  BlockGeometry<T,D> const& getBlockGeometry(int locIC) const;

  /// Read and write access to a single extended block geometry
  template <typename BLOCK = BlockGeometry<T,D>>
  BLOCK& getBlock(int locIC);
  /// Read only access to a single extended block geometry
  template <typename BLOCK = BlockGeometry<T,D>>
  const BLOCK& getBlock(int locIC) const;

  /// Returns the statistics object
  SuperGeometryStatistics<T,D>& getStatistics();
  /// Read and write access to the statistic status flag, update needed = true
  bool& getStatisticsStatus();
  /// Read only access to the statistic status flag, update needed = true
  bool const& getStatisticsStatus() const;
  /// Updates the super geometry at the boundaries if needed and afterwards the statisics if needed
  void updateStatistics(bool verbose=true);

  /// Executes an outer cleaning: Sets all material numbers which are not
  /// bulk-materials to 0 if there is no neighbour from bulkMaterials
  int clean(bool verbose=true, std::vector<int> bulkMaterials={1});
  /// Removes not needed fluid cells from the outer domain
  int outerClean(bool verbose=true, std::vector<int> bulkMaterials={1});
  /// inner cleaning for all boundary types
  int innerClean(bool verbose=true);
  /// inner cleaning for specific boundary types
  int innerClean(int material, bool verbose=true);
  /// check for errors (searches for all outer voxels (=0) with an inner voxel (=1) as a direct neighbour)
  bool checkForErrors(bool verbose=true);

  /// reset all cell materials inside of a domain to 0
  void reset(IndicatorF<T,D>& domain);

  /// replace one material with another
  void rename(int fromM, int toM);
  /// replace one material that fulfills an indicator functor condition with another
  void rename(int fromM, int toM, FunctorPtr<IndicatorF<T,D>>&& condition);
  /// replace one material with another respecting an offset (overlap)
  void rename(int fromM, int toM, LatticeR<D> offset);
  /// renames all voxels of material fromM to toM if the number of voxels given by testDirection is of material testM
  void rename(int fromM, int toM, int testM, std::vector<int> testDirection);
  /// renames all boundary voxels of material fromBcMat to toBcMat if two neighbour voxel in the direction of the discrete normal are fluid voxel with material fluidM in the region where the indicator function is fulfilled
  void rename(int fromBcMat, int toBcMat, int fluidMat, IndicatorF<T,D>& condition);
  /// renames all boundary voxels of material fromBcMat to toBcMat if two neighbour voxel in the direction of the discrete normal are fluid voxel with material fluidM in the region where the indicator function is fulfilled
    void rename(int fromBcMat, int toBcMat, int fluidMat, FunctorPtr<IndicatorF<T,D>>&& condition);


  /// Prints some information about the super geometry
  void print();

  /**
   * Returns a material indicator using the given vector of materials
   *
   * \param  materials Materials to be indicated
   * \returns          Unique ownership of the constructed indicator.
   *                   May be stored or passed directly to e.g. defineDynamics
   **/
  std::unique_ptr<SuperIndicatorF<T,D>> getMaterialIndicator(std::vector<int>&& materials);
  /**
   * Returns a material indicator using a single material number
   *
   * \param material Material to be indicated
   * \returns        Unique ownership of the constructed indicator.
   *                 May be stored or passed directly to e.g. defineDynamics
   **/
  std::unique_ptr<SuperIndicatorF<T,D>> getMaterialIndicator(int material);

  void communicate() override
  {
    if (_communicationNeeded) {
      _communicator->communicate();
      _communicationNeeded = false;
    }
  }

};

} // namespace olb

#endif
