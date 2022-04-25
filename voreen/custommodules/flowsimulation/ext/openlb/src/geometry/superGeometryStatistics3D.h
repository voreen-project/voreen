/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause
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
 * Representation of a statistic for a parallel 3D geometry -- header file.
 */

#ifndef SUPER_GEOMETRY_STATISTICS_3D_H
#define SUPER_GEOMETRY_STATISTICS_3D_H


#include <map>
#include <string>
#include <vector>

#include "communication/mpiManager.h"
#include "geometry/superGeometry.h"
#include "io/ostreamManager.h"


// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a statistic for a parallel 3D geometry
/** A super geometry statistic computes different integral
 * values, like total number of different materials,
 * materials of any kind, min./max. physical position, of an
 * underlying super geometry.
 *
 * This class is not intended to be derived from.
 */

template<typename T, unsigned D>
class SuperGeometry;

template<typename T>
class SuperGeometryStatistics3D {
private:
  /// Points to the underlying data from which the statistics is taken
  const SuperGeometry<T,3>* _superGeometry;
  /// Specifies if an update is needed
  bool _statisticsUpdateNeeded;

  /// Number of different material numbers
  int _nMaterials;
  /// Mapping a material number to the number of this kind found in the super geometry
  std::map<int, std::size_t> _material2n;
  /// Mapping a material number to the min. physical position in each space direction
  std::map<int, olb::Vector<T, 3> > _material2min;
  /// Mapping a material number to the max. physical position in each space direction
  std::map<int, olb::Vector<T, 3> > _material2max;
  /// Componentwise min extension over all material numbers. Default is 0.
  olb::Vector<T, 3> _minOverMaterial;
  /// Componentwise maximal extension over all material numbers. Default is 0.
  olb::Vector<T, 3> _maxOverMaterial;

  /// class specific cout
  mutable OstreamManager clout;

  const SuperGeometryStatistics3D<T>* const_this;

public:

  /// Constructor
  SuperGeometryStatistics3D(SuperGeometry<T,3>* superGeometry);
  /// Copy constructor
  SuperGeometryStatistics3D(SuperGeometryStatistics3D const& rhs);
  /// Copy assignment
  SuperGeometryStatistics3D<T>& operator=(SuperGeometryStatistics3D const& rhs);

  /// Read and write access to a flag, which indicates if an uptate is needed (=true)
  bool& getStatisticsStatus();
  /// Read only access to a flag, which indicates if an uptate is needed (=true)
  bool const & getStatisticsStatus() const;

  /// Updates the statistics if it is really needed
  void update(bool verbose=false);

  /// Returns the number of different materials
  int getNmaterials();
  int getNmaterials() const;
  /// Returns the number of voxels for a given material number
  std::size_t getNvoxel(int material);
  std::size_t getNvoxel(int material) const;
  /// Returns the number of voxels with material!=0
  std::size_t getNvoxel();
  std::size_t getNvoxel() const;
  /// Returns the min. phys position in each direction corresponding to material number
  olb::Vector<T, 3>getMinPhysR(int material);
  olb::Vector<T, 3>getMinPhysR(int material) const;
  /// Returns the min. phys position in each direction corresponding to all non-zero material numbers
  olb::Vector<T, 3>getMinPhysR();
  olb::Vector<T, 3>getMinPhysR() const;
  /// Returns the max. phys position in each direction corresponding to material number
  olb::Vector<T, 3> getMaxPhysR(int material);
  olb::Vector<T, 3> getMaxPhysR(int material) const;
  /// Returns the max. phys position in each direction corresponding to all non-zero material numbers
  olb::Vector<T, 3> getMaxPhysR();
  olb::Vector<T, 3> getMaxPhysR() const;
  /// Returns the phys extend as length in each direction
  olb::Vector<T, 3> getPhysExtend(int material);
  std::vector<T> getPhysExtend(int material) const;
  /// Returns the phys radius as length in each direction
  olb::Vector<T, 3> getPhysRadius(int material);
  olb::Vector<T, 3> getPhysRadius(int material) const;
  /// Returns the center position
  olb::Vector<T, 3> getCenterPhysR(int material);
  olb::Vector<T, 3> getCenterPhysR(int material) const;
  /// Returns the boundary type which is characterized by a discrte normal (c.f. Zimny)
  olb::Vector<int, 3> getType(int iC, int iX, int iY, int iZ);
  olb::Vector<int, 3> getType(int iC, int iX, int iY, int iZ) const;

  /// Returns normal that points into the fluid for paraxial surfaces
  olb::Vector<T, 3> computeNormal (int material);
  olb::Vector<T, 3> computeNormal (int material) const;
  /// Returns discrete normal with norm maxNorm that points into the fluid for paraxial surfaces
  /// maxNorm=1.1 implies only normals parallel to the axises
  olb::Vector<int, 3> computeDiscreteNormal (int material, T maxNorm = 1.1);
  olb::Vector<int, 3> computeDiscreteNormal (int material, T maxNorm = 1.1) const;
  /// Returns util::sqrt( maxX^2 + maxY^2 + maxZ^2 ) max over a certain material number
  T computeMaxPhysDistance(int material ) const;
  /// Returns util::sqrt( maxX^2 + maxY^2 + maxZ^2 ) max over all material numbers
  T computeMaxPhysDistance() const;
  /// Prints some statistic information, i.e. the number of voxels and min. max. physical position for each different material
  void print();
  void print() const;
};

} // namespace olb

#endif

