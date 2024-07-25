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
 * Representation of a statistic for a parallel 2D geometry -- header file.
 */

#ifndef SUPER_GEOMETRY_STATISTICS_2D_H
#define SUPER_GEOMETRY_STATISTICS_2D_H


#include <map>
#include <string>
#include <vector>

#include "geometry/superGeometry.h"
#include "io/ostreamManager.h"


// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a statistic for a parallel 2D geometry
/** A super geomety statistic computes different integral
 * values, like total number of different materials,
 * materials of any kind, min./max. physical position, of an
 * underlying super geoemtry.
 *
 * This class is not intended to be derived from.
 */

template<typename T, unsigned D>
class SuperGeometry;

template<typename T>
class SuperGeometryStatistics2D {
private:
  /// Points to the underlying data from which the statistics is taken
  const SuperGeometry<T,2>* _superGeometry;
  /// Specifies if an update is needed
  bool _statisticsUpdateNeeded;

  /// Number of different material numbers
  int _nMaterials;
  /// Mapping a material number to the number of this kind found in the super geometry
  std::map<int, int> _material2n;
  /// Mapping a material number to the min. physical position in each space direction
  std::map<int, std::vector<T> > _material2min;
  /// Mapping a material number to the max. physical position in each space direction
  std::map<int, std::vector<T> > _material2max;

  /// class specific cout
  mutable OstreamManager clout;

  const SuperGeometryStatistics2D<T>* const_this;

public:

  /// Constructor
  SuperGeometryStatistics2D(SuperGeometry<T,2>* superGeometry);
  /// Copy constructor
  SuperGeometryStatistics2D(SuperGeometryStatistics2D const& rhs);
  /// Copy assignment
  SuperGeometryStatistics2D<T>& operator=(SuperGeometryStatistics2D const& rhs);

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
  int getNvoxel(int material);
  int getNvoxel(int material) const;
  /// Returns the number of voxels with material!=0
  int getNvoxel();
  int getNvoxel() const;
  /// Returns the min. phys position in each direction
  std::vector<T> getMinPhysR(int material);
  std::vector<T> getMinPhysR(int material) const;
  /// Returns the max. phys position in each direction
  std::vector<T> getMaxPhysR(int material);
  std::vector<T> getMaxPhysR(int material) const;
  /// Returns the phys extend as length in each direction
  std::vector<T> getPhysExtend(int material);
  std::vector<T> getPhysExtend(int material) const;
  /// Returns the phys radius as length in each direction
  std::vector<T> getPhysRadius(int material);
  std::vector<T> getPhysRadius(int material) const;
  /// Returns the center position
  std::vector<T> getCenterPhysR(int material);
  std::vector<T> getCenterPhysR(int material) const;
  /// Returns the boundary type which is characterized by a discrte normal (c.f. Zimny)
  std::vector<int> getType(int iC, int iX, int iY);
  std::vector<int> getType(int iC, int iX, int iY) const;

  /// Returns normal that points into the fluid for paraxial surfaces
  std::vector<T> computeNormal (int material);
  std::vector<T> computeNormal (int material) const;
  /// Returns discrete normal with norm maxNorm that points into the fluid for paraxial surfaces
  /// maxNorm=1.1 implies only normals parallel to the axises
  std::vector<int> computeDiscreteNormal (int material, T maxNorm = 1.1);
  std::vector<int> computeDiscreteNormal (int material, T maxNorm = 1.1) const;

  /// Prints some statistic information, i.e. the number of voxels and min. max. physical position for each different material
  void print();
  void print() const;
};

} // namespace olb

#endif

