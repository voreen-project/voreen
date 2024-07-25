/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause
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
 * Input in STL format -- header file.
 */

#ifndef STL_READER_H
#define STL_READER_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <limits>

#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry3D.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"
#include "octree.h"
#include "core/vector.h"


// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
class Octree;

template<typename T, unsigned D>
struct Vertex {
  /// Constructor constructs
  Vertex() : coords() {};
  Vertex( Vector<T,D> coords_ ) : coords(coords_) {};
  /// Operator= equals
  Vertex<T,D>& operator=(Vertex<T,D> const& rhs)
  {
    coords = rhs.coords;
    return *this;
  };
  /// CopyConstructor copies
  Vertex(Vertex<T,D> const& rhs):coords(rhs.coords) {};

  /// Point coordinates in SI units
  Vector<T,D> coords;

  /// Get dimension 
  int getDim() const
  {
    return D;
  }
  
};

template<typename T>
struct STLtriangle {
  /** Test intersection between ray and triangle
   * \param pt Raypoint
   * \param dir Direction
   * \param q Point of intersection (if intersection occurs)
   * \param alpha Explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
   *            q = pt + alpha * dir
   * \param rad It's complicated. Imagine you have a sphere with radius rad moving a long the ray. Then q becomes the first point of the sphere to touch the triangle.
   */
  bool testRayIntersect(const Vector<T,3>& pt,const Vector<T,3>& dir, Vector<T,3>& q, T& alpha, const T& rad = T(), bool print = false);
  Vector<T,3> closestPtPointTriangle(const Vector<T,3>& pt) const;

  /// A triangle contains 3 Points
  std::vector<Vertex<T,3> > point;

  /// normal of triangle
  Vector<T,3> normal;

  /// variables explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
  /// precomputed for speedup
  Vector<T,3> uBeta, uGamma;
  T d, kBeta, kGamma;

public:
  /// Constructor constructs
  STLtriangle():point(3, Vertex<T,3>()), normal(T()), uBeta(T()), uGamma(T()), d(T()), kBeta(T()), kGamma(T()) {};
  /// CopyConstructor copies
  STLtriangle(STLtriangle<T> const& tri):point(tri.point), normal(tri.normal), uBeta(tri.uBeta), uGamma(tri.uGamma), d(tri.d), kBeta(tri.kBeta), kGamma(tri.kGamma) {};
  /// Operator= equals
  STLtriangle<T>& operator=(STLtriangle<T> const& tri)
  {
    point = tri.point;
    normal = tri.normal;
    uBeta = tri.uBeta;
    uGamma = tri.uGamma;
    d = tri.d;
    kBeta = tri.kBeta;
    kGamma = tri.kGamma;
    return *this;
  };

  ~STLtriangle() {};

  /// Initializes triangle and precomputes member variables.
  void init();
  /// Return write access to normal
  inline Vector<T,3>& getNormal()
  {
    return normal;
  }
  /// Return read access to normal
  inline const Vector<T,3>& getNormal() const
  {
    return normal;
  }
  /// Returns center
  Vector<T,3> getCenter();
  /// Returns Pt0-Pt1
  std::vector<T> getE0();
  /// Returns Pt0-Pt2
  std::vector<T> getE1();
  /// Check whether a point is inside a triangle
  bool isPointInside(const PhysR<T,3>& pt) const;
};

template<typename T>
class STLmesh {
  /// Computes distance squared betwenn p1 and p2
  T distPoints(Vertex<T,3>& p1, Vertex<T,3>& p2);
  /// Filename
  const std::string _fName;
  /// Vector of Triangles
  std::vector<STLtriangle<T> > _triangles;
  /// Min and Max points of axis aligned bounding box coordinate in SI units
  Vector<T,3> _min, _max;
  /// largest squared length of edge of all triangles
  T _maxDist2;
  /// OstreamManager
  mutable OstreamManager clout;

public:
  /**
   * Constructs a new STLmesh from a file
   * \param Filename - Filename
   * \param stlSize - Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   */
  STLmesh(std::string, T stlSize = 1.);

  /**
   * Constructs a new STLmesh from a file
   * \param Filename - Filename
   * \param stlSize - Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   */
  STLmesh(const std::vector<std::vector<T>> meshPoints, T stlSize = 1.);

  /// Returns reference to a triangle
  inline STLtriangle<T>& getTri(unsigned int i)
  {
    return _triangles[i];
  }
  /// Returns reference to all triangles
  inline std::vector<STLtriangle<T> >& getTriangles()
  {
    return _triangles;
  }
  /// Returns number of triangles
  inline unsigned int triangleSize() const
  {
    return _triangles.size();
  }
  /// Returns _min
  inline Vector<T,3>& getMin()
  {
    return _min;
  };
  /// Returns _max
  inline Vector<T,3>& getMax()
  {
    return _max;
  };
  /// Returns maxDist squared
  inline float maxDist2() const
  {
    return _maxDist2;
  }
  /// Prints console output
  void print(bool full = false);
  /// Writes STL mesh in Si units
  void write(std::string fName);
  /// Compute intersection between Ray and set of triangles; returns true if intersection is found
  bool testRayIntersect(const std::set<unsigned int>& tris, const Vector<T,3>& pt,const Vector<T,3>& dir, Vector<T,3>& q, T& alpha);
};


template<typename T>
class STLreader : public IndicatorF3D<T> {
private:
  /*
   *  Old indicate function (slower, more stable)
   *  Define three rays (X-, Y-, Z-direction) for each leaf and count intersections
   *  with STL for each ray. Odd number of intersection means inside (Majority vote).
   */
  void indicate1();
  /*
   *  New indicate function (faster, less stable)
   *  Define ray in Z-direction for each Voxel in XY-layer. Indicate all nodes on the fly.
   */
  void indicate2();
  /*
  *  New indicate function (faster, less stable)
  *  Define ray in X-direction for each Voxel in YZ-layer. Indicate all nodes on the fly.
  */
  void indicate2_Xray();
  /*
  *  New indicate function (faster, less stable)
  *  Define ray in Y-direction for each Voxel in XZ-layer. Indicate all nodes on the fly.
  */
  void indicate2_Yray();
  /*
   *  Double ray approach: two times (X-, Y-, Z-direction) for each leaf.
   *  Could be use to deal with double layer triangles and face intersections.
   */

  void indicate3();

  /// Finds normal for points on the surface (do not use for points that aren't on the surface!)
  Vector<T,3> findNormalOnSurface(const PhysR<T,3>& pt);

  /// Finds surface normal
  Vector<T,3> evalSurfaceNormal(const Vector<T,3>& origin);

  /// Size of the smallest voxel
  T _voxelSize;
  /// Factor to get Si unit (m), i.e. "0.001" means mm
  T _stlSize;
  /// Overlap increases Octree radius by _overlap
  T _overlap;
  /// Pointer to tree
  Octree<T>* _tree;
  /// The filename
  const std::string _fName;
  /// The mesh
  STLmesh<T> _mesh;
  /// Variable for output
  bool _verbose;
  /// The OstreamManager
  mutable OstreamManager clout;

public:
  /**
   * Constructs a new STLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */

  STLreader(const std::string fName, T voxelSize, T stlSize=1, int method=2,
            bool verbose = false, T overlap=0., T max=0.);
  /**
   * Constructs a new STLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */
  STLreader(const std::vector<std::vector<T>> meshPoints, T voxelSize, T stlSize=1, int method=2,
            bool verbose = false, T overlap=0., T max=0.);

  ~STLreader() override;
  /// Returns whether node is inside or not.
  bool operator() (bool output[], const T input[]) override;

  /// Computes distance to closest triangle intersection
  bool distance(T& distance,const Vector<T,3>& origin, const Vector<T,3>& direction, int iC=-1) override;

  /// Computes signed distance to closest triangle in direction of the surface normal
  T signedDistance(const Vector<T,3>& input) override;

  /// Finds and returns normal of the closest surface (triangle)
  Vector<T,3> surfaceNormal(const Vector<T,3>& pos, const T meshSize=0) override;

  /// Prints console output
  void print();

  /// Writes STL mesh in Si units
  void writeSTL(std::string stlName="");

  /// Writes Octree
  void writeOctree();

  /// Rearranges normals of triangles to point outside of geometry
  void setNormalsOutside();

  /// Every octree leaf intersected by the STL will be part of the inside nodes.
  /// Artificially enlarges all details that would otherwise be cut off by the voxelSize.
  void setBoundaryInsideNodes();

  /// Returns tree
  inline Octree<T>* getTree() const
  {
    return _tree;
  };


  /// Returns mesh
  inline STLmesh<T>& getMesh()
  {
    return _mesh;
  };
};

}  // namespace olb

#endif
