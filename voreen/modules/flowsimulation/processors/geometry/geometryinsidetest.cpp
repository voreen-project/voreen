/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "geometryinsidetest.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"


namespace olb {

/*  The code within this very "olb" namespace is adapted from the OpenLB library.
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

namespace util
{
    
template<typename T>
bool nearZero(const T &val, float epsilon = std::numeric_limits<T>::epsilon()) {
    return std::abs(val) < epsilon;
}
    
}

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
    bool testRayIntersect(const tgt::Vector3<T> &pt, const tgt::Vector3<T> &dir, tgt::Vector3<T> &q, T &alpha, const T &rad = T()) {
        float rn = 0.;
        tgt::Vector3<T> testPt = pt + rad * normal;

        for (int i = 0; i < 3; i++) {
            rn += dir[i] * normal[i];
        }

        if (util::nearZero(rn)) {
            return false;
        }
        alpha = d - testPt[0] * normal[0] - testPt[1] * normal[1] - testPt[2] * normal[2];
        //  alpha -= testPt[i] * normal[i];
        alpha /= rn;

        if (alpha < -std::numeric_limits<T>::epsilon()) {
            return false;
        }
        for (int i = 0; i < 3; i++) {
            q[i] = testPt[i] + alpha * dir[i];
        }
        T beta = kBeta;
        for (int i = 0; i < 3; i++) {
            beta += uBeta[i] * q[i];
        }

        // Schnittpunkt q in der Ebene?
        if (beta < -std::numeric_limits<T>::epsilon()) {
            return false;
        }
        T gamma = kGamma;
        for (int i = 0; i < 3; i++) {
            gamma += uGamma[i] * q[i];
        }
        if (gamma < -std::numeric_limits<T>::epsilon()) {
            return false;
        }
        if (1. - beta - gamma < -std::numeric_limits<T>::epsilon()) {
            return false;
        }

        return true;
    }

    tgt::Vector3<T> closestPtPointTriangle(const tgt::Vector3<T> &pt) const {
        const T nEps = -std::numeric_limits<T>::epsilon();
        const T Eps = std::numeric_limits<T>::epsilon();

        tgt::Vector3<T> ab = point[1] - point[0];
        tgt::Vector3<T> ac = point[2] - point[0];
        tgt::Vector3<T> bc = point[2] - point[1];

        T snom = (pt - point[0]) * ab;
        T sdenom = (pt - point[1]) * (point[0] - point[1]);

        T tnom = (pt - point[0]) * ac;
        T tdenom = (pt - point[2]) * (point[0] - point[2]);

        if (snom < nEps && tnom < nEps) {
            return point[0];
        }

        T unom = (pt - point[1]) * bc;
        T udenom = (pt - point[2]) * (point[1] - point[2]);

        if (sdenom < nEps && unom < nEps) {
            return point[1];
        }
        if (tdenom < nEps && udenom < nEps) {
            return point[2];
        }

        T vc = normal * tgt::cross(point[0] - pt, point[1] - pt);

        if (vc < nEps && snom > Eps && sdenom > Eps) {
            return point[0] + snom / (snom + sdenom) * ab;
        }

        T va = normal * tgt::cross(point[1] - pt, point[2] - pt);

        if (va < nEps && unom > Eps && udenom > Eps) {
            return point[1] + unom / (unom + udenom) * bc;
        }

        T vb = normal * tgt::cross(point[2] - pt, point[0] - pt);

        if (vb < nEps && tnom > Eps && tdenom > Eps) {
            return point[0] + tnom / (tnom + tdenom) * ac;
        }

        T u = va / (va + vb + vc);
        T v = vb / (va + vb + vc);
        T w = 1. - u - v;

        return u * point[0] + v * point[1] + w * point[2];
    }

    /// A triangle contains 3 Points
    std::vector<tgt::Vector3<T>> point;

    /// normal of triangle
    tgt::Vector3<T> normal;

    /// variables explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
    /// precomputed for speedup
    tgt::Vector3<T> uBeta, uGamma;
    T d, kBeta, kGamma;

public:
    /// Constructor constructs
    STLtriangle(): point(3, tgt::Vector3<T>()), normal(T()), uBeta(T()), uGamma(T()), d(T()), kBeta(T()), kGamma(T()) {
    };
    /// CopyConstructor copies
    STLtriangle(STLtriangle<T> const &tri): point(tri.point), normal(tri.normal), uBeta(tri.uBeta),
                                            uGamma(tri.uGamma), d(tri.d), kBeta(tri.kBeta), kGamma(tri.kGamma) {
    };
    /// Operator= equals
    STLtriangle<T> &operator=(STLtriangle<T> const &tri) {
        point = tri.point;
        normal = tri.normal;
        uBeta = tri.uBeta;
        uGamma = tri.uGamma;
        d = tri.d;
        kBeta = tri.kBeta;
        kGamma = tri.kGamma;
        return *this;
    };

    ~STLtriangle() {
    };

    /// Initializes triangle and precomputes member variables.
    void init() {
        tgt::Vector3<T> A = point[0];
        tgt::Vector3<T> B = point[1];
        tgt::Vector3<T> C = point[2];
        tgt::Vector3<T> b, c;
        T bb = 0., bc = 0., cc = 0.;

        for (int i = 0; i < 3; i++) {
            b[i] = B[i] - A[i];
            c[i] = C[i] - A[i];
            bb += b[i] * b[i];
            bc += b[i] * c[i];
            cc += c[i] * c[i];
        }

        normal[0] = b[1] * c[2] - b[2] * c[1];
        normal[1] = b[2] * c[0] - b[0] * c[2];
        normal[2] = b[0] * c[1] - b[1] * c[0];

        normal = tgt::normalize(normal);

        T D = 1.0 / (cc * bb - bc * bc);
        T bbD = bb * D;
        T bcD = bc * D;
        T ccD = cc * D;

        kBeta = 0.;
        kGamma = 0.;
        d = 0.;

        for (int i = 0; i < 3; i++) {
            uBeta[i] = b[i] * ccD - c[i] * bcD;
            uGamma[i] = c[i] * bbD - b[i] * bcD;
            kBeta -= A[i] * uBeta[i];
            kGamma -= A[i] * uGamma[i];
            d += A[i] * normal[i];
        }
    }

    /// Return write access to normal
    tgt::Vector3<T> &getNormal() {
        return normal;
    }

    /// Return read access to normal
    const tgt::Vector3<T> &getNormal() const {
        return normal;
    }

    /// Returns center
    tgt::Vector3<T> getCenter() {
        return (point[0] + point[1] + point[2]) / 3.0f;
    }

    /// Returns Pt0-Pt1
    std::vector<T> getE0() {
        return point[0] - point[1];
    }

    /// Returns Pt0-Pt2
    std::vector<T> getE1() {
        return point[0] - point[2];
    }

    /// Check whether a point is inside a triangle
    bool isPointInside(const tgt::Vector3<T> &pt) const {
        // tests with T=double and T=float show that the epsilon must be increased
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(10);

        const T beta = pt * uBeta + kBeta;
        const T gamma = pt * uGamma + kGamma;

        // check if approximately equal
        if (util::nearZero(norm(pt - (point[0] + beta * (point[1] - point[0]) + gamma * (point[2] - point[0]))),
                     epsilon)) {
            const T alpha = T(1) - beta - gamma;
            return (beta >= T(0) || util::nearZero(beta, epsilon))
                   && (gamma >= T(0) || util::nearZero(gamma, epsilon))
                   && (alpha >= T(0) || util::nearZero(alpha, epsilon));
        }
        return false;
    }
};

template<typename T>
class STLmesh {
    /// Computes distance squared betwenn p1 and p2
    T distPoints(tgt::Vector3<T> &p1, tgt::Vector3<T> &p2) {
        return tgt::distanceSq(p1, p2);
    }

    /// Vector of Triangles
    std::vector<STLtriangle<T> > _triangles;
    /// Min and Max points of axis aligned bounding box coordinate in SI units
    tgt::Vector3<T> _min, _max;
    /// largest squared length of edge of all triangles
    T _maxDist2;

public:
    /**
     * Constructs a new STLmesh from a file
     * \param Filename - Filename
     * \param stlSize - Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
     */
    STLmesh(const std::vector<std::vector<T> > meshPoints, T stlSize = 1.)
        : _min(T()),
          _max(T()),
          _maxDist2(0) {
        _triangles.reserve(meshPoints.size() / 3);
        for (size_t i = 0; i < meshPoints.size() / 3; i++) {
            STLtriangle<T> tri;
            tri.point[0][0] = meshPoints[i * 3 + 0][0];
            tri.point[0][1] = meshPoints[i * 3 + 0][1];
            tri.point[0][2] = meshPoints[i * 3 + 0][2];

            tri.point[1][0] = meshPoints[i * 3 + 1][0];
            tri.point[1][1] = meshPoints[i * 3 + 1][1];
            tri.point[1][2] = meshPoints[i * 3 + 1][2];

            tri.point[2][0] = meshPoints[i * 3 + 2][0];
            tri.point[2][1] = meshPoints[i * 3 + 2][1];
            tri.point[2][2] = meshPoints[i * 3 + 2][2];
            for (int k = 0; k < 3; k++) {
                tri.point[0][k] *= stlSize;
                tri.point[1][k] *= stlSize;
                tri.point[2][k] *= stlSize;
            }
            if (i == 0) {
                _min *= T();
                _max *= T();

                _min[0] = tri.point[0][0];
                _min[1] = tri.point[0][1];
                _min[2] = tri.point[0][2];

                _max[0] = tri.point[0][0];
                _max[1] = tri.point[0][1];
                _max[2] = tri.point[0][2];

                _min[0] = std::min(_min[0], (T) tri.point[1][0]);
                _min[1] = std::min(_min[1], (T) tri.point[1][1]);
                _min[2] = std::min(_min[2], (T) tri.point[1][2]);

                _max[0] = std::max(_max[0], (T) tri.point[1][0]);
                _max[1] = std::max(_max[1], (T) tri.point[1][1]);
                _max[2] = std::max(_max[2], (T) tri.point[1][2]);

                _min[0] = std::min(_min[0], (T) tri.point[2][0]);
                _min[1] = std::min(_min[1], (T) tri.point[2][1]);
                _min[2] = std::min(_min[2], (T) tri.point[2][2]);

                _max[0] = std::max(_max[0], (T) tri.point[2][0]);
                _max[1] = std::max(_max[1], (T) tri.point[2][1]);
                _max[2] = std::max(_max[2], (T) tri.point[2][2]);
            } else {
                _min[0] = std::min(_min[0], (T) tri.point[0][0]);
                _min[1] = std::min(_min[1], (T) tri.point[0][1]);
                _min[2] = std::min(_min[2], (T) tri.point[0][2]);

                _max[0] = std::max(_max[0], (T) tri.point[0][0]);
                _max[1] = std::max(_max[1], (T) tri.point[0][1]);
                _max[2] = std::max(_max[2], (T) tri.point[0][2]);

                _min[0] = std::min(_min[0], (T) tri.point[1][0]);
                _min[1] = std::min(_min[1], (T) tri.point[1][1]);
                _min[2] = std::min(_min[2], (T) tri.point[1][2]);

                _max[0] = std::max(_max[0], (T) tri.point[1][0]);
                _max[1] = std::max(_max[1], (T) tri.point[1][1]);
                _max[2] = std::max(_max[2], (T) tri.point[1][2]);

                _min[0] = std::min(_min[0], (T) tri.point[2][0]);
                _min[1] = std::min(_min[1], (T) tri.point[2][1]);
                _min[2] = std::min(_min[2], (T) tri.point[2][2]);

                _max[0] = std::max(_max[0], (T) tri.point[2][0]);
                _max[1] = std::max(_max[1], (T) tri.point[2][1]);
                _max[2] = std::max(_max[2], (T) tri.point[2][2]);
            }

            tri.init();
            _triangles.push_back(tri);

            _maxDist2 = std::max(distPoints(tri.point[0], tri.point[1]),
                                 _maxDist2);
            _maxDist2 = std::max(distPoints(tri.point[2], tri.point[1]),
                                 _maxDist2);
            _maxDist2 = std::max(distPoints(tri.point[0], tri.point[2]),
                                 _maxDist2);
        }
    }

    /// Returns reference to a triangle
    STLtriangle<T> &getTri(unsigned int i) {
        return _triangles[i];
    }

    /// Returns reference to all triangles
    std::vector<STLtriangle<T> > &getTriangles() {
        return _triangles;
    }

    /// Returns number of triangles
    unsigned int triangleSize() const {
        return _triangles.size();
    }

    /// Returns _min
    tgt::Vector3<T> &getMin() {
        return _min;
    };
    /// Returns _max
    tgt::Vector3<T> &getMax() {
        return _max;
    };
    /// Returns maxDist squared
    float maxDist2() const {
        return _maxDist2;
    }

    /// Compute intersection between Ray and set of triangles; returns true if intersection is found
    bool testRayIntersect(const std::set<unsigned int> &tris, const tgt::Vector3<T> &pt, const tgt::Vector3<T> &dir,
                          tgt::Vector3<T> &q, T &alpha) {
        std::set<unsigned int>::iterator it = tris.begin();
        for (; it != tris.end(); ++it) {
            if (_triangles[*it].testRayIntersect(pt, dir, q, alpha) && alpha < 1) {
                return true;
            }
        }
        return false;
    }
};

template<typename T>
class Octree {
public:
    /*
     * Constructs Octree containing triangles of an STLmesh.
     * \param center Centerpoint
     * \param rad Radius
     * \param mesh STLmesh
     * \param maxDepth Maximal depth of tree
     * \param overlap Triangles within rad+overlap are added to this Octree
     */
    Octree(tgt::Vector3<T> center, T rad, STLmesh<T> *mesh, short maxDepth, T overlap = 0., Octree<T> *parent = nullptr)
        : _center(center), _radius(rad), _mesh(mesh), _maxDepth(maxDepth), _isLeaf(false), _boundaryNode(false),
          _inside(false), _parent(parent), _child(nullptr) {
        findTriangles(overlap);
        if (_triangles.size() > 0 && 0 < _maxDepth) {
            _child = new Octree<T> *[8];

            tgt::Vector3<T> tmpCenter = _center;
            T tmpRad = _radius / 2.;
            tmpCenter[0] = _center[0] - tmpRad;
            tmpCenter[1] = _center[1] - tmpRad;
            tmpCenter[2] = _center[2] + tmpRad;
            _child[0] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] + tmpRad;
            tmpCenter[1] = _center[1] - tmpRad;
            tmpCenter[2] = _center[2] + tmpRad;
            _child[1] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] - tmpRad;
            tmpCenter[1] = _center[1] - tmpRad;
            tmpCenter[2] = _center[2] - tmpRad;
            _child[2] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] + tmpRad;
            tmpCenter[1] = _center[1] - tmpRad;
            tmpCenter[2] = _center[2] - tmpRad;
            _child[3] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] - tmpRad;
            tmpCenter[1] = _center[1] + tmpRad;
            tmpCenter[2] = _center[2] + tmpRad;
            _child[4] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] + tmpRad;
            tmpCenter[1] = _center[1] + tmpRad;
            tmpCenter[2] = _center[2] + tmpRad;
            _child[5] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] - tmpRad;
            tmpCenter[1] = _center[1] + tmpRad;
            tmpCenter[2] = _center[2] - tmpRad;
            _child[6] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

            tmpCenter[0] = _center[0] + tmpRad;
            tmpCenter[1] = _center[1] + tmpRad;
            tmpCenter[2] = _center[2] - tmpRad;
            _child[7] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);
        } else {
            _isLeaf = true;
            if (_triangles.size() > 0) {
                _boundaryNode = true;
            }
        }
    }

    /// Destructor destructs
    ~Octree() {
        if (_maxDepth != 0 && !_isLeaf) {
            for (int i = 0; i < 8; i++) {
                delete _child[i];
            }
            delete[] _child;
        }
    }

    /// Find the node containing the first param with remaining maxDepth
    Octree<T> *find(const tgt::Vector3<T> &pt, const int &maxDepth = 0) {
        if (_isLeaf || maxDepth == _maxDepth) {
            if (std::abs(_center[0] - pt[0]) < _radius + std::numeric_limits<T>::epsilon() &&
                std::abs(_center[1] - pt[1]) < _radius + std::numeric_limits<T>::epsilon() &&
                std::abs(_center[2] - pt[2]) < _radius + std::numeric_limits<T>::epsilon()) {
                return this;
            } else {
                //throw std::runtime_error("[Octree->find] Point outside of geometry.");
                return nullptr;
            }
        } else {
            if (pt[0] < _center[0]) {
                if (pt[1] < _center[1]) {
                    if (pt[2] < _center[2]) {
                        return _child[2]->find(pt, maxDepth);
                    } else {
                        return _child[0]->find(pt, maxDepth);
                    }
                } else {
                    if (pt[2] < _center[2]) {
                        return _child[6]->find(pt, maxDepth);
                    } else {
                        return _child[4]->find(pt, maxDepth);
                    }
                }
            } else {
                if (pt[1] < _center[1]) {
                    if (pt[2] < _center[2]) {
                        return _child[3]->find(pt, maxDepth);
                    } else {
                        return _child[1]->find(pt, maxDepth);
                    }
                } else {
                    if (pt[2] < _center[2]) {
                        return _child[7]->find(pt, maxDepth);
                    } else {
                        return _child[5]->find(pt, maxDepth);
                    }
                }
            }
        }
    }

    /// Test intersection of ray with all triangles in Octree
    /// returns number of intersections
    int testIntersection(const tgt::Vector3<T> &pt, const tgt::Vector3<T> &dir) {
        int intersections = 0;
        tgt::Vector3<T> q;
        std::vector<tgt::Vector3<T>> qs;
        T a;

        for (unsigned k = 0; k < _triangles.size(); ++k) {
            if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, dir, q, a)) {
                if (std::abs(_center[0] - q[0]) <= _radius + std::numeric_limits<T>::epsilon() + 1 / 1000. * _radius
                    && std::abs(_center[1] - q[1]) <= _radius + std::numeric_limits<T>::epsilon() + 1 / 1000. *
                    _radius && std::abs(_center[2] - q[2]) <= _radius + std::numeric_limits<T>::epsilon() + 1 /
                    1000. * _radius) {
                    bool newpoint = true;
                    for (unsigned i = 0; i < qs.size(); i++) {
                        newpoint = (!util::nearZero(q[0] - qs[i][0]) || !util::nearZero(q[1] - qs[i][1]) || !util::nearZero(
                                        q[2] - qs[i][2]));
                    }
                    if (newpoint) {
                        qs.push_back(q);
                        intersections++;
                    }
                }
            }
        }

        return intersections;
    }

    /// Test intersection of ray with all triangles in Octree
    /// q contains point of closest intersection to pt in direction direction
    bool closestIntersection(const tgt::Vector3<T> &pt, const tgt::Vector3<T> &direction, tgt::Vector3<T> &q, T &a) {
        STLtriangle<T> tri;
        return closestIntersection(pt, direction, q, a, tri, 0.);
    }

    /// Test intersection of ray with all triangles in Octree
    /// q contains point of closest intersection to pt in direction direction
    /// tri contains triangle with closest intersection
    bool closestIntersection(const tgt::Vector3<T> &pt, const tgt::Vector3<T> &direction, tgt::Vector3<T> &q, T &a,
                             STLtriangle<T> &tri, const T &rad = T()) {
        a = std::numeric_limits<T>::infinity();
        T alpha = T();
        tgt::Vector3<T> qtmp;
        bool found = false;

        for (unsigned int k = 0; k < _triangles.size(); ++k) {
            if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, direction, qtmp, alpha, rad)) {
                if (alpha < a) {
                    a = alpha;
                    q = qtmp;
                    found = true;
                    tri = _mesh->getTri(_triangles[k]);
                }
            }
        }
        return found;
    }

    /// Test intersection of sphere moving along ray with radius rad
    /// q contains point of closest intersection to pt in direction direction
    /// tri contains triangle with closest intersection
    bool closestIntersectionSphere(const tgt::Vector3<T> &pt, const T &rad, const tgt::Vector3<T> &direction, tgt::Vector3<T> &q,
                                   T &a, STLtriangle<T> &tri) {
        a = std::numeric_limits<T>::infinity();
        T alpha = T();
        std::vector<T> qtmp(3, T());
        bool found = false;
        for (unsigned int k = 0; k < _triangles.size(); ++k) {
            if (_mesh->getTri(_triangles[k]).testMovingSphereIntersect(pt, rad, direction, qtmp, alpha)) {
                if (alpha < a) {
                    a = alpha;
                    q = qtmp;
                    found = true;
                    tri = _mesh->getTri(_triangles[k]);
                }
            }
        }
        return found;
    }

    /// It's complicated. Computes intersections of a ray with triangles inside this Octree. Sets _inside depending on value of rayInside and changes rayInside depending on the number of intersections. Also takes into account if the intersections happen before or after the center.
    void checkRay(const tgt::Vector3<T> &pt, const tgt::Vector3<T> &dir, unsigned short &rayInside) {
        unsigned short left = 0, right = 0;
        tgt::Vector3<T> dirNormed(dir);
        dirNormed = normalize(dirNormed);
        dirNormed *= _radius * 2.;
        tgt::Vector3<T> q;
        std::vector<tgt::Vector3<T>> qs;
        T a = 1.;

        for (unsigned int k = 0; k < _triangles.size(); ++k) {
            if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, dirNormed, q, a, 0.) && a < 1.) {
                bool newpoint = true;
                for (unsigned int i = 0; i < qs.size(); i++) {
                    newpoint &= (!util::nearZero(q[0] - qs[i][0]) || !util::nearZero(q[1] - qs[i][1]) || !util::nearZero(
                                     q[2] - qs[i][2]));
                }
                if (newpoint) {
                    qs.push_back(q);
                    if (a < .5) {
                        left++;
                    } else {
                        right++;
                    }
                }
            }
        }
        rayInside += left;
        rayInside %= 2;
        setInside(rayInside);
        rayInside += right;
        rayInside %= 2;
    }

    /// Computes intersection of ray with Octree boundaries
    void intersectRayNode(const tgt::Vector3<T> &pt, const tgt::Vector3<T> &dir, tgt::Vector3<T> &s) {
        T t, d;
        s *= T();
        //Plane Normals outside

        if (dir[0] > 0.) {
            // n = {1, 0, 0}
            d = _center[0] + _radius;
            t = (d - pt[0]) / dir[0];
            s = pt + t * dir;
            if (std::abs(s[1] - _center[1]) < _radius && std::abs(s[2] - _center[2]) < _radius) {
                return;
            }
        } else if (dir[0] < 0.) {
            // n = {-1, 0, 0}
            d = _center[0] - _radius;
            t = (d - pt[0]) / dir[0];
            s = pt + t * dir;
            if (std::abs(s[1] - _center[1]) < _radius && std::abs(s[2] - _center[2]) < _radius) {
                return;
            }
        }

        if (dir[1] > 0.) {
            d = _center[1] + _radius;
            t = (d - pt[1]) / dir[1];
            s = pt + t * dir;
            if (std::abs(s[0] - _center[0]) < _radius && std::abs(s[2] - _center[2]) < _radius) {
                return;
            }
        } else if (dir[1] < 0.) {
            // n = {0, 0, -1}
            d = _center[1] - _radius;
            t = (d - pt[1]) / dir[1];
            s = pt + t * dir;
            if (std::abs(s[0] - _center[0]) < _radius && std::abs(s[2] - _center[2]) < _radius) {
                return;
            }
        }

        if (dir[2] > 0.) {
            // n = {0, 0, 1}
            d = _center[2] + _radius;
            t = (d - pt[2]) / dir[2];
            s = pt + t * dir;
            if (std::abs(s[0] - _center[0]) < _radius && std::abs(s[1] - _center[1]) < _radius) {
                return;
            }
        } else if (dir[2] < 0.) {
            // n = {0, 0, -1}
            d = _center[2] - _radius;
            t = (d - pt[2]) / dir[2];
            s = pt + t * dir;
            if (std::abs(s[0] - _center[0]) < _radius && std::abs(s[1] - _center[1]) < _radius) {
                return;
            }
        }
    }

    /// Computes all centerpoints of Octree
    void getCenterpoints(std::vector<std::vector<T> > &pts) {
        if (_isLeaf) {
            pts.push_back(_center);
        } else {
            for (int i = 0; i < 8; i++) {
                _child[i]->getCenterpoints(pts);
            }
        }
    }

    /// Collectes all leafs
    void getLeafs(std::vector<Octree<T> *> &pts) {
        if (_isLeaf) {
            pts.push_back(this);
        } else {
            for (int i = 0; i < 8; i++) {
                _child[i]->getLeafs(pts);
            }
        }
    }

    /// Return status of _isLeaf;
    bool isLeaf() {
        return _isLeaf;
    }

    /// Sets Inside
    void setInside(bool ins) {
        _inside = ins;
    };
    /// Gets Inside
    bool getInside() {
        return _inside;
    };
    /// Gets _boundarNode
    bool getBoundaryNode() {
        return _boundaryNode;
    };
    /// Gets Maxdepth
    int getMaxdepth() const {
        return _maxDepth;
    };
    /// Gets numbers of triangles contained by this Octree
    const std::vector<unsigned int> &getTriangles() const {
        return _triangles;
    };
    /// Gets centerpoint
    const tgt::Vector3<T> &getCenter() const {
        return _center;
    };
    /// Gets radius
    const T getRadius() const {
        return _radius;
    };

    /// Returns set of indices of all triangles in nodes containing a line.
    void trianglesOnLine(const tgt::Vector3<T> &pt1, const tgt::Vector3<T> &pt2, std::set<unsigned int> &tris) {
        tris.clear();
        std::vector<T> line = pt2 - pt1;
        std::vector<T> s = pt1;
        T lineNorm2 = line[0] * line[0] + line[1] * line[1] + line[2] * line[2];
        T dist2 = T();
        Octree<T> *node = nullptr;
        int it = 0;
        while (dist2 < lineNorm2 && it < 50) {
            node = find(s);
            tris.insert(node->_triangles.begin(), node->_triangles.end());
            node->intersectRayNode(s, line, s);
            for (int i = 0; i < 3; i++) {
                s[i] = s[i] + line[i] * _radius * 0.001 /* *node->getRadius()*/;
            }
            it++;
            dist2 = (pt1[0] - s[0]) * (pt1[0] - s[0]) + (pt1[1] - s[1]) * (pt1[1] - s[1]) + (pt1[2] - s[2]) * (
                        pt1[2] - s[2]);
        }
    }

    /// Returns reference to _mesh
    STLmesh<T> *getMesh() {
        return _mesh;
    }

protected:
    ///_vector _triangles contains number of triangles
    std::vector<unsigned int> _triangles;
    tgt::Vector3<T> _center;
    T _radius;
    STLmesh<T> *_mesh;
    short _maxDepth;
    bool _isLeaf;
    bool _boundaryNode;
    bool _inside;
    Octree<T> *_parent;
    Octree<T> **_child;

    void findTriangles(T overlap = 0.) {
        if (_parent == nullptr) {
            _triangles.reserve(_mesh->triangleSize());
            for (unsigned int i = 0; i < _mesh->triangleSize(); ++i) {
                if (AABBTri(_mesh->getTri(i))) {
                    _triangles.push_back(i);
                }
            }
        } else {
            std::vector<unsigned int>::iterator it;
            for (it = _parent->_triangles.begin(); it != _parent->_triangles.end(); ++it) {
                if (AABBTri(_mesh->getTri(*it), overlap)) {
                    _triangles.push_back(*it);
                }
            }
        }
    }

    bool AABBTri(const STLtriangle<T> &tri, T overlap = 0.) {
        std::vector<T> v0(3, T()), v1(3, T()), v2(3, T()), f0(3, T()), f1(3, T()), f2(3, T()), e(3, T());

        /* Test intersection cuboids - triangle
        * Intersection test after Christer Ericson - Real time Collision Detection p.
        * TestTriangleAABB p.171 */
        tgt::Vector3<T> c(_center);
        T eps = std::numeric_limits<T>::epsilon();

        for (int j = 0; j < 3; j++) {
            v0[j] = tri.point[0][j] - _center[j];
            v1[j] = tri.point[1][j] - _center[j];
            v2[j] = tri.point[2][j] - _center[j];
            e[j] = _radius * 1.01 + overlap; // + std::numeric_limits<T>::epsilon(); // *1.01;
        }
        for (int j = 0; j < 3; j++) {
            f0[j] = v1[j] - v0[j];
            f1[j] = v2[j] - v1[j];
            f2[j] = v0[j] - v2[j];
        }
        T p0 = T(), p1 = T(), r = T();
        //test a00
        p0 = v0[2] * v1[1] - v0[1] * v1[2];
        p1 = v2[2] * v1[1] - v2[2] * v0[1] + v0[2] * v2[1] - v1[2] * v2[1];
        r = e[1] * std::abs(f0[2]) + e[2] * std::abs(f0[1]);
        T mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a01
        p0 = v0[1] * v1[2] - v0[1] * v2[2] - v0[2] * v1[1] + v0[2] * v2[1];
        p1 = -v1[1] * v2[2] + v1[2] * v2[1];
        r = e[1] * std::abs(f1[2]) + e[2] * std::abs(f1[1]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a02
        p0 = v0[1] * v2[2] - v0[2] * v2[1];
        p1 = v0[1] * v1[2] - v0[2] * v1[1] + v1[1] * v2[2] - v1[2] * v2[1];
        r = e[1] * std::abs(f2[2]) + e[2] * std::abs(f2[1]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a10
        p0 = v0[0] * v1[2] - v0[2] * v1[0];
        p1 = v0[0] * v2[2] - v0[2] * v2[0] - v1[0] * v2[2] + v1[2] * v2[0];
        r = e[0] * std::abs(f0[2]) + e[2] * std::abs(f0[0]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a11
        p0 = -v0[0] * v1[2] + v0[0] * v2[2] + v0[2] * v1[0] - v0[2] * v2[0];
        p1 = v1[0] * v2[2] - v1[2] * v2[0];
        r = (T) (e[0] * std::abs(f1[2]) + e[2] * std::abs(f1[0]));
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a12
        p0 = -v0[0] * v2[2] + v0[2] * v2[0];
        p1 = -v0[0] * v1[2] + v0[2] * v1[0] - v1[0] * v2[2] + v1[2] * v2[0];
        r = e[0] * std::abs(f2[2]) + e[2] * std::abs(f2[0]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a20
        p0 = -v0[0] * v1[1] + v0[1] * v1[0];
        p1 = -v0[0] * v2[1] + v0[1] * v2[0] + v1[0] * v2[1] - v1[1] * v2[0];
        r = e[0] * std::abs(f0[1]) + e[1] * std::abs(f0[0]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a21
        p0 = v0[0] * v1[1] - v0[0] * v2[1] - v0[1] * v1[0] + v0[1] * v2[0];
        p1 = -v1[0] * v2[1] + v1[1] * v2[0];
        r = e[0] * std::abs(f1[1]) + e[1] * std::abs(f1[0]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        // test a22
        p0 = v0[0] * v2[1] - v0[1] * v2[0];
        p1 = v0[0] * v1[1] - v0[1] * v1[0] + v1[0] * v2[1] - v1[1] * v2[0];
        r = e[0] * std::abs(f2[1]) + e[1] * std::abs(f2[0]);
        mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
        if (mmm > r + eps) {
            return false;
        }

        if (std::max(std::max(v0[0], v1[0]), v2[0]) < -e[0] || std::min(std::min(v0[0], v1[0]), v2[0]) > e[0]) {
            return false;
        }
        if (std::max(std::max(v0[1], v1[1]), v2[1]) < -e[1] || std::min(std::min(v0[1], v1[1]), v2[1]) > e[1]) {
            return false;
        }
        if (std::max(std::max(v0[2], v1[2]), v2[2]) < -e[2] || std::min(std::min(v0[2], v1[2]), v2[2]) > e[2]) {
            return false;
        }

        /* Test intersection cuboids - triangle plane*/
        r = e[0] * std::abs(tri.normal[0]) + e[1] * std::abs(tri.normal[1]) + e[2] * std::abs(tri.normal[2]);
        T s = tri.normal[0] * c[0] + tri.normal[1] * c[1] + tri.normal[2] * c[2] - tri.d;
        return (std::abs(s) <= r);
    }
};

template<typename T>
void indicate(Octree<T>& tree, T spacing) {

    auto* mesh = tree.getMesh();

    std::vector<Octree<T>*> leafs;
    tree.getLeafs(leafs);
    auto it = leafs.begin();
    tgt::Vector3<T> dir, pt, s;

    int intersections = 0;
    int inside = 0;
    Octree<T>* node = nullptr;
    float step = 1. / 1000. * spacing;
    for (; it != leafs.end(); ++it) {
        inside = 0;

        pt = (*it)->getCenter();
        intersections = 0;
        s = pt;  // + step;

        /// X+ dir
        dir[0] = 1;
        dir[1] = 0;
        dir[2] = 0;
        while (s[0] < mesh->getMax()[0] + std::numeric_limits<float>::epsilon()) {
            node = tree.find(s, (*it)->getMaxdepth());
            intersections += node->testIntersection(pt, dir);
            node->intersectRayNode(pt, dir, s);
            s = s + step * dir;
        }
        inside += (intersections % 2);

        /// Y+ Test
        intersections = 0;
        s = pt;  // + step;
        dir[0] = 0;
        dir[1] = 1;
        dir[2] = 0;
        while (s[1] < mesh->getMax()[1] + std::numeric_limits<float>::epsilon()) {
            node = tree.find(s, (*it)->getMaxdepth());
            intersections += node->testIntersection(pt, dir);
            node->intersectRayNode(pt, dir, s);
            s = s + step * dir;
        }
        inside += (intersections % 2);

        /// Z+ Test
        intersections = 0;
        s = pt;  // + step;
        dir[0] = 0;
        dir[1] = 0;
        dir[2] = 1;
        while (s[2] < mesh->getMax()[2] + std::numeric_limits<float>::epsilon()) {
            node = tree.find(s, (*it)->getMaxdepth());
            intersections += node->testIntersection(pt, dir);
            node->intersectRayNode(pt, dir, s);
            s = s + step * dir;
        }
        inside += (intersections % 2);
        (*it)->setInside(inside > 1);
    }
}

}




namespace {

template<typename T>
std::vector<std::vector<float>> getMeshPoints(const T* geometry) {

    tgtAssert(geometry, "geometry is null");
    tgtAssert(geometry->getPrimitiveType() == GL_TRIANGLES, "geometry has unsupported primitive type");

    std::vector<std::vector<float>> meshPoints;

    // Transform all vertex position to world space.
    auto addPoint = [&](tgt::vec3 point) {
        point = geometry->getTransformationMatrix() * point;
        meshPoints.push_back({point.x, point.y, point.z});
    };

    auto vertices = geometry->getVertices();

    if (geometry->usesIndexedDrawing()) {
        for (auto& index : geometry->getIndices()) {
            addPoint(vertices[index].pos_);
        }
    }
    else {
        for (auto& vertex: vertices) {
            addPoint(vertex.pos_);
        }
    }

    return meshPoints;
}

}





namespace voreen {

const std::string GeometryInsideTest::loggerCat_("voreen.flowsimulation.GeometryInsideTest");

GeometryInsideTest::GeometryInsideTest()
    : Processor()
    , inport_(Port::INPORT, "geometryinsidetest.inport", "")
    , outport_(Port::OUTPORT, "geometryinsidetest.outport", "ID Volume Output", false)
    , dimensions_("dimensions", "Dimensions", 256, 32, 1024)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(dimensions_);
    dimensions_.setTracking(false);
}

Processor* GeometryInsideTest::create() const {
    return new GeometryInsideTest();
}

bool GeometryInsideTest::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    if(!dynamic_cast<const GlMeshGeometryBase*>(inport_.getData())) {
        setNotReadyErrorMessage("Invalid input");
        return false;
    }

    return true;
}

void GeometryInsideTest::process() {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Test input data.
    const GlMeshGeometryBase* inputGeometry = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData());
    tgtAssert(inputGeometry, "Invalid input");

    if(inputGeometry->getNumVertices() == 0) {
        LERROR("Geometry is empty!");
        outport_.setData(nullptr);
        return;
    }

    if (inputGeometry->getPrimitiveType() != GL_TRIANGLES) {
        LERROR("Only triangles are supported");
        outport_.setData(nullptr);
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gather mesh points
    std::vector<std::vector<float>> meshPoints;

    if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32Simple*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32Color*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32Normal*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32TexCoord*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32ColorNormal*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32NormalTexCoord*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32ColorTexCoord*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else if (auto* geom = dynamic_cast<const GlMeshGeometryUInt32ColorNormalTexCoord*>(inputGeometry)) {
        meshPoints = getMeshPoints(geom);
    }
    else {
        LERROR("Unsupported input geometry type");
        outport_.setData(nullptr);
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Setup mesh and octree
    olb::STLmesh<float> mesh(meshPoints, 1.0f);

    const tgt::vec3 extension = (mesh.getMax() - mesh.getMin());

    const float max = tgt::max(extension);
    const float spacing = max / (dimensions_.get() - 1);

    int j = 0;
    for (; spacing * std::pow(2, j) < max; j++);

    const tgt::vec3 center = (mesh.getMin() + mesh.getMax()) / 2.0f - spacing / 4.0f;
    const float radius = spacing * std::pow(2, j - 1);

    olb::Octree<float> tree(center, radius, &mesh, j);
    indicate(tree, spacing);

    auto inside = [&] (const tgt::vec3& input)
    {
        const int OUTSIDE = 0;
        const int INSIDE = 255;

        float coords = tree.getRadius();
        tgt::vec3 c(tree.getCenter());
        if (c[0] - coords < input[0] && input[0] < c[0] + coords && c[1] - coords < input[1]
            && input[1] < c[1] + coords && c[2] - coords < input[2] && input[2] < c[2] + coords) {
            return tree.find(input)->getInside() ? INSIDE : OUTSIDE;
        }
        return OUTSIDE;
    };


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Write volume.
    tgt::svec3 dim(dimensions_.get());
    VolumeRAM_UInt8* idVolume = new VolumeRAM_UInt8(dim);

    tgt::vec3 offset = mesh.getMin();
    Volume* outputVolume = new Volume(idVolume, tgt::vec3(spacing), offset);

    auto voxelToWorldMatrix = outputVolume->getVoxelToWorldMatrix();

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for(int z=0; z<dim.z; z++) {
        for(int y=0; y<dim.y; y++) {
            for(int x=0; x<dim.x; x++) {
                tgt::vec3 pos = voxelToWorldMatrix * tgt::vec3(x, y, z);
                idVolume->voxel(x, y, z) = inside(pos);
            }
        }
    }

    outport_.setData(outputVolume);
}

}