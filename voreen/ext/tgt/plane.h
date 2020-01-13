/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_PLANE_H
#define TGT_PLANE_H

#include <algorithm> // for std::swap
#include <vector>

#include "tgt/assert.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/line.h"

namespace tgt {

/**
 * @brief a plane class
 * The plane is stored in hessian normal form: <br>
 * n*x + d = 0; <br>
 * n is the unit normal of this plane pointing from the origin to this plane <br>
 * (-d) >= 0 is the distance from the origin. <br>
 * <br>
 * Note: You can transform this plane like this: <br>
 * myPlane.toVec4() = myMatrix4 * myPlane.toVec4(); <br>
 * <br>
 * Perhaps you want to do this: <br>
 * mat4 mv = tgt::getModelViewMatrix(); <br>
 * mat4 inv; <br>
 * mv.invert(inv); <br>
 * myPlane.toVec4() = inv * myPlane.toVec4(); <br>
*/
template<class T>
struct Plane {
    typedef T ElemType;

    /// The normal of the plane.
    Vector3<T> n;
    /// The negative distance from the orign.
    T d;

    /// trivial default constructor
    Plane() {}

    /// init with unit normal [a, b, c] and distance _d <br>
    /// BEWARE: It is your job to ensure that [a, b, c] is a unit vector
    Plane(T a, T b, T c, T _d)
        : n(a, b, c)
        , d(_d)
    {}

    /// creates a plane with a unit length normal and a distance <br>
    /// BEWARE: It is your job to ensure that [a, b, c] is a unit vector
    Plane(const Vector3<T>& _n, T _d)
        : n(_n)
        , d(_d)
    {}

    /**
     * Creates a plane defined by three points.
     */
    Plane(const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) {
        n = normalize( cross(v2 - v1, v3 - v1) );
        d = dot(v1, n);
    }

    /// creates a plane with one from another type
    template<class U>
    Plane(const Plane<U>& p)
        : n(p.n)
        , d(T(p.d))
    {}

    const Vector4<T>& toVec4() const {
        return *reinterpret_cast<Vector4<T>*>(const_cast<Plane<T>*>(this)->n.elem);
    }

    Vector4<T>& toVec4() {
        return *reinterpret_cast<Vector4<T>*>(n.elem);
    }

    /**
     * Returns the signed distance to this plane to a given point. <br>
     * If result < 0 -> v lies in the same half space of the plane like the origin. <br>
     * If result > 0 -> v lies in the other half space of the plane like the origin.
     * @return the signed distance from the plane to the point
    */
    T distance(const Vector3<T>& v) const {
        return dot(n, v) - d;
    }

    /**
     * Returns whether this plane intersects with line \p l.
     * When there is an intersection the parameter to compute
     * the point will be returned in \p f.
     * If \p l represents an edge between v1 and v2 that were given to \p l
     * via the constructor an intersection takes place when \p f is in [0, 1]
    */
    bool intersect(const Line3<T>& l, T& f) const;

    /**
     * Returns whether this plane intersects with \p p.
     * When there is an intersection it will be computed in \p l
    */
    bool intersect(const Plane<T>& p, Line3<T>& l) const;

    /**
     * clips a given polygon to this plane
     * @param vin the input polygon
     * @param vout the output polygon
    */
    void clipPolygon(const std::vector< Vector3<T> >& vin, std::vector< Vector3<T> >& vout) const;

    /**
     * Clips a given polygon to this plane.
     * This is a second version with c-arrays. Can be faster, but beware:
     * Buffer overruns are likly!!!
     * @param vin the input polygon
     * @param numIn number of input vertices
     * @param vout the output polygon
     * @param numOut number of output vertices
    */
    void clipPolygon(const Vector3<T>* vin, size_t numIn, Vector3<T>* vout, size_t& numOut) const;

    static void createCubeVertices(const Vector3<T>& llf, const Vector3<T>& urb, Vector3<T> cubeVertices[8]);

    static void createCubePlanes(const Vector3<T> cv[8], Plane<T> cubePlanes[6]);

    Plane<T> transform(tgt::mat4 m) const;

    /**
     * Creates planes in the aabPlanes array from a axis aligned box given via llf and urb. The normal points to inside.
     * @param llf lower left front of the input axis aligned box
     * @param urb upper right back of the input axis aligned box
     * @param cubePlanes the result is written to this array which must be able to hold at least 6 planes
     */
    static void createAABPlanes(const Vector3<T>& llf, const Vector3<T>& urb, Plane<T> aabPlanes[6]) {
        Vector3<T> cv[8];// cube vertices
        createCubeVertices(llf, urb, cv);

        // cube planes
        createCubePlanes(cv, aabPlanes);
    }

    /**
     * This method creates a polygon from this plane.
     * It is ensured that this polygon is huge enough
     * so it may be clipped later against the axis aligned box
     * defined by llf and urb. <br>
     * Note that you can use this method with std::vector, too: <br>
     * std::vector\<vec3\> vout(4);
     * createQuad(llf, urb, &vout[0]);
     * @param llf lower left front of the input axis aligned box
     * @param urb upper right back of the input axis aligned box
     * @param v the result is written to this array which must be able to hold at least 4 vertices
     */
    void createQuad(Vector3<T> llf, Vector3<T> urb, Vector3<T>* v) const;

    /**
     * This routine clips a polygon against 6 planes
     * @param aabPlanes against these 6 planes is clipped
     * @param vin Points to the input vertices. <br>
     *      BEWARE: vin is changed after this call. <br>
     *      AND EVEN WORSE: The size should be enough. Buffer overruns are likely!!!
     * @param numIn number of vertices in vin
     * @param vout Points to the output vertices. This is the result. <br>
     *      BEWARE: The size should be enough. Buffer overruns are likely!!!
     * @param numOut Number of output vertices. This a result, too.
     */
    static void clipAAB(Plane<T> aabPlanes[6], Vector3<T>* vin, size_t numIn, Vector3<T>* vout, size_t& numOut);

    /**
     * This routine clips a polygon against 6 planes
     * @param aabPlanes against these 6 planes is clipped
     * @param vin The input vertices. <br>
     *      BEWARE: vin is changed after this call.
     * @param vout Points to the output vertices. This is the result,
     */
    static void clipAAB(Plane<T> aabPlanes[6], std::vector< Vector3<T> >& vin, std::vector< Vector3<T> >& vout);

    /**
     * This method cilps this plane against an axis aligned box defined by llf and urb.
     * The output is written to vout. <br>
     * BEWARE: It is likely that you can gain performance by precomputing the planes of the AAB once,
     * then calculate a proper polygon with createQuad and finally call a static clipAAB version
     * which takes the planes an the polygon as input.
     * @param llf lower left front of the input axis aligned box
     * @param urb upper right back of the input axis aligned box
     * @param vout the result polygon
     */
    void clipAAB(const Vector3<T>& llf, const Vector3<T>& urb, std::vector< Vector3<T> >& vout) const;
};

template<class T>
bool operator==(const Plane<T>& x, const Plane<T>& y){
    return ( (x.n == y.n) &&  (x.d == y.d) );
}

template<class T>
bool operator!=(const Plane<T>& x, const Plane<T>& y){
    return !(x == y);
}

// non inline implementation

template<class T>
bool Plane<T>::intersect(const Line3<T>& l, T& f) const {
    T f0 = distance(l.getStart());
    T f1 = -dot(n, l.getDirection());

    if (std::abs(f1) < 0.00001f)
        return false;
    else {
        f = f0 / f1;
        return true;
    }
}

template<class T>
bool Plane<T>::intersect(const Plane<T>& p, Line3<T>& l) const {
    Vector3<T> pn = p.n;
    T n00 = dot( n,  n);
    T n01 = dot( n, pn);
    T n11 = dot(pn, pn);
    T det = n00 * n11 - n01 * n01;

    if (std::abs(det) < 0.00001f)
        return false;
    else {
        T inv_det = 1.0f/det;
        T c0 = (-n11 *   -d + n01 * -p.d) * inv_det;
        T c1 = (-n00 * -p.d + n01 *   -d) * inv_det;
        l.setStart(n * c0 + pn * c1);
        l.setDirection(cross(n, pn));
        return true;
    }
}

template<class T>
void Plane<T>::clipPolygon(const std::vector< Vector3<T> >& vin, std::vector< Vector3<T> >& vout) const {
    // temporary array for distances
    std::vector<T> vd;
    vd.reserve(vin.size());

    for (size_t i = 0; i < vin.size(); ++i)
        vd.push_back( distance(vin[i]) );

    for (size_t i1 = 0, i2 = 1; i1 < vin.size(); ++i1, ++i2) {
        // wrap i2 around
        i2 = ( i2 == vin.size() ) ? 0 : i2;

        // inside-inside
        if (vd[i1] < T(0) && vd[i2] < T(0))
            vout.push_back( vin[i2] );
        // inside-outside OR outside-inside
        else if (vd[i1] * vd[i2] < 0.f) {
            T scale = vd[i1] / (vd[i1] - vd[i2]);
            // interpolate vertex
            Vector3<T> newv = vin[i1] + (vin[i2] - vin[i1]) * scale;
            vout.push_back(newv);

            // outside-inside
            if (vd[i1] > T(0) && vd[i2] < T(0))
                vout.push_back( vin[i2] );
        }
        // do nothing for outside-outside
    }
}

template<class T>
void Plane<T>::clipPolygon(const Vector3<T>* vin, size_t numIn, Vector3<T>* vout, size_t& numOut) const {
    // temporary array for distances
    T* vd = new T[numIn];
    numOut = 0;

    for (size_t i = 0; i < numIn; ++i)
        vd[i] = distance(vin[i]);

    for (size_t i1 = 0, i2 = 1; i1 < numIn; ++i1, ++i2) {
        // wrap i2 around
        i2 = ( i2 == numIn ) ? 0 : i2;

        // inside-inside
        if (vd[i1] < T(0) && vd[i2] < T(0)) {
            vout[numOut] = vin[i2];
            ++numOut;
        }
        // inside-outside OR outside-inside
        else if (vd[i1] * vd[i2] < 0.f) {
            T scale = vd[i1] / (vd[i1] - vd[i2]);
            // interpolate vertex
            vout[numOut] = vin[i1] + (vin[i2] - vin[i1]) * scale;
            ++numOut;

            // outside-inside
            if (vd[i1] > T(0) && vd[i2] < T(0)) {
                vout[numOut] = vin[i2];
                ++numOut;
            }
        }
        // do nothing for outside-outside
    }

    // clean up
    delete[] vd;
}

template<class T>
void Plane<T>::createCubeVertices(const Vector3<T>& llf, const Vector3<T>& urb, Vector3<T> cubeVertices[8]) {
    cubeVertices[0] = Vector3<T>(llf.x, llf.y, urb.z);// llb 0
    cubeVertices[1] = Vector3<T>(urb.x, llf.y, urb.z);// lrb 1
    cubeVertices[2] = Vector3<T>(urb.x, urb.y, urb.z);// urb 2
    cubeVertices[3] = Vector3<T>(llf.x, urb.y, urb.z);// ulb 3

    cubeVertices[4] = Vector3<T>(llf.x, llf.y, llf.z);// llf 4
    cubeVertices[5] = Vector3<T>(urb.x, llf.y, llf.z);// lrf 5
    cubeVertices[6] = Vector3<T>(urb.x, urb.y, llf.z);// urf 6
    cubeVertices[7] = Vector3<T>(llf.x, urb.y, llf.z);// ulf 7
}

template<class T>
void Plane<T>::createCubePlanes(const Vector3<T> cv[8], Plane<T> cubePlanes[6]) {
    cubePlanes[0] = Plane<T>( cv[2], cv[1], cv[0] );// front
    cubePlanes[1] = Plane<T>( cv[0], cv[3], cv[7] );// left
    cubePlanes[2] = Plane<T>( cv[4], cv[5], cv[6] );// back
    cubePlanes[3] = Plane<T>( cv[6], cv[2], cv[1] );// right
    cubePlanes[4] = Plane<T>( cv[5], cv[1], cv[0] );// bottom
    cubePlanes[5] = Plane<T>( cv[7], cv[3], cv[2] );// top
}

template<class T>
Plane<T> Plane<T>::transform(tgt::mat4 m) const {
    tgt::mat4 mInv;
    m.invert(mInv);
    tgt::mat4 mInvTr = transpose(mInv);
    vec3 newNorm = normalize((mInvTr * vec4(n, T(0))).xyz());
    float newDist = dot(newNorm, (m * vec4(d * n, T(1))).xyz());

    return Plane<T>(newNorm, newDist);
}

template<class T>
void Plane<T>::createQuad(Vector3<T> llf, Vector3<T> urb, Vector3<T>* v) const {
// copy of llf and urb is needed here
    T param;
    Line3<T> l;

    // scale to ensure that a polygon is created that is huge enough
    urb.x += T(5);
    urb.y += T(5);
    urb.z -= T(5);
    llf.x -= T(5);
    llf.y -= T(5);
    llf.z += T(5);

    // check whether the plane's normal lies in the x-y plane
    if (abs(n).z < T(0.000001)) {
        Vector3<T> llb(llf.x, llf.y, urb.z);
        Vector3<T> urf(urb.x, urb.y, llf.z);

        // each found vector must be scaled hence a overlapping polygon is guaranteed
        l = Line3<T>(urb, urb + Vector3<T>(T(1), T(0), T(0)));
        intersect(l, param);
        v[0] = l.getFromParam(param);

        l = Line3<T>(urf, urf + Vector3<T>(T(1), T(0), T(0)));
        intersect(l, param);
        v[1] = l.getFromParam(param);

        l = Line3<T>(llf, llf + Vector3<T>(T(1), T(0), T(0)));
        intersect(l, param);
        v[2] = l.getFromParam(param);

        l = Line3<T>(llb, llb + Vector3<T>(T(1), T(0), T(0)));
        intersect(l, param);
        v[3] = l.getFromParam(param);
    }
    else {
        Vector3<T> lrf(urb.x, llf.y, llf.z);
        Vector3<T> ulb(llf.x, urb.y, urb.z);

        // each found vector must be scaled hence a overlapping polygon is guaranteed
        l = Line3<T>(ulb, ulb + Vector3<T>(T(0), T(0), T(1)));
        intersect(l, param);
        v[0] = l.getFromParam(param);

        l = Line3<T>(urb, urb + Vector3<T>(T(0), T(0), T(1)));
        intersect(l, param);
        v[1] = l.getFromParam(param);

        l = Line3<T>(lrf, lrf + Vector3<T>(T(0), T(0), T(1)));
        intersect(l, param);
        v[2] = l.getFromParam(param);

        l = Line3<T>(llf, llf + Vector3<T>(T(0), T(0), T(1)));
        intersect(l, param);
        v[3] = l.getFromParam(param);
    }
}

template<class T>
void Plane<T>::clipAAB(Plane<T> aabPlanes[6], Vector3<T>* vin, size_t numIn, Vector3<T>* vout, size_t& numOut) {
    for (size_t i = 0; i < 6; ++i) {
        aabPlanes[i].clipPolygon(vin, numIn, vout, numOut);

        if (i == 5)
            return;

        // swap
        Vector3<T>* tmp = vin;
        vin = vout;
        vout = tmp;
        numIn = numOut;
    }
}

template<class T>
void Plane<T>::clipAAB(Plane<T> aabPlanes[6], std::vector< Vector3<T> >& vin, std::vector< Vector3<T> >& vout) {
    for (size_t i = 0; i < 6; ++i) {
        aabPlanes[i].clipPolygon(vin, vout);

        if (i == 5)
            return;

        std::swap(vin, vout);
        vout.clear();
    }
}

template<class T>
void Plane<T>::clipAAB(const Vector3<T>& llf, const Vector3<T>& urb, std::vector< Vector3<T> >& vout) const {
    // 1. compute AAB Planes
    Plane<T> aabPlanes[6];
    createAABPlanes(llf, urb, aabPlanes);

    // 2. calculate a polygon which is huge enough
    std::vector< Vector3<T> > vin(4);
    createQuad(llf, urb, &vin[0]);

    // 3. clip against all 6 planes of the AAB
    clipAAB(aabPlanes, vin, vout);
}

/// ostream-operator
template <class T>
std::ostream& operator << (std::ostream& s, const Plane<T>& p) {
    return (s << p.n << " * x + " << p.d << " = 0");
}

// typedefs for easy usage
// defined types with a small 'p' like other absolute math classes in tgt
typedef Plane<float>    plane;
typedef Plane<double>  dplane;

} // namespace tgt

#endif //TGT_PLANE_H
