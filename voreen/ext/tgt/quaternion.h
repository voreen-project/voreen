/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#ifndef TGT_QUATERNION_H
#define TGT_QUATERNION_H

#include <cmath>

#include "tgt/assert.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"

/**
    This is your standard Quaternion class.  It was originally
    implemented to control the mouse-view of a simple camera-class and
    therefore, at the moment, lacks many features that make Quaternions
    very useful for a lot of other things.
*/

namespace tgt {

template<class T>
struct Quaternion {
    typedef T ElemType;

    union {
        struct {
            T x;
            T y;
            T z;
            T w;
        };

        T elem[4];
    };

    /// Constructor, the standard is the identity-Quaternion.
    Quaternion (T X = 0, T Y = 0, T Z = 0, T W = 1)
                : x(X), y(Y), z(Z), w(W) {}

    /// Constructor expecting a 3d-vector and a scalar value
    Quaternion (Vector3<T> V, T W)
                : x(V.x), y(V.y), z(V.z), w(W) {}

    /// Init from array with size 4
    Quaternion (const T* t) {
        for (size_t i = 0; i < 4; ++i)
            elem[i] = t[i];
    }

    /// Init from matrix
    Quaternion (const Matrix4<T>& m) {
        *this = generateQuatFromMatrix(m);
    }

    /// calculate the absolute length of the Quaternion
    T abs() const {
        return std::sqrt(x * x + y * y +
                         z * z + w * w);
    }

    /// return the vector consisting of the imaginary parts of the Quaternion
    Vector3<T> vec() const {
        return Vector3<T>(x, y, z);
    }

    /// normalizes this Quaternion
    void normalize() {
        T ab = abs();

        if (ab != 0) {
            x /= ab;
            y /= ab;
            z /= ab;
            w /= ab;
        }
    }

    /// conjugates this Quaternion
    void conjugate() {
        x = -x;
        y = -y;
        z = -z;
    }

    /// inverts this Quaternion
    void invert() {
        T absq = x*x + y*y + z*z + w*w;
        if (absq != 0) {
            conjugate();
            *this = *this / absq;
        }
    }

    /// Calculates the Quaternion dot-product of this Quaternion
    /// with another Quaternion.  The *-operator is overloaded for this purpose.
    Quaternion operator * (const Quaternion& q) const {
        Quaternion res;

        res.x = w*q.x + x*q.w + y*q.z - z*q.y;
        res.y = w*q.y - x*q.z + y*q.w + z*q.x;
        res.z = w*q.z + x*q.y - y*q.x + z*q.w;
        res.w = w*q.w - x*q.x - y*q.y - z*q.z;

        return res;
    }

    /// Calculates the Quaternion dot-product of this Quaternion
    /// with another Quaternion.  The *-operator is overloaded for this purpose.
    Quaternion operator *= (const Quaternion& q) {
        *this = *this * q;
        return *this;
    }

    /// Quaternion addition
    Quaternion operator + (const Quaternion& q) const {
        Quaternion res;

        res.x = x + q.x;
        res.y = y + q.y;
        res.z = z + q.z;
        res.w = w + q.w;

        return res;
    }

    /// Quaternion subtraction
    Quaternion operator - (const Quaternion& q) const {
        Quaternion res;

        res.x = x - q.x;
        res.y = y - q.y;
        res.z = z - q.z;
        res.w = w - q.w;

        return res;
    }

    /// divide Quaternion by a float
    Quaternion operator / (const T f) const {
        Quaternion res;

        res.x = x / f;
        res.y = y / f;
        res.z = z / f;
        res.w = w / f;

        return res;
    }

    /// multiply Quaternion with a float
    Quaternion operator * (const T f) const {
        Quaternion res;

        res.x = x * f;
        res.y = y * f;
        res.z = z * f;
        res.w = w * f;

        return res;
    }

    /// multiply Quaternion with a vector
    Vector3<T> operator * (const Vector3<T> v) const {
        Vector3<T> tmp1, tmp2;
        Vector3<T> qvector = vec();

        tmp1 = cross(qvector, v);
        tmp2 = cross(qvector, tmp1);
        tmp1 = (T)(2 * w) * tmp1;
        tmp2 = (T)2 * tmp2;

        return v + tmp1 + tmp2;
    }

    /// equality operator
    bool operator == (const Quaternion& q) {
        return(x == q.x && y == q.y && z == q.z && w == q.w);
    }

    /// calculate the t-th power of this quaternion
    Quaternion pow(const T t) {
        using namespace std; // use overloaded cos and sin
        if (abs() != 1 || t < 0 || t > 0) return 0;

        double phi = acos(w);
        Vector3<T> u = vec() / sin(phi);

        double factor = sin(phi*t);
        Quaternion res = Quaternion(factor*u.x, factor*u.y, factor*u.z, cos(phi*t));
        return res;
    }

    /// calculate the natural logarithm of this quaternion
    Vector3<T> log() const {
        using namespace std; // use overloaded cos and sin
        if (abs() != 1) return 0;

        double phi = acos(w);
        Vector3<T> u = vec() / sin(phi);

        return u*phi;
    }

    /// create a Quaternion from a a given angle and axis representing that rotation
    static Quaternion createQuat(T angle, const Vector3<T>& axis) {
        Quaternion res;

        Vector3<T> newAxis = tgt::normalize(axis);

        res.x = newAxis.x * sinf(angle / 2);
        res.y = newAxis.y * sinf(angle / 2);
        res.z = newAxis.z * sinf(angle / 2);
        res.w = std::cos(angle / 2);

        return res;
    }

    /// rotate a vector about an axis using quaterions
    static Vector3<T> rotate(const Vector3<T>& target, T angle, const Vector3<T>& axis) {
        return rotate(target, createQuat(angle, axis));
    }

    /// rotate a vector about an axis using quaterions
    static Vector3<T> rotate(const Vector3<T>& target, Quaternion rot) {

        Quaternion rotconj = rot;
        rotconj.conjugate();

        // create a target quaternion so we can use quaternion-multiplication
        Quaternion targ(target.x, target.y, target.z, 0.f);

        // We calculate r * t * r^-1.  This yields our rotated vector.  Technically,
        // we would have to invert rot, but since it is a unit quaternion the inverted
        // quaternion is equal to the conjugated one, so we use conjugation which is cheaper.
        // Also note that the returned vector is automatically of length 1 since we only use
        // use unit-quaternions.
        Quaternion res = (rot * targ) * rotconj;

        return vec3(res.x, res.y, res.z);
    }
};

/**
    Multiply a float with a quaternion
*/
template<class T>
inline Quaternion<T> operator * (T f, const Quaternion<T>& q) {
    Quaternion<T> res;

    res.x = q.x * f;
    res.y = q.y * f;
    res.z = q.z * f;
    res.w = q.w * f;

    return res;
}

/**
    Calculates the Quaternion dot-product of two Quaternions.
*/
template<class T>
inline Quaternion<T> multQuat(const Quaternion<T>& q1, const Quaternion<T>& q2) {
    Quaternion<T> res;

    res.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
    res.y = q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x;
    res.z = q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w;
    res.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;

    return res;
}

/**
    e^quaternion given as:
    Qexp(v*a) = [cos(a),vsin(a)]
*/
template<class T>
Quaternion<T> expQuat(const Quaternion<T>& q) {
    float a = std::sqrt(q.x*q.x + q.y*q.y + q.y*q.y);
    float sina = sinf(a);
    Quaternion<T> ret;

    ret.w = cosf(a);

    if (a > 0) {
        ret.x = sina * q.x / a;
        ret.y = sina * q.y / a;
        ret.z = sina * q.z / a;
    }
    else {
        ret.x = q.x;
        ret.y = q.y;
        ret.z = q.z; //= ret.y = ret.z = 0;
    }

    return ret;
}

/**
    Returns the logarithm of a Quaternion.
*/
template<class T>
Quaternion<T> logQuat(const Quaternion<T>& q) {
    float a = acosf(q.w);
    float sina = sinf(a);
    Quaternion<T> ret;

    ret.w = 0;

    if (sina > 0) {
            ret.x = a*q.x/sina;
            ret.y = a*q.y/sina;
            ret.z = a*q.z/sina;
    }
    else {
            ret.x=ret.y=ret.z=0;
    }

    return ret;
}

/**
     Returns a Quaternion in normalized form.
*/
template<class T>
inline Quaternion<T> normalize(const Quaternion<T>& q) {
    T ab = q.abs();
    Quaternion<T> res = q;

    if (ab != 0) {
        res.x /= ab;
        res.y /= ab;
        res.z /= ab;
        res.w /= ab;
    }

    return res;
}

/**
    Returns a Quaternion in conjugated form.
*/
template<class T>
inline Quaternion<T> conjugate(const Quaternion<T>& q) {
    Quaternion<T> res = q;

    res.x *= -1;
    res.y *= -1;
    res.z *= -1;

    return res;
}

/**
    Return the inverted form of a Quaternion.
*/
template<class T>
inline Quaternion<T> invert(const Quaternion<T>& q) {
    T absq = q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
    if (absq != 0) {
        Quaternion<T> res = conjugate(q);
        return res / absq;
    }

    return q;
}

/**
    Lerp is the simplest form of Quaternion interpolation; the two given Quaternions are simply linearly
    interpolated.  Lerp may not yield as nice results as slerp, but it is faster and in certain situations
    sufficient.
*/
template<class T>
inline Quaternion<T> lerpQuat(const Quaternion<T>& a, const Quaternion<T>& b, T t) {
    tgtAssert(!((t < 0.0) || (t > 1.0)), "running-parameter must be between 0 and 1");

    return normalize((b - a)*t + a);
}

/**
    Slerp (Spherical Linear Interpolation) returns a interpolated Quaternion between the two given
    Quaternions, depending on the interpolation variable t.  This is based on on the chapter concerning
    Quaternions in Real Time Rendering and subsequently on Ken Shoemakes paper "Animating rotation with
    quaternion curves".
*/
template<class T>
Quaternion<T> slerpQuat(const Quaternion<T>& a, const Quaternion<T>& b, T t, bool inv = true) {
    using namespace std; // use overloaded cos and sin
    tgtAssert(!((t < T(0.0)) || (t > T(1.0))), "running-parameter must be between 0 and 1");

    Quaternion<T> c = b;
    T tmp = a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z;

    if (inv && tmp < T(0.0)) {
        tmp = -tmp;
        c = -1.f*c;
    }

    if (tmp < T(1.0) && tmp > T(-1.0)) {
        T ang = acos(tmp);
        return normalize(((a * sin(ang * (T(1.0) - t))) + (c * sin(ang * t))) / sin(ang));
    }
    else {
        return lerpQuat(a, c, t);
    }
}

/**
    This function returns an interpolated Quaternion between p and q, using a and b as "Control-
    Quaternions" and t as the interpolation variable.  This is based on the chapter concerning
    Quaternions in Real Time Rendering.
*/
template<class T>
inline Quaternion<T> squadQuat(const Quaternion<T>& p, const Quaternion<T>& q,
                               const Quaternion<T>& a, const Quaternion<T>& b, double t)
{
    tgtAssert(!((t < 0.0) || (t > 1.0)), "running-parameter must be between 0 and 1");

    return slerpQuat<T>(slerpQuat<T>(p, q, static_cast<T>(t)), slerpQuat<T>(a, b, static_cast<T>(t)), static_cast<T>(2.0*t*(1.0 - t)), false);
}

/**
    Generate "Control-Quaternions" that can be used for the squad-Algorithm to interpolate between
    given Quaternions.  This is based on the chapter concerning Quaternions in Real Time Rendering.
*/
template<class T>
inline Quaternion<T> splineQuat(const Quaternion<T>& pa, const Quaternion<T>& p, const Quaternion<T>& pb) {
    Quaternion<T> pinvert = conjugate(p);

    return p * expQuat((logQuat(pa*pinvert) + logQuat(pb*pinvert)) * static_cast<T>(-0.25));
}

/**
    Generate a Quaternion that represents the rotation from the source to the destination
    vector.
*/
template<class T>
Quaternion<T> generateQuaternionFromTo(const Vector3<T> src, const Vector3<T> dest) {
    Quaternion<T> q;

    Vector3<T> v1 = normalize(src);
    Vector3<T> v2 = normalize(dest);

    Vector3<T> crs = cross(v1, v2);

    T v1v2dot = dot(v1, v2);

    if (v1v2dot == 1) {// the vectors are identical
        return Quaternion<T>(0, 0, 0, 1); // ... so we return a rotation that does nothing
    }

    if (v1v2dot+1 < T(1e-11)) {
        // Pick a reasonable rotation axis if vectors face in (almost) the opposite direction
        // for better numerical stability
        Vector3<T> axis = cross(Vector3<T>(1, 0, 0), src);
        if (length(axis) == 0) {
            axis = cross(Vector3<T>(0, 1, 0), src);
        }
        crs = normalize(axis);
    }
    // Let a be the angle between v1 and v2
    // crs' length is 1*sin(a)
    q.x = crs.x;
    q.y = crs.y;
    q.z = crs.z;
    // v1v2dot is cos(a)
    q.w = T(1) + v1v2dot;
    // After normalization q.w will be cos(a/2) and q.xyz's length will be sin(a/2),
    // which is exactly what we want! Do the math!

    return normalize(q);
}

/**
    Generates a rotational matrix from a unit quaternion
*/
template<class T>
Matrix4<T> generateMatrixFromQuat(const Quaternion<T>& q) {
    T xx      = q.x * q.x;
    T xy      = q.x * q.y;
    T xz      = q.x * q.z;
    T xw      = q.x * q.w;
    T yy      = q.y * q.y;
    T yz      = q.y * q.z;
    T yw      = q.y * q.w;
    T zz      = q.z * q.z;
    T zw      = q.z * q.w;

    Matrix4<T> mat;

    mat.elem[0]  = 1 - 2 * ( yy + zz );
    mat.elem[1]  =     2 * ( xy - zw );
    mat.elem[2]  =     2 * ( xz + yw );
    mat.elem[4]  =     2 * ( xy + zw );
    mat.elem[5]  = 1 - 2 * ( xx + zz );
    mat.elem[6]  =     2 * ( yz - xw );
    mat.elem[8]  =     2 * ( xz - yw );
    mat.elem[9]  =     2 * ( yz + xw );
    mat.elem[10] = 1 - 2 * ( xx + yy );
    mat.elem[3]  = mat.elem[7] = mat.elem[11] = mat.elem[12] = mat.elem[13] = mat.elem[14] = 0;
    mat.elem[15] = 1;

    return transpose(mat);
}
/**
    Generates a unit quaternion from a rotational matrix
*/
template<class T>
Quaternion<T> generateQuatFromMatrix(const Matrix4<T>& mat) {
    Quaternion<T> q;
    T t = mat[0][0] + mat[1][1] + mat[2][2];

    if (t > T(0)) {
        T s = T(2.0) * std::sqrt(t + T(1));
        q.w = T(0.25) * s;
        q.x = ( mat[1][2] - mat[2][1] ) / s;
        q.y = ( mat[2][0] - mat[0][2] ) / s;
        q.z = ( mat[0][1] - mat[1][0] ) / s;
    } else if ((mat[0][0] > mat[1][1]) && (mat[0][0] > mat[2][2])) {
        T s = std::sqrt( T(1.0) + mat[0][0] - mat[1][1] - mat[2][2] ) * T(2.0);
        q.w = (mat[1][2] - mat[2][1] ) / s;
        q.x = T(0.25) * s;
        q.y = (mat[1][0] + mat[0][1] ) / s;
        q.z = (mat[2][0] + mat[0][2] ) / s;
    } else if (mat[1][1] > mat[2][2]) {
        T s = std::sqrt( T(1.0) + mat[1][1] - mat[0][0] - mat[2][2] ) * T(2.0);
        q.w = (mat[2][0] - mat[0][2] ) / s;
        q.x = (mat[1][0] + mat[0][1] ) / s;
        q.y = T(0.25) * s;
        q.z = (mat[2][1] + mat[1][2] ) / s;
    } else {
        T s = std::sqrt( T(1.0) + mat[2][2] - mat[0][0] - mat[1][1] ) * T(2.0);
        q.w = (mat[0][1] - mat[1][0] ) / s;
        q.x = (mat[0][2] + mat[2][0] ) / s;
        q.y = (mat[1][2] + mat[2][1] ) / s;
        q.z = T(0.25) * s;
    }

    return q;
}

template<class T>
void generateAxisAngleFromQuat(const Quaternion<T>& qu, T& angle, Vector3<T>& axis) {

    Quaternion<T> q = qu;

    if (q.w > 1)
        q.normalize();

    angle = 2 * std::acos(q.w);
    T s = std::sqrt(1 - q.w * q.w);

    if (s < 0.001) {
        axis.x = q.x;
        axis.y = q.y;
        axis.z = q.z;
    } else {
        axis.x = q.x / s; // normalise axis
        axis.y = q.y / s;
        axis.z = q.z / s;
    }
}

/**
    ostream-operator for output
*/
template<typename T>
std::ostream& operator << (std::ostream& s, const Quaternion<T>& q) {
    return (s << "{ " << q.x << " " << q.y << " " << q.z << " " << q.w << " }" << std::endl);
}

// just a typedef for easier use, since most of the time float-Quaternions are used
typedef Quaternion<float>   quat;
typedef Quaternion<double>  dquat;

} // namespace tgt

#endif //_QUATERNION_H
